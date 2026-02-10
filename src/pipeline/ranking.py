"""Candidate ranking using worst-metric-rank and diversity selection.

Implements the BoltzGen-validated ranking approach:
1. Worst-metric-rank: For each candidate, rank across all metrics independently,
   then use the worst (highest) weighted rank as the quality key.
2. Greedy maximin diversity selection: Iteratively pick candidates that maximize
   a combination of quality and sequence diversity.

This replaces the broken composite score which put 30% weight on pDockQ
(always 0.0 from Boltz-2) and used interface metrics only as binary filters.
"""

from dataclasses import dataclass, field
from typing import Optional


# Default metric weights following BoltzGen pattern.
# Lower weight = metric rank is divided by a larger number = less penalized.
# Weight 1 = high importance, weight 2 = lower importance.
DEFAULT_METRIC_WEIGHTS = {
    "iptm": 1,              # Best available interface quality signal
    "ptm": 1,               # Fold confidence
    "interface_area": 1,    # Binding interface size
    "humanness": 1,         # Critical for therapeutics
    "num_contacts": 2,      # Correlated with area, lower importance
    "plddt": 2,             # Typically high, less discriminating
    "proteinmpnn_ll": 1,   # Inverse folding NLL (lower = better)
    "antifold_ll": 2,      # Antibody-specific inverse folding NLL (lower = better)
    "protenix_iptm": 1,    # Cross-validation interface quality
}


@dataclass
class RankedCandidate:
    """A candidate with per-metric ranks and final ranking."""

    candidate_id: str
    sequence: str
    sequence_vl: Optional[str] = None

    # Raw metrics
    iptm: float = 0.0
    ptm: float = 0.0
    plddt: float = 0.0
    interface_area: float = 0.0
    num_contacts: int = 0
    humanness: float = 0.0

    # Validation metrics (from step 04a)
    proteinmpnn_ll: Optional[float] = None
    antifold_ll: Optional[float] = None
    protenix_iptm: Optional[float] = None

    # Per-metric ranks (1 = best)
    metric_ranks: dict[str, int] = field(default_factory=dict)

    # Weighted ranks (rank / weight)
    weighted_ranks: dict[str, float] = field(default_factory=dict)

    # quality_key = max of weighted ranks (lower = better)
    quality_key: float = float("inf")

    # Final rank after sorting by quality_key
    final_rank: int = 0

    # Reference back to CandidateScore (for output mapping)
    _score_ref: Optional[object] = field(default=None, repr=False)

    def get_binding_sequence(self) -> str:
        """Get the full binding sequence for identity comparison."""
        if self.sequence_vl:
            return self.sequence + self.sequence_vl
        return self.sequence


def worst_metric_rank(
    candidates: list[RankedCandidate],
    metric_weights: Optional[dict[str, int]] = None,
) -> list[RankedCandidate]:
    """Rank candidates using worst-metric-rank approach.

    For each metric, rank all candidates (best=1). Divide each rank by the
    metric's importance weight. The quality_key is the maximum weighted rank
    across all metrics. Sort ascending (lower quality_key = better).

    Args:
        candidates: List of RankedCandidate objects with raw metrics populated.
        metric_weights: Dict mapping metric name to importance weight.
            Higher weight = more important = rank divided by smaller number.

    Returns:
        Sorted list of candidates (best first) with ranks populated.
    """
    if not candidates:
        return []

    weights = metric_weights or DEFAULT_METRIC_WEIGHTS

    # Define metric accessors.
    # Validation metrics may be None — candidates without a score get worst rank.
    metric_accessors = {
        "iptm": lambda c: c.iptm,
        "ptm": lambda c: c.ptm,
        "plddt": lambda c: c.plddt,
        "interface_area": lambda c: c.interface_area,
        "num_contacts": lambda c: c.num_contacts,
        "humanness": lambda c: c.humanness,
        "proteinmpnn_ll": lambda c: c.proteinmpnn_ll,
        "antifold_ll": lambda c: c.antifold_ll,
        "protenix_iptm": lambda c: c.protenix_iptm,
    }

    # Metrics where lower values are better (NLL scores)
    lower_is_better = {"proteinmpnn_ll", "antifold_ll"}

    n = len(candidates)

    # Rank candidates per metric
    for metric_name, accessor in metric_accessors.items():
        if metric_name not in weights:
            continue

        weight = weights[metric_name]

        # Check if any candidates have this metric
        has_metric = [c for c in candidates if accessor(c) is not None]
        if not has_metric:
            # No candidates have this metric — skip entirely
            continue

        # Sort: best = rank 1. Candidates with None get worst rank.
        # For most metrics: higher = better (sort descending)
        # For NLL metrics: lower = better (sort ascending)
        descending = metric_name not in lower_is_better
        sorted_by_metric = sorted(
            candidates,
            key=lambda c, acc=accessor: (0 if acc(c) is None else 1, acc(c) if acc(c) is not None else 0),
            reverse=descending,
        )

        for rank_idx, candidate in enumerate(sorted_by_metric):
            rank = rank_idx + 1
            candidate.metric_ranks[metric_name] = rank
            candidate.weighted_ranks[metric_name] = rank / weight

    # Compute quality_key = max weighted rank
    for candidate in candidates:
        if candidate.weighted_ranks:
            candidate.quality_key = max(candidate.weighted_ranks.values())
        else:
            candidate.quality_key = float("inf")

    # Sort by quality_key ascending (lower = better), tiebreak by ipTM descending
    candidates.sort(key=lambda c: (c.quality_key, -c.iptm))

    # Assign final ranks
    for i, candidate in enumerate(candidates):
        candidate.final_rank = i + 1

    return candidates


def pairwise_sequence_identity(seq1: str, seq2: str) -> float:
    """Compute pairwise sequence identity between two sequences.

    Simple positional comparison, normalized by the length of the longer sequence.

    Args:
        seq1: First amino acid sequence.
        seq2: Second amino acid sequence.

    Returns:
        Fraction identity (0.0 to 1.0).
    """
    if not seq1 or not seq2:
        return 0.0

    max_len = max(len(seq1), len(seq2))
    min_len = min(len(seq1), len(seq2))

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)

    return matches / max_len if max_len > 0 else 0.0


def diversity_select(
    candidates: list[RankedCandidate],
    n_select: int,
    alpha: float = 0.001,
) -> list[RankedCandidate]:
    """Greedy maximin diversity selection.

    Pick highest quality first, then iteratively pick the candidate maximizing:
        (1 - alpha) * quality_score + alpha * (1 - max_identity_to_selected)

    Quality score is normalized so rank 1 = 1.0, last rank = 0.0.

    Args:
        candidates: Pre-ranked list of RankedCandidate (sorted by quality_key).
        n_select: Number of candidates to select.
        alpha: Diversity weight (0 = quality only, 1 = diversity only).
            Default 0.001 matches BoltzGen's approach.

    Returns:
        Selected candidates in selection order.
    """
    if not candidates:
        return []
    if n_select >= len(candidates):
        return list(candidates)

    n = len(candidates)

    # Normalize quality: best (rank 1) = 1.0, worst (rank n) = 0.0
    quality_scores = {}
    for c in candidates:
        if n > 1:
            quality_scores[c.candidate_id] = 1.0 - (c.final_rank - 1) / (n - 1)
        else:
            quality_scores[c.candidate_id] = 1.0

    # Pre-compute binding sequences for identity comparison
    binding_seqs = {c.candidate_id: c.get_binding_sequence() for c in candidates}

    selected = []
    remaining = list(candidates)

    # Pick best quality first
    selected.append(remaining.pop(0))

    while len(selected) < n_select and remaining:
        best_score = -1.0
        best_idx = 0

        for idx, candidate in enumerate(remaining):
            # Max identity to any already-selected candidate
            max_identity = max(
                pairwise_sequence_identity(
                    binding_seqs[candidate.candidate_id],
                    binding_seqs[s.candidate_id],
                )
                for s in selected
            )

            diversity_score = 1.0 - max_identity
            quality = quality_scores[candidate.candidate_id]
            combined = (1.0 - alpha) * quality + alpha * diversity_score

            if combined > best_score:
                best_score = combined
                best_idx = idx

        selected.append(remaining.pop(best_idx))

    return selected
