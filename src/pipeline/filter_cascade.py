"""Multi-stage filtering cascade for candidate selection.

This module implements the filtering pipeline that takes
raw design candidates and filters them through multiple
quality gates to select final candidates.
"""

from dataclasses import dataclass, field
from typing import Optional, Any
from enum import Enum


class FilterResult(Enum):
    """Result of a filter check."""

    PASS = "pass"
    FAIL = "fail"
    SOFT_FAIL = "soft_fail"  # Flagged but not rejected


@dataclass
class CandidateScore:
    """Scores and metrics for a design candidate."""

    candidate_id: str
    sequence: str
    sequence_vl: Optional[str] = None
    binder_type: str = "vhh"  # "vhh" or "scfv"
    source: str = "unknown"  # "denovo" or "optimized"

    # Binding metrics
    iptm: Optional[float] = None
    pdockq: Optional[float] = None
    interface_area: Optional[float] = None
    num_contacts: Optional[int] = None

    # Epitope annotation
    epitope_class: str = "unknown"  # "OKT3-like" or "novel_epitope"
    okt3_overlap: float = 0.0

    # Humanness
    oasis_score_vh: Optional[float] = None
    oasis_score_vl: Optional[float] = None
    oasis_score_mean: Optional[float] = None
    has_back_mutation_variant: bool = False

    # Liabilities (positions)
    deamidation_sites: list[int] = field(default_factory=list)
    isomerization_sites: list[int] = field(default_factory=list)
    glycosylation_sites: list[int] = field(default_factory=list)
    oxidation_sites: list[int] = field(default_factory=list)
    unpaired_cys: int = 0

    # CDR-specific liability counts (for targeted filtering)
    cdr_deamidation_count: int = 0
    cdr_isomerization_count: int = 0
    cdr_glycosylation_count: int = 0
    cdr_oxidation_count: int = 0

    # Developability
    cdr_h3_length: Optional[int] = None
    net_charge: Optional[float] = None
    isoelectric_point: Optional[float] = None
    hydrophobic_patches: int = 0
    cdr_positions: Optional[dict[str, tuple[int, int]]] = None
    full_sequence: Optional[str] = None

    # Filter results
    filter_results: dict[str, FilterResult] = field(default_factory=dict)
    risk_flags: list[str] = field(default_factory=list)

    # Composite score (computed)
    composite_score: float = 0.0
    rank: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "candidate_id": self.candidate_id,
            "sequence": self.sequence,
            "sequence_vl": self.sequence_vl,
            "binder_type": self.binder_type,
            "source": self.source,
            "binding": {
                "iptm": self.iptm,
                "pdockq": self.pdockq,
                "pdockq_note": "Structural confidence, NOT affinity predictor",
                "interface_area": self.interface_area,
                "num_contacts": self.num_contacts,
            },
            "epitope": {
                "epitope_class": self.epitope_class,
                "okt3_overlap": self.okt3_overlap,
            },
            "humanness": {
                "oasis_score_vh": self.oasis_score_vh,
                "oasis_score_vl": self.oasis_score_vl,
                "oasis_score_mean": self.oasis_score_mean,
                "has_back_mutation_variant": self.has_back_mutation_variant,
            },
            "liabilities": {
                "deamidation_sites": self.deamidation_sites,
                "isomerization_sites": self.isomerization_sites,
                "glycosylation_sites": self.glycosylation_sites,
                "oxidation_sites": self.oxidation_sites,
                "unpaired_cys": self.unpaired_cys,
                "cdr_deamidation_count": self.cdr_deamidation_count,
                "cdr_isomerization_count": self.cdr_isomerization_count,
                "cdr_glycosylation_count": self.cdr_glycosylation_count,
                "cdr_oxidation_count": self.cdr_oxidation_count,
            },
            "developability": {
                "cdr_h3_length": self.cdr_h3_length,
                "net_charge": self.net_charge,
                "isoelectric_point": self.isoelectric_point,
                "hydrophobic_patches": self.hydrophobic_patches,
            },
            "filter_results": {k: v.value for k, v in self.filter_results.items()},
            "risk_flags": self.risk_flags,
            "composite_score": self.composite_score,
            "rank": self.rank,
        }


class FilterCascade:
    """Multi-stage filtering cascade for candidate selection.

    Filters are applied in order:
    1. Binding quality (calibrated thresholds)
    2. Humanness
    3. Sequence liabilities
    4. Developability properties
    5. Aggregation propensity
    """

    def __init__(self, config: Optional[Any] = None):
        """Initialize filter cascade.

        Args:
            config: PipelineConfig or FilteringConfig.
        """
        self.config = config
        self._thresholds = self._get_thresholds()

    def _get_thresholds(self) -> dict:
        """Get filter thresholds from config or defaults."""
        if self.config is None:
            return {
                "min_pdockq": 0.5,
                "min_interface_area": 800.0,
                "min_contacts": 10,
                "min_oasis_score": 0.8,
                "cdr_h3_length_range": (8, 20),
                "net_charge_range": (-2, 4),
                "pi_range": (6.0, 9.0),
                "max_hydrophobic_patches": 2,
                "max_oxidation_sites": 2,
            }

        # Handle both PipelineConfig and FilteringConfig
        if hasattr(self.config, "get_effective_thresholds"):
            effective = self.config.get_effective_thresholds()
            filtering = self.config.filtering
        else:
            effective = {
                "min_pdockq": self.config.min_pdockq,
                "min_interface_area": self.config.min_interface_area,
                "min_contacts": self.config.min_contacts,
            }
            filtering = self.config

        return {
            "min_pdockq": effective.get("min_pdockq", 0.5),
            "min_interface_area": effective.get("min_interface_area", 800.0),
            "min_contacts": effective.get("min_contacts", 10),
            "min_oasis_score": getattr(filtering, "min_oasis_score", 0.8),
            "cdr_h3_length_range": getattr(filtering, "cdr_h3_length_range", (8, 20)),
            "net_charge_range": getattr(filtering, "net_charge_range", (-2, 4)),
            "pi_range": getattr(filtering, "pi_range", (6.0, 9.0)),
            "max_hydrophobic_patches": getattr(filtering, "max_hydrophobic_patches", 2),
            "max_oxidation_sites": getattr(filtering, "max_oxidation_sites", 2),
            "allow_deamidation_cdr": getattr(filtering, "allow_deamidation_cdr", False),
            "allow_isomerization_cdr": getattr(filtering, "allow_isomerization_cdr", False),
            "allow_glycosylation_cdr": getattr(filtering, "allow_glycosylation_cdr", False),
        }

    def filter_binding(self, candidate: CandidateScore) -> FilterResult:
        """Filter by binding quality metrics.

        Interface area is the primary hard filter. pDockQ is only checked
        if the threshold is non-zero AND the candidate has a non-zero value
        (Boltz-2 does not produce pDockQ — it's always 0.0).

        Soft fails if interface_area or num_contacts are missing,
        as this indicates incomplete binding evidence.
        """
        # Track if we have incomplete binding data
        has_incomplete_data = False

        # Interface area — primary hard filter
        if candidate.interface_area is not None:
            if candidate.interface_area < self._thresholds["min_interface_area"]:
                return FilterResult.FAIL
        else:
            has_incomplete_data = True

        # Contact count — hard filter
        if candidate.num_contacts is not None:
            if candidate.num_contacts < self._thresholds["min_contacts"]:
                return FilterResult.FAIL
        else:
            has_incomplete_data = True

        # pDockQ — only apply if threshold is non-zero AND candidate has non-zero value.
        # Boltz-2 does not produce pDockQ (always 0.0), so this avoids rejecting
        # every candidate when calibrated_min_pdockq is also 0.0.
        min_pdockq = self._thresholds["min_pdockq"]
        if min_pdockq > 0 and candidate.pdockq is not None and candidate.pdockq > 0:
            if candidate.pdockq < min_pdockq:
                return FilterResult.FAIL

        # Soft-fail if binding evidence is incomplete
        if has_incomplete_data:
            return FilterResult.SOFT_FAIL

        return FilterResult.PASS

    def filter_humanness(self, candidate: CandidateScore) -> FilterResult:
        """Filter by humanness score."""
        score = candidate.oasis_score_mean
        if score is None:
            score = candidate.oasis_score_vh

        if score is None:
            return FilterResult.SOFT_FAIL  # Can't assess, flag but don't reject

        if score < self._thresholds["min_oasis_score"]:
            return FilterResult.FAIL

        return FilterResult.PASS

    def filter_liabilities(self, candidate: CandidateScore) -> FilterResult:
        """Filter by sequence liabilities.

        When allow_*_cdr is False, only CDR liabilities are rejected.
        Framework region liabilities are allowed as they're less likely to affect binding.
        """
        # Hard filters - check CDR-specific counts when configured
        if not self._thresholds["allow_deamidation_cdr"] and candidate.cdr_deamidation_count > 0:
            return FilterResult.FAIL

        if not self._thresholds["allow_isomerization_cdr"] and candidate.cdr_isomerization_count > 0:
            return FilterResult.FAIL

        if not self._thresholds["allow_glycosylation_cdr"] and candidate.cdr_glycosylation_count > 0:
            return FilterResult.FAIL

        if candidate.unpaired_cys > 0:
            return FilterResult.FAIL

        # Soft filter: oxidation (check total, as even framework oxidation matters for stability)
        if len(candidate.oxidation_sites) > self._thresholds["max_oxidation_sites"]:
            return FilterResult.SOFT_FAIL

        return FilterResult.PASS

    def filter_developability(self, candidate: CandidateScore) -> FilterResult:
        """Filter by developability properties."""
        # CDR-H3 length
        if candidate.cdr_h3_length is not None:
            min_len, max_len = self._thresholds["cdr_h3_length_range"]
            if not (min_len <= candidate.cdr_h3_length <= max_len):
                return FilterResult.SOFT_FAIL

        # Net charge
        if candidate.net_charge is not None:
            min_charge, max_charge = self._thresholds["net_charge_range"]
            if not (min_charge <= candidate.net_charge <= max_charge):
                return FilterResult.SOFT_FAIL

        # Isoelectric point
        if candidate.isoelectric_point is not None:
            min_pi, max_pi = self._thresholds["pi_range"]
            if not (min_pi <= candidate.isoelectric_point <= max_pi):
                return FilterResult.SOFT_FAIL

        # Hydrophobic patches
        if candidate.hydrophobic_patches > self._thresholds["max_hydrophobic_patches"]:
            return FilterResult.SOFT_FAIL

        return FilterResult.PASS

    def filter_aggregation(self, candidate: CandidateScore) -> FilterResult:
        """Filter by aggregation propensity indicators.

        Checks for:
        1. High aromatic content in CDRs (>20% suggests aggregation risk)
        2. Consecutive aromatic residues in CDRs (2+ in a row is problematic)
        """
        from src.utils.constants import AROMATIC

        sequence = candidate.full_sequence or (candidate.sequence or "")
        if not candidate.full_sequence and candidate.sequence_vl:
            sequence += candidate.sequence_vl

        if not sequence:
            return FilterResult.SOFT_FAIL

        # Check aromatic content in CDRs when positions are available.
        cdr_positions = candidate.cdr_positions or {}
        if cdr_positions:
            cdr_sequences = [
                sequence[start:end + 1]
                for start, end in cdr_positions.values()
                if start < len(sequence)
            ]
            cdr_sequence = "".join(cdr_sequences)
            if cdr_sequence:
                aromatic_count = sum(1 for aa in cdr_sequence if aa in AROMATIC)
                aromatic_fraction = aromatic_count / len(cdr_sequence)
                if aromatic_fraction > 0.20:  # >20% aromatic in CDRs
                    return FilterResult.SOFT_FAIL

                # Check for consecutive aromatics within each CDR
                for cdr_seq in cdr_sequences:
                    consecutive = 0
                    max_consecutive = 0
                    for aa in cdr_seq:
                        if aa in AROMATIC:
                            consecutive += 1
                            max_consecutive = max(max_consecutive, consecutive)
                        else:
                            consecutive = 0
                    if max_consecutive >= 2:
                        return FilterResult.SOFT_FAIL
        else:
            # Fallback: use full sequence if CDR positions are unavailable.
            aromatic_count = sum(1 for aa in sequence if aa in AROMATIC)
            aromatic_fraction = aromatic_count / len(sequence) if sequence else 0
            if aromatic_fraction > 0.15:
                return FilterResult.SOFT_FAIL

            consecutive = 0
            max_consecutive = 0
            for aa in sequence:
                if aa in AROMATIC:
                    consecutive += 1
                    max_consecutive = max(max_consecutive, consecutive)
                else:
                    consecutive = 0
            if max_consecutive >= 3:
                return FilterResult.SOFT_FAIL

        return FilterResult.PASS

    def run_all_filters(self, candidate: CandidateScore) -> CandidateScore:
        """Run all filters on a candidate.

        Args:
            candidate: Candidate to filter.

        Returns:
            Updated candidate with filter results.
        """
        # Run each filter
        candidate.filter_results["binding"] = self.filter_binding(candidate)
        candidate.filter_results["humanness"] = self.filter_humanness(candidate)
        candidate.filter_results["liabilities"] = self.filter_liabilities(candidate)
        candidate.filter_results["developability"] = self.filter_developability(candidate)
        candidate.filter_results["aggregation"] = self.filter_aggregation(candidate)

        # Add risk flags
        for filter_name, result in candidate.filter_results.items():
            if result == FilterResult.SOFT_FAIL:
                candidate.risk_flags.append(f"{filter_name}_soft_fail")

        return candidate

    def passes_hard_filters(self, candidate: CandidateScore) -> bool:
        """Check if candidate passes all hard filters."""
        for result in candidate.filter_results.values():
            if result == FilterResult.FAIL:
                return False
        return True

    def compute_composite_score(
        self,
        candidate: CandidateScore,
        weights: Optional[dict] = None,
    ) -> float:
        """Compute composite score for ranking.

        Args:
            candidate: Candidate to score.
            weights: Optional custom weights.

        Returns:
            Composite score (0-1).
        """
        if weights is None:
            weights = {
                # NOTE: "structural_confidence" uses pDockQ, which measures confidence
                # in the PREDICTED STRUCTURE, NOT binding affinity. High pDockQ means
                # the model is confident the complex folds correctly, not that binding
                # is strong. This is appropriate for filtering unreliable predictions.
                "structural_confidence": 0.30,
                "humanness": 0.25,
                "liabilities": 0.25,
                "developability": 0.20,
            }

        score = 0.0

        # Structural confidence score (normalized pDockQ, 0-1 range)
        # NOTE: pDockQ is structural confidence, NOT affinity prediction.
        # We use it to filter low-confidence predictions, not to rank by binding.
        if candidate.pdockq is not None:
            score += weights["structural_confidence"] * min(candidate.pdockq, 1.0)

        # Humanness score
        humanness = candidate.oasis_score_mean
        if humanness is None:
            humanness = candidate.oasis_score_vh
        if humanness is None:
            humanness = 0.5
        score += weights["humanness"] * humanness

        # Liability score (inverse of liability count)
        total_liabilities = (
            len(candidate.deamidation_sites) +
            len(candidate.isomerization_sites) +
            len(candidate.glycosylation_sites) +
            len(candidate.oxidation_sites)
        )
        liability_score = max(0, 1 - total_liabilities * 0.1)
        score += weights["liabilities"] * liability_score

        # Developability score (based on flags)
        dev_flags = sum(1 for r in candidate.filter_results.values() if r == FilterResult.SOFT_FAIL)
        dev_score = max(0, 1 - dev_flags * 0.2)
        score += weights["developability"] * dev_score

        candidate.composite_score = score
        return score

    def filter_candidates(
        self,
        candidates: list[CandidateScore],
        compute_scores: bool = True,
    ) -> tuple[list[CandidateScore], list[CandidateScore]]:
        """Filter a list of candidates.

        Args:
            candidates: List of candidates to filter.
            compute_scores: If True, compute composite scores.

        Returns:
            Tuple of (passing_candidates, failing_candidates).
        """
        passing = []
        failing = []

        for candidate in candidates:
            self.run_all_filters(candidate)

            if compute_scores:
                self.compute_composite_score(candidate)

            if self.passes_hard_filters(candidate):
                passing.append(candidate)
            else:
                failing.append(candidate)

        # Sort passing by composite score
        passing.sort(key=lambda c: c.composite_score, reverse=True)

        # Assign ranks
        for i, candidate in enumerate(passing):
            candidate.rank = i + 1

        return passing, failing


def run_filter_cascade(
    candidates: list[CandidateScore],
    config: Optional[Any] = None,
    min_candidates: int = 10,
) -> tuple[list[CandidateScore], dict]:
    """Run filter cascade with fallback for insufficient candidates.

    Args:
        candidates: List of candidates to filter.
        config: Pipeline or filtering configuration.
        min_candidates: Minimum number of candidates required.

    Returns:
        Tuple of (selected_candidates, filter_stats).
    """
    cascade = FilterCascade(config)

    # Get fallback config options
    relax_soft_first = True
    max_relaxation = 0.1
    if config is not None:
        if hasattr(config, "filtering"):
            relax_soft_first = getattr(config.filtering, "relax_soft_filters_first", True)
            max_relaxation = getattr(config.filtering, "max_threshold_relaxation", 0.1)
        else:
            relax_soft_first = getattr(config, "relax_soft_filters_first", True)
            max_relaxation = getattr(config, "max_threshold_relaxation", 0.1)

    # First pass
    passing, failing = cascade.filter_candidates(candidates)

    stats = {
        "total_input": len(candidates),
        "passing_first_pass": len(passing),
        "failing_first_pass": len(failing),
        "relaxations_applied": [],
    }

    # Fallback if insufficient candidates
    if len(passing) < min_candidates:
        print(f"Warning: Only {len(passing)} candidates passed. Applying fallback...")

        if relax_soft_first:
            # Phase 1: Accept candidates that only have soft failures
            for candidate in failing:
                if len(passing) >= min_candidates:
                    break
                hard_fails = [k for k, v in candidate.filter_results.items() if v == FilterResult.FAIL]

                if not hard_fails:
                    candidate.risk_flags.append("ACCEPTED_VIA_FALLBACK_SOFT")
                    passing.append(candidate)
                    stats["relaxations_applied"].append(f"{candidate.candidate_id}: soft filter relaxation")

        # Phase 2: If still insufficient, relax thresholds (up to max_relaxation)
        if len(passing) < min_candidates and max_relaxation > 0:
            # Re-run with relaxed thresholds
            relaxed_thresholds = cascade._thresholds.copy()
            relaxed_thresholds["min_pdockq"] *= (1 - max_relaxation)
            relaxed_thresholds["min_interface_area"] *= (1 - max_relaxation)
            relaxed_thresholds["min_oasis_score"] *= (1 - max_relaxation)

            for candidate in failing:
                if len(passing) >= min_candidates:
                    break
                if candidate in passing:
                    continue

                # Check if candidate passes with relaxed thresholds
                passes_relaxed = True

                if candidate.pdockq is not None and candidate.pdockq < relaxed_thresholds["min_pdockq"]:
                    passes_relaxed = False
                if candidate.interface_area is not None and candidate.interface_area < relaxed_thresholds["min_interface_area"]:
                    passes_relaxed = False

                humanness = candidate.oasis_score_mean or candidate.oasis_score_vh
                if humanness is not None and humanness < relaxed_thresholds["min_oasis_score"]:
                    passes_relaxed = False

                if passes_relaxed:
                    candidate.risk_flags.append(f"ACCEPTED_VIA_FALLBACK_THRESHOLD_{int(max_relaxation*100)}PCT")
                    passing.append(candidate)
                    stats["relaxations_applied"].append(f"{candidate.candidate_id}: threshold relaxation ({int(max_relaxation*100)}%)")

        # Re-sort and re-rank
        passing.sort(key=lambda c: c.composite_score, reverse=True)
        for i, candidate in enumerate(passing):
            candidate.rank = i + 1

    stats["final_passing"] = len(passing)
    stats["used_fallback"] = len(stats["relaxations_applied"]) > 0

    return passing, stats
