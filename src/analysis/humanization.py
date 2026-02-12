"""Post-hoc humanization of antibody candidates using HuDiff.

Identifies near-miss candidates that failed humanness filtering and
generates humanized variants using HuDiff (diffusion-based framework
regeneration). CDR sequences are preserved.

Functions follow the module pattern (no classes) per project convention.
"""

from typing import Optional


def identify_near_misses(
    candidates: list,
    min_score: float = 0.70,
    max_score: float = 0.80,
) -> list:
    """Identify candidates that narrowly failed humanness filtering.

    Looks for candidates with humanness scores in the [min_score, max_score)
    range — close enough to threshold that framework humanization could push
    them past 0.80.

    Args:
        candidates: List of candidate dicts from step 04a output.
        min_score: Lower bound for humanness (below this, unlikely to recover).
        max_score: Upper bound (at or above this, already passes).

    Returns:
        List of near-miss candidate dicts.
    """
    near_misses = []
    for c in candidates:
        # Get humanness score from various possible locations
        humanness = _get_humanness_score(c)
        if humanness is None:
            continue

        if min_score <= humanness < max_score:
            near_misses.append(c)

    return near_misses


def _get_humanness_score(candidate: dict) -> Optional[float]:
    """Extract humanness score from a candidate dict.

    Checks multiple locations where humanness may be stored.
    """
    # Pre-computed in step 04a quick filter
    score = candidate.get("_oasis_score_mean")
    if score is not None:
        return float(score)

    # From validation_scores
    vs = candidate.get("validation_scores", {})
    score = vs.get("humanness_mean")
    if score is not None:
        return float(score)

    # Try computing from VH/VL scores
    vh_score = candidate.get("_oasis_score_vh")
    vl_score = candidate.get("_oasis_score_vl")
    if vh_score is not None:
        if vl_score is not None:
            return (float(vh_score) + float(vl_score)) / 2.0
        return float(vh_score)

    return None


def prepare_hudiff_input(candidate: dict) -> dict:
    """Prepare a candidate for HuDiff humanization.

    Extracts the right sequences and binder type information.

    Args:
        candidate: Candidate dict from pipeline.

    Returns:
        Dict with keys needed by HuDiff Modal functions.
    """
    design_id = candidate.get("design_id", candidate.get("name", "unknown"))
    binder_type = candidate.get("binder_type", "vhh")

    # Extract sequences
    vh_seq = candidate.get("sequence") or candidate.get("vh", "")
    vl_seq = candidate.get("sequence_vl") or candidate.get("vl")

    return {
        "design_id": design_id,
        "binder_type": binder_type,
        "sequence": vh_seq,
        "sequence_vl": vl_seq,
    }


def parse_hudiff_output(
    hudiff_result: dict,
    parent_candidate: dict,
) -> list:
    """Convert HuDiff output variants into candidate dicts.

    Each variant gets a new design_id and preserves parent metadata.

    Args:
        hudiff_result: Result dict from batch_humanize with 'variants' list.
        parent_candidate: Original candidate that was humanized.

    Returns:
        List of new candidate dicts (one per variant).
    """
    parent_id = parent_candidate.get("design_id", parent_candidate.get("name", "unknown"))
    binder_type = parent_candidate.get("binder_type", "vhh")
    variants = hudiff_result.get("variants", [])

    results = []
    for i, variant in enumerate(variants):
        new_candidate = create_humanized_candidate(
            parent=parent_candidate,
            new_vh=variant.get("vh") or variant.get("sequence", ""),
            new_vl=variant.get("vl"),
            variant_index=i,
        )
        results.append(new_candidate)

    return results


def create_humanized_candidate(
    parent: dict,
    new_vh: str,
    new_vl: Optional[str] = None,
    variant_index: int = 0,
) -> dict:
    """Create a new candidate entry from a humanized variant.

    Preserves parent metadata (target, source info) while replacing
    sequences with humanized versions. Clears stale structural predictions
    since framework changes require re-prediction.

    Args:
        parent: Original parent candidate dict.
        new_vh: New VH/VHH sequence from HuDiff.
        new_vl: New VL sequence from HuDiff (None for VHH).
        variant_index: Index of this variant among siblings.

    Returns:
        New candidate dict ready for structure prediction.
    """
    parent_id = parent.get("design_id", parent.get("name", "unknown"))
    binder_type = parent.get("binder_type", "vhh")

    new_id = f"{parent_id}_hudiff_{variant_index}"

    # Compute mutations list
    parent_vh = parent.get("sequence") or parent.get("vh", "")
    mutations = _compute_mutations(parent_vh, new_vh, "VH")
    if new_vl and parent.get("sequence_vl"):
        mutations.extend(_compute_mutations(
            parent.get("sequence_vl") or parent.get("vl", ""),
            new_vl, "VL"
        ))

    candidate = {
        "design_id": new_id,
        "sequence": new_vh,
        "vh": new_vh,
        "binder_type": binder_type,
        "source": "hudiff_humanized",
        "parent_design_id": parent_id,
        "humanization_mutations": mutations,
        "humanization_num_mutations": len(mutations),
        # Preserve target info from parent
        "target_name": parent.get("target_name"),
        "target_pdb": parent.get("target_pdb"),
        # Preserve BoltzGen metrics from parent for reference
        "parent_boltzgen_rank": parent.get("boltzgen_rank"),
        "parent_iptm": (parent.get("structure_prediction") or {}).get("iptm"),
    }

    if new_vl:
        candidate["sequence_vl"] = new_vl
        candidate["vl"] = new_vl

    return candidate


def _compute_mutations(
    original: str,
    mutated: str,
    chain_label: str,
) -> list:
    """Compute list of mutations between two sequences.

    Args:
        original: Original amino acid sequence.
        mutated: Mutated amino acid sequence.
        chain_label: Label for the chain (e.g., "VH", "VL").

    Returns:
        List of mutation dicts with position, original, mutated, chain.
    """
    mutations = []
    min_len = min(len(original), len(mutated))

    for i in range(min_len):
        if original[i] != mutated[i]:
            mutations.append({
                "position": i + 1,  # 1-indexed
                "original": original[i],
                "mutated": mutated[i],
                "chain": chain_label,
            })

    # Handle length differences
    if len(mutated) > len(original):
        for i in range(len(original), len(mutated)):
            mutations.append({
                "position": i + 1,
                "original": "-",
                "mutated": mutated[i],
                "chain": chain_label,
            })
    elif len(original) > len(mutated):
        for i in range(len(mutated), len(original)):
            mutations.append({
                "position": i + 1,
                "original": original[i],
                "mutated": "-",
                "chain": chain_label,
            })

    return mutations
