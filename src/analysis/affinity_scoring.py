"""Affinity proxy scoring using inverse folding models.

ProteinMPNN (MIT) and AntiFold (BSD-3) provide log-likelihood scores as
affinity proxies for de novo antibody designs. These are NOT direct affinity
predictors, but higher log-likelihoods correlate with better binding
(ProteinMPNN: Spearman r=0.27-0.41 on AbBiBench).

Results are informational — used in reports, not for filtering/ranking.

Usage:
    from src.analysis.affinity_scoring import batch_score_affinity
    results = batch_score_affinity("data/outputs/structures/cif/", design_ids, binder_types)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class AffinityResult:
    """Result from inverse folding affinity scoring."""

    design_id: str
    proteinmpnn_ll: Optional[float] = None  # Log-likelihood (higher = better fit)
    antifold_ll: Optional[float] = None  # Log-likelihood (higher = better fit)
    error: Optional[str] = None


def score_proteinmpnn(cif_path: str, design_id: str) -> AffinityResult:
    """Score a complex using ProteinMPNN log-likelihood.

    ProteinMPNN scores how well the binder sequence "fits" the predicted
    complex structure. Higher log-likelihood suggests better structural
    complementarity at the interface.

    Args:
        cif_path: Path to CIF structure file from Boltz-2.
        design_id: Identifier for the design.

    Returns:
        AffinityResult with proteinmpnn_ll or error.
    """
    if not Path(cif_path).exists():
        return AffinityResult(design_id=design_id, error=f"CIF file not found: {cif_path}")

    try:
        from proteinmpnn import score_complex  # noqa: F401
    except ImportError:
        return AffinityResult(
            design_id=design_id,
            error="ProteinMPNN not installed. Install with: pip install proteinmpnn",
        )

    try:
        # ProteinMPNN scoring mode: compute log-likelihood of the binder
        # sequence conditioned on the complex structure
        result = score_complex(cif_path)

        # Extract mean log-likelihood across binder residues
        ll = result.get("mean_log_likelihood", result.get("score"))
        if ll is None:
            return AffinityResult(design_id=design_id, error="No log-likelihood in ProteinMPNN output")

        return AffinityResult(design_id=design_id, proteinmpnn_ll=float(ll))

    except Exception as e:
        return AffinityResult(design_id=design_id, error=f"ProteinMPNN error: {e}")


def score_antifold(cif_path: str, design_id: str, binder_type: str = "vhh") -> AffinityResult:
    """Score a complex using AntiFold log-likelihood.

    AntiFold is an antibody-specific inverse folding model trained on
    paired antibody structures. It supports nanobodies via --nanobody_chain.

    Args:
        cif_path: Path to CIF structure file.
        design_id: Identifier for the design.
        binder_type: "vhh" for nanobody, "scfv"/"fab" for paired.

    Returns:
        AffinityResult with antifold_ll or error.
    """
    if not Path(cif_path).exists():
        return AffinityResult(design_id=design_id, error=f"CIF file not found: {cif_path}")

    try:
        import antifold  # noqa: F401
    except ImportError:
        return AffinityResult(
            design_id=design_id,
            error="AntiFold not installed. Install with: pip install antifold",
        )

    try:
        # AntiFold supports nanobody scoring with explicit flag
        if binder_type == "vhh":
            from antifold import score_nanobody
            result = score_nanobody(cif_path)
        else:
            from antifold import score_antibody
            result = score_antibody(cif_path)

        ll = result.get("mean_log_likelihood", result.get("score"))
        if ll is None:
            return AffinityResult(design_id=design_id, error="No log-likelihood in AntiFold output")

        return AffinityResult(design_id=design_id, antifold_ll=float(ll))

    except Exception as e:
        return AffinityResult(design_id=design_id, error=f"AntiFold error: {e}")


def batch_score_affinity(
    cif_dir: str,
    design_ids: list[str],
    binder_types: list[str],
    run_proteinmpnn: bool = True,
    run_antifold: bool = True,
) -> list[AffinityResult]:
    """Batch affinity scoring for multiple designs.

    Runs both ProteinMPNN and AntiFold on each design, merging results.
    Graceful error handling: returns error message if a tool is not installed.

    Args:
        cif_dir: Directory containing CIF files (named {design_id}.cif).
        design_ids: List of design identifiers.
        binder_types: Corresponding binder types ("vhh" or "scfv").
        run_proteinmpnn: Whether to run ProteinMPNN scoring.
        run_antifold: Whether to run AntiFold scoring.

    Returns:
        List of AffinityResult objects with merged scores.
    """
    cif_base = Path(cif_dir)
    results = []

    for design_id, binder_type in zip(design_ids, binder_types):
        cif_path = cif_base / f"{design_id}.cif"
        merged = AffinityResult(design_id=design_id)
        errors = []

        if run_proteinmpnn:
            mpnn_result = score_proteinmpnn(str(cif_path), design_id)
            if mpnn_result.proteinmpnn_ll is not None:
                merged.proteinmpnn_ll = mpnn_result.proteinmpnn_ll
            elif mpnn_result.error:
                errors.append(mpnn_result.error)

        if run_antifold:
            af_result = score_antifold(str(cif_path), design_id, binder_type)
            if af_result.antifold_ll is not None:
                merged.antifold_ll = af_result.antifold_ll
            elif af_result.error:
                errors.append(af_result.error)

        if errors and merged.proteinmpnn_ll is None and merged.antifold_ll is None:
            merged.error = "; ".join(errors)

        results.append(merged)

    return results
