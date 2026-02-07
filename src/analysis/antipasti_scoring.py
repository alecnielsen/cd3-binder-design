"""DEPRECATED: ANTIPASTI affinity prediction wrapper.

ANTIPASTI has been replaced by ProteinMPNN + AntiFold scoring in
`src/analysis/affinity_scoring.py`. ANTIPASTI was not integrated because:
- Requires R dependency (rpy2), adding complexity
- No VHH/nanobody support
- Performance degrades on predicted (vs experimental) structures

Use `src.analysis.affinity_scoring` instead:
    from src.analysis.affinity_scoring import batch_score_affinity

See also: scripts/05b_validate_candidates.py for the validation pipeline step.
"""

from dataclasses import dataclass
from typing import Optional

# Re-export for any code that imports from here
from src.analysis.affinity_scoring import AffinityResult, batch_score_affinity  # noqa: F401


@dataclass
class AntipastiResult:
    """DEPRECATED: Use AffinityResult from src.analysis.affinity_scoring instead."""

    design_id: str
    log10_kd: Optional[float] = None
    kd_nm: Optional[float] = None
    error: Optional[str] = None


def score_affinity(cif_path: str, design_id: str) -> AntipastiResult:
    """DEPRECATED: Use score_proteinmpnn() or score_antifold() from affinity_scoring."""
    return AntipastiResult(
        design_id=design_id,
        error="DEPRECATED: Use src.analysis.affinity_scoring instead. "
              "See src/analysis/affinity_scoring.py for ProteinMPNN + AntiFold scoring.",
    )
