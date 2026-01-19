"""Humanness scoring for antibody sequences using BioPhi.

Evaluates how "human-like" an antibody sequence is, which correlates
with reduced immunogenicity risk.
"""

from dataclasses import dataclass
from typing import Optional
import warnings


@dataclass
class HumannessReport:
    """Humanness analysis results."""

    sequence: str
    chain_type: str
    oasis_score: Optional[float]  # 0-1, higher is more human; None if scoring unavailable
    oasis_percentile: Optional[float] = None
    closest_human_germline: Optional[str] = None
    germline_identity: Optional[float] = None
    humanization_suggestions: Optional[list[dict]] = None

    @property
    def is_human_like(self) -> bool:
        """Check if sequence passes typical humanness threshold."""
        if self.oasis_score is None:
            return False  # Can't assess without score
        return self.oasis_score >= 0.8

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "oasis_score": self.oasis_score,
            "oasis_percentile": self.oasis_percentile,
            "closest_human_germline": self.closest_human_germline,
            "germline_identity": self.germline_identity,
            "is_human_like": self.is_human_like,
        }


@dataclass
class PairedHumannessReport:
    """Humanness analysis for VH/VL pair."""

    vh_report: HumannessReport
    vl_report: Optional[HumannessReport]
    mean_score: Optional[float]  # None if scoring unavailable

    @property
    def is_human_like(self) -> bool:
        """Check if both chains pass humanness threshold."""
        if self.mean_score is None:
            return False  # Can't assess without score
        if self.vl_report:
            return self.vh_report.is_human_like and self.vl_report.is_human_like
        return self.vh_report.is_human_like

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        result = {
            "vh": self.vh_report.to_dict(),
            "mean_score": self.mean_score,
            "is_human_like": self.is_human_like,
        }
        if self.vl_report:
            result["vl"] = self.vl_report.to_dict()
        return result


def score_humanness(
    sequence: str,
    chain_type: str = "H",
) -> HumannessReport:
    """Score humanness of an antibody sequence using BioPhi OASis.

    Args:
        sequence: Amino acid sequence (VH or VL).
        chain_type: 'H' for heavy chain, 'L' for light chain.

    Returns:
        HumannessReport with OASis score and related metrics.
    """
    try:
        from biophi.humanization.methods.sapiens import Sapiens
        from biophi.humanization.methods.oasis import OASis

        # Initialize OASis scorer
        oasis = OASis()
        chain_type_oasis = "H" if chain_type.upper() in ["H", "VH"] else "L"

        # Get OASis humanness score
        score = oasis.score_humanness(sequence, chain_type=chain_type_oasis)

        # Get closest germline
        try:
            germline_info = oasis.get_closest_germline(
                sequence, chain_type=chain_type_oasis
            )
            closest_germline = germline_info.get("germline")
            germline_identity = germline_info.get("identity")
        except Exception:
            closest_germline = None
            germline_identity = None

        return HumannessReport(
            sequence=sequence,
            chain_type=chain_type,
            oasis_score=score,
            closest_human_germline=closest_germline,
            germline_identity=germline_identity,
        )

    except ImportError:
        warnings.warn(
            "BioPhi not installed. Install with: pip install biophi. "
            "Humanness scoring will be skipped (soft-fail)."
        )
        # Return None for oasis_score to trigger soft-fail in filter cascade
        # This allows the pipeline to continue without BioPhi installed
        return HumannessReport(
            sequence=sequence,
            chain_type=chain_type,
            oasis_score=None,  # None triggers soft-fail, not hard-fail
        )
    except Exception as e:
        warnings.warn(f"Humanness scoring failed: {e}. Will soft-fail.")
        return HumannessReport(
            sequence=sequence,
            chain_type=chain_type,
            oasis_score=None,  # None triggers soft-fail
        )


def score_humanness_pair(
    vh_sequence: str,
    vl_sequence: Optional[str] = None,
) -> PairedHumannessReport:
    """Score humanness of a VH/VL pair.

    Args:
        vh_sequence: Heavy chain variable region sequence.
        vl_sequence: Light chain variable region sequence (optional for VHH).

    Returns:
        PairedHumannessReport with scores for both chains.
    """
    vh_report = score_humanness(vh_sequence, chain_type="H")

    if vl_sequence:
        vl_report = score_humanness(vl_sequence, chain_type="L")
        # Handle None scores (BioPhi unavailable)
        if vh_report.oasis_score is not None and vl_report.oasis_score is not None:
            mean_score = (vh_report.oasis_score + vl_report.oasis_score) / 2
        elif vh_report.oasis_score is not None:
            mean_score = vh_report.oasis_score
        elif vl_report.oasis_score is not None:
            mean_score = vl_report.oasis_score
        else:
            mean_score = None
    else:
        vl_report = None
        mean_score = vh_report.oasis_score  # May be None

    return PairedHumannessReport(
        vh_report=vh_report,
        vl_report=vl_report,
        mean_score=mean_score,
    )


def get_humanization_suggestions(
    sequence: str,
    chain_type: str = "H",
) -> list[dict]:
    """Get suggestions for humanizing a sequence using Sapiens.

    Args:
        sequence: Amino acid sequence to humanize.
        chain_type: 'H' for heavy, 'L' for light.

    Returns:
        List of mutation suggestions with positions and scores.
    """
    try:
        from biophi.humanization.methods.sapiens import Sapiens

        sapiens = Sapiens()
        chain_type_sapiens = "H" if chain_type.upper() in ["H", "VH"] else "L"

        # Get humanization suggestions
        result = sapiens.humanize(sequence, chain_type=chain_type_sapiens)

        suggestions = []
        if hasattr(result, "mutations"):
            for mut in result.mutations:
                suggestions.append({
                    "position": mut.position,
                    "original": mut.original,
                    "suggested": mut.suggested,
                    "score_change": getattr(mut, "score_change", None),
                })

        return suggestions

    except ImportError:
        warnings.warn("BioPhi not installed for humanization suggestions")
        return []
    except Exception as e:
        warnings.warn(f"Humanization suggestions failed: {e}")
        return []


def filter_by_humanness(
    sequences: list[str],
    chain_type: str = "H",
    min_score: float = 0.8,
) -> list[tuple[str, HumannessReport]]:
    """Filter sequences by humanness score.

    Args:
        sequences: List of sequences to filter.
        chain_type: Chain type for all sequences.
        min_score: Minimum OASis score to pass filter.

    Returns:
        List of (sequence, report) tuples that pass the filter.
    """
    results = []

    for seq in sequences:
        report = score_humanness(seq, chain_type)
        if report.oasis_score >= min_score:
            results.append((seq, report))

    return results
