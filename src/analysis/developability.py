"""Combined developability assessment for antibody sequences.

Integrates multiple metrics:
- Sequence liabilities
- Humanness
- Physicochemical properties
- Aggregation propensity
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np

from src.analysis.liabilities import LiabilityScanner, LiabilityReport
from src.analysis.humanness import score_humanness_pair, PairedHumannessReport
from src.utils.constants import (
    KYTE_DOOLITTLE,
    AROMATIC,
    HYDROPHOBIC,
    CHARGED_POS,
    CHARGED_NEG,
    DEFAULT_FILTER_THRESHOLDS,
)


@dataclass
class PhysicochemicalProperties:
    """Physicochemical properties of a sequence."""

    length: int
    molecular_weight: float
    net_charge: float
    isoelectric_point: float
    gravy: float  # Grand average of hydropathy
    aromatic_fraction: float
    aliphatic_index: float

    def to_dict(self) -> dict:
        return {
            "length": self.length,
            "molecular_weight": self.molecular_weight,
            "net_charge": self.net_charge,
            "isoelectric_point": self.isoelectric_point,
            "gravy": self.gravy,
            "aromatic_fraction": self.aromatic_fraction,
            "aliphatic_index": self.aliphatic_index,
        }


@dataclass
class AggregationMetrics:
    """Aggregation propensity metrics."""

    hydrophobic_patches: int
    max_patch_length: int
    aromatic_clusters: int
    cdr_hydrophobicity: Optional[float] = None

    def to_dict(self) -> dict:
        return {
            "hydrophobic_patches": self.hydrophobic_patches,
            "max_patch_length": self.max_patch_length,
            "aromatic_clusters": self.aromatic_clusters,
            "cdr_hydrophobicity": self.cdr_hydrophobicity,
        }


@dataclass
class DevelopabilityReport:
    """Complete developability assessment."""

    sequence: str
    chain_type: str
    liability_report: LiabilityReport
    humanness_report: Optional[PairedHumannessReport]
    physicochemical: PhysicochemicalProperties
    aggregation: AggregationMetrics
    cdr_h3_length: Optional[int] = None
    composite_score: float = 0.0
    flags: list[str] = field(default_factory=list)
    passes_filters: bool = True

    def to_dict(self) -> dict:
        return {
            "chain_type": self.chain_type,
            "liabilities": self.liability_report.to_dict(),
            "humanness": self.humanness_report.to_dict() if self.humanness_report else None,
            "physicochemical": self.physicochemical.to_dict(),
            "aggregation": self.aggregation.to_dict(),
            "cdr_h3_length": self.cdr_h3_length,
            "composite_score": self.composite_score,
            "flags": self.flags,
            "passes_filters": self.passes_filters,
        }


class DevelopabilityAssessor:
    """Assess antibody developability using multiple metrics."""

    # Amino acid molecular weights (Da)
    MW = {
        "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
        "Q": 146.2, "E": 147.1, "G": 75.1, "H": 155.2, "I": 131.2,
        "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
        "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1,
    }

    # pKa values for isoelectric point calculation
    PKA = {
        "N_term": 9.69,
        "C_term": 2.34,
        "D": 3.86, "E": 4.25, "H": 6.00,
        "C": 8.33, "Y": 10.07, "K": 10.53, "R": 12.48,
    }

    def __init__(self, thresholds: Optional[dict] = None):
        """Initialize assessor with filtering thresholds.

        Args:
            thresholds: Dict of threshold values. Uses defaults if not provided.
        """
        self.thresholds = {**DEFAULT_FILTER_THRESHOLDS, **(thresholds or {})}

    def calculate_physicochemical(self, sequence: str) -> PhysicochemicalProperties:
        """Calculate physicochemical properties of a sequence."""
        sequence = sequence.upper()
        length = len(sequence)

        # Molecular weight
        mw = sum(self.MW.get(aa, 110) for aa in sequence) - (length - 1) * 18.015

        # Net charge at pH 7
        pos_charge = sum(1 for aa in sequence if aa in CHARGED_POS)
        neg_charge = sum(1 for aa in sequence if aa in CHARGED_NEG)
        net_charge = pos_charge - neg_charge

        # Isoelectric point (simplified Henderson-Hasselbalch)
        pi = self._calculate_pi(sequence)

        # GRAVY (Grand Average of Hydropathy)
        gravy = np.mean([KYTE_DOOLITTLE.get(aa, 0) for aa in sequence])

        # Aromatic fraction
        aromatic = sum(1 for aa in sequence if aa in AROMATIC) / length

        # Aliphatic index
        aliphatic = self._calculate_aliphatic_index(sequence)

        return PhysicochemicalProperties(
            length=length,
            molecular_weight=mw,
            net_charge=net_charge,
            isoelectric_point=pi,
            gravy=gravy,
            aromatic_fraction=aromatic,
            aliphatic_index=aliphatic,
        )

    def _calculate_pi(self, sequence: str) -> float:
        """Calculate isoelectric point using bisection method."""
        # Count titratable residues
        counts = {aa: sequence.count(aa) for aa in "DECHKRY"}

        def charge_at_ph(ph: float) -> float:
            # N-terminus
            charge = 1 / (1 + 10 ** (ph - self.PKA["N_term"]))
            # C-terminus
            charge -= 1 / (1 + 10 ** (self.PKA["C_term"] - ph))
            # Acidic residues (D, E)
            for aa in "DE":
                charge -= counts.get(aa, 0) / (1 + 10 ** (self.PKA[aa] - ph))
            # Basic residues (H, K, R)
            for aa in "HKR":
                charge += counts.get(aa, 0) / (1 + 10 ** (ph - self.PKA[aa]))
            # Cysteine and Tyrosine
            for aa in "CY":
                charge -= counts.get(aa, 0) / (1 + 10 ** (self.PKA[aa] - ph))
            return charge

        # Bisection to find pI
        low, high = 0.0, 14.0
        while high - low > 0.01:
            mid = (low + high) / 2
            if charge_at_ph(mid) > 0:
                low = mid
            else:
                high = mid

        return (low + high) / 2

    def _calculate_aliphatic_index(self, sequence: str) -> float:
        """Calculate aliphatic index."""
        length = len(sequence)
        if length == 0:
            return 0

        a = sequence.count("A") / length * 100
        v = sequence.count("V") / length * 100
        i = sequence.count("I") / length * 100
        l = sequence.count("L") / length * 100

        return a + 2.9 * v + 3.9 * (i + l)

    def calculate_aggregation_metrics(
        self,
        sequence: str,
        cdr_sequences: Optional[list[str]] = None,
    ) -> AggregationMetrics:
        """Calculate aggregation propensity metrics."""
        sequence = sequence.upper()

        # Count hydrophobic patches (stretches of 5+ hydrophobic residues)
        patches = 0
        max_patch = 0
        current_patch = 0

        for aa in sequence:
            if aa in HYDROPHOBIC:
                current_patch += 1
                if current_patch >= 5:
                    if current_patch == 5:
                        patches += 1
                    max_patch = max(max_patch, current_patch)
            else:
                current_patch = 0

        # Count aromatic clusters (consecutive aromatic residues)
        aromatic_clusters = 0
        prev_aromatic = False
        for aa in sequence:
            is_aromatic = aa in AROMATIC
            if is_aromatic and prev_aromatic:
                aromatic_clusters += 1
            prev_aromatic = is_aromatic

        # CDR hydrophobicity
        cdr_hydro = None
        if cdr_sequences:
            all_cdr = "".join(cdr_sequences)
            if all_cdr:
                cdr_hydro = np.mean([KYTE_DOOLITTLE.get(aa, 0) for aa in all_cdr])

        return AggregationMetrics(
            hydrophobic_patches=patches,
            max_patch_length=max_patch,
            aromatic_clusters=aromatic_clusters,
            cdr_hydrophobicity=cdr_hydro,
        )

    def assess(
        self,
        vh_sequence: str,
        vl_sequence: Optional[str] = None,
        include_humanness: bool = True,
        cdr_positions: Optional[dict] = None,
    ) -> DevelopabilityReport:
        """Perform complete developability assessment.

        Args:
            vh_sequence: Heavy chain variable region sequence.
            vl_sequence: Light chain variable region (optional for VHH).
            include_humanness: Whether to run humanness scoring.
            cdr_positions: Pre-computed CDR positions.

        Returns:
            DevelopabilityReport with all metrics.
        """
        # Determine sequence type
        chain_type = "VHH" if vl_sequence is None else "Fv"
        full_sequence = vh_sequence + (vl_sequence or "")

        # Liability scanning
        scanner = LiabilityScanner(cdr_positions)
        liability_report = scanner.scan(full_sequence)

        # Humanness scoring
        humanness_report = None
        if include_humanness:
            humanness_report = score_humanness_pair(vh_sequence, vl_sequence)

        # Physicochemical properties
        physicochemical = self.calculate_physicochemical(full_sequence)

        # CDR sequences for aggregation assessment
        cdr_sequences = None
        if cdr_positions:
            cdr_sequences = [
                full_sequence[start:end+1]
                for start, end in cdr_positions.values()
                if start < len(full_sequence)
            ]

        aggregation = self.calculate_aggregation_metrics(full_sequence, cdr_sequences)

        # Extract CDR-H3 length if available
        cdr_h3_length = None
        if cdr_positions and "H3" in cdr_positions:
            h3_start, h3_end = cdr_positions["H3"]
            cdr_h3_length = h3_end - h3_start + 1

        # Compute flags and filter status
        flags = []
        passes = True

        # Check liabilities
        if liability_report.cdr_liabilities > 0:
            flags.append(f"CDR liabilities: {liability_report.cdr_liabilities}")
            passes = False

        if liability_report.unpaired_cysteines > 0:
            flags.append("Unpaired cysteine")
            passes = False

        # Check humanness
        if humanness_report and humanness_report.mean_score < self.thresholds["min_oasis_score"]:
            flags.append(f"Low humanness: {humanness_report.mean_score:.2f}")

        # Check physicochemical
        if not (self.thresholds["net_charge_min"] <= physicochemical.net_charge <= self.thresholds["net_charge_max"]):
            flags.append(f"Extreme charge: {physicochemical.net_charge}")

        if not (self.thresholds["pi_min"] <= physicochemical.isoelectric_point <= self.thresholds["pi_max"]):
            flags.append(f"Extreme pI: {physicochemical.isoelectric_point:.1f}")

        # Check aggregation
        if aggregation.hydrophobic_patches > self.thresholds["max_hydrophobic_patches"]:
            flags.append(f"Hydrophobic patches: {aggregation.hydrophobic_patches}")

        # Check CDR-H3 length
        if cdr_h3_length:
            if not (self.thresholds["cdr_h3_length_min"] <= cdr_h3_length <= self.thresholds["cdr_h3_length_max"]):
                flags.append(f"Unusual CDR-H3 length: {cdr_h3_length}")

        # Compute composite score
        composite_score = self._compute_composite_score(
            liability_report, humanness_report, physicochemical, aggregation
        )

        return DevelopabilityReport(
            sequence=full_sequence,
            chain_type=chain_type,
            liability_report=liability_report,
            humanness_report=humanness_report,
            physicochemical=physicochemical,
            aggregation=aggregation,
            cdr_h3_length=cdr_h3_length,
            composite_score=composite_score,
            flags=flags,
            passes_filters=passes,
        )

    def _compute_composite_score(
        self,
        liability: LiabilityReport,
        humanness: Optional[PairedHumannessReport],
        physicochemical: PhysicochemicalProperties,
        aggregation: AggregationMetrics,
    ) -> float:
        """Compute composite developability score (0-1, higher is better)."""
        scores = []
        weights = []

        # Liability score (25% weight)
        liability_score = max(0, 1 - liability.cdr_liabilities * 0.2 - liability.total_liabilities * 0.05)
        if liability.unpaired_cysteines > 0:
            liability_score *= 0.5
        scores.append(liability_score)
        weights.append(0.25)

        # Humanness score (25% weight)
        if humanness:
            scores.append(humanness.mean_score)
            weights.append(0.25)

        # Physicochemical score (25% weight)
        physico_score = 1.0
        # Penalize extreme charge
        if abs(physicochemical.net_charge) > 5:
            physico_score -= 0.2
        # Penalize extreme pI
        if physicochemical.isoelectric_point < 5 or physicochemical.isoelectric_point > 10:
            physico_score -= 0.2
        # Penalize high hydrophobicity
        if physicochemical.gravy > 0:
            physico_score -= min(0.3, physicochemical.gravy * 0.1)
        scores.append(max(0, physico_score))
        weights.append(0.25)

        # Aggregation score (25% weight)
        agg_score = 1.0
        agg_score -= aggregation.hydrophobic_patches * 0.15
        agg_score -= aggregation.aromatic_clusters * 0.1
        scores.append(max(0, agg_score))
        weights.append(0.25)

        # Weighted average
        total_weight = sum(weights)
        return sum(s * w for s, w in zip(scores, weights)) / total_weight


def assess_developability(
    vh_sequence: str,
    vl_sequence: Optional[str] = None,
    thresholds: Optional[dict] = None,
) -> DevelopabilityReport:
    """Convenience function for developability assessment.

    Args:
        vh_sequence: Heavy chain variable region.
        vl_sequence: Light chain variable region (optional).
        thresholds: Custom filtering thresholds.

    Returns:
        DevelopabilityReport with all metrics.
    """
    assessor = DevelopabilityAssessor(thresholds)
    return assessor.assess(vh_sequence, vl_sequence)
