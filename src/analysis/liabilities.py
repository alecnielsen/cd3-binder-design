"""Sequence liability detection for antibody developability assessment.

Identifies sequence motifs associated with:
- Deamidation (NG, NS, NT, ND, NH)
- Isomerization (DG, DS, DT, DD, DH, DN)
- N-glycosylation (N-X-S/T where X != P)
- Oxidation (exposed M, W)
- Unpaired cysteines
"""

import re
from dataclasses import dataclass
from typing import Optional

from src.utils.constants import (
    DEAMIDATION_MOTIFS,
    ISOMERIZATION_MOTIFS,
    OXIDATION_RESIDUES,
    GLYCOSYLATION_PATTERN,
)


@dataclass
class LiabilitySite:
    """A single liability site in a sequence."""

    motif: str
    position: int
    liability_type: str
    in_cdr: bool = False
    cdr_name: Optional[str] = None
    severity: str = "medium"  # low, medium, high

    def __str__(self) -> str:
        cdr_info = f" (in {self.cdr_name})" if self.in_cdr else ""
        return f"{self.liability_type}: {self.motif} at position {self.position}{cdr_info}"


@dataclass
class LiabilityReport:
    """Complete liability analysis for a sequence."""

    sequence: str
    deamidation_sites: list[LiabilitySite]
    isomerization_sites: list[LiabilitySite]
    glycosylation_sites: list[LiabilitySite]
    oxidation_sites: list[LiabilitySite]
    unpaired_cysteines: int
    total_liabilities: int
    cdr_liabilities: int  # Hard liabilities only (deamidation, isomerization, glycosylation)
    cdr_oxidation_count: int = 0  # Soft filter: oxidation in CDRs (flagged but not rejected)

    @property
    def has_cdr_liabilities(self) -> bool:
        """Check if any hard liabilities are in CDR regions.

        Note: Oxidation is a soft filter and not included here.
        """
        return self.cdr_liabilities > 0

    @property
    def is_clean(self) -> bool:
        """Check if sequence has no critical liabilities.

        Note: CDR oxidation sites are soft filters and don't affect this check.
        """
        return self.cdr_liabilities == 0 and self.unpaired_cysteines == 0

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "deamidation_sites": [
                {"motif": s.motif, "position": s.position, "in_cdr": s.in_cdr}
                for s in self.deamidation_sites
            ],
            "isomerization_sites": [
                {"motif": s.motif, "position": s.position, "in_cdr": s.in_cdr}
                for s in self.isomerization_sites
            ],
            "glycosylation_sites": [
                {"motif": s.motif, "position": s.position, "in_cdr": s.in_cdr}
                for s in self.glycosylation_sites
            ],
            "oxidation_sites": [
                {"motif": s.motif, "position": s.position, "in_cdr": s.in_cdr}
                for s in self.oxidation_sites
            ],
            "unpaired_cysteines": self.unpaired_cysteines,
            "total_liabilities": self.total_liabilities,
            "cdr_liabilities": self.cdr_liabilities,
            "cdr_oxidation_count": self.cdr_oxidation_count,
            "is_clean": self.is_clean,
        }


class LiabilityScanner:
    """Scan antibody sequences for developability liabilities."""

    def __init__(self, cdr_positions: Optional[dict[str, tuple[int, int]]] = None):
        """Initialize scanner.

        Args:
            cdr_positions: Dict mapping CDR names to (start, end) positions.
                          If None, CDR-specific analysis is disabled.
        """
        self.cdr_positions = cdr_positions or {}

    def _is_in_cdr(self, position: int) -> tuple[bool, Optional[str]]:
        """Check if a position falls within a CDR region."""
        for cdr_name, (start, end) in self.cdr_positions.items():
            if start <= position <= end:
                return True, cdr_name
        return False, None

    def find_deamidation_sites(self, sequence: str) -> list[LiabilitySite]:
        """Find potential deamidation sites (NG, NS, NT, ND, NH motifs)."""
        sites = []
        sequence = sequence.upper()

        for motif in DEAMIDATION_MOTIFS:
            for match in re.finditer(motif, sequence):
                pos = match.start()
                in_cdr, cdr_name = self._is_in_cdr(pos)
                severity = "high" if in_cdr else "medium"

                sites.append(
                    LiabilitySite(
                        motif=motif,
                        position=pos,
                        liability_type="deamidation",
                        in_cdr=in_cdr,
                        cdr_name=cdr_name,
                        severity=severity,
                    )
                )

        return sites

    def find_isomerization_sites(self, sequence: str) -> list[LiabilitySite]:
        """Find potential isomerization sites (DG, DS, DT, DD, DH, DN motifs)."""
        sites = []
        sequence = sequence.upper()

        for motif in ISOMERIZATION_MOTIFS:
            for match in re.finditer(motif, sequence):
                pos = match.start()
                in_cdr, cdr_name = self._is_in_cdr(pos)
                severity = "high" if in_cdr else "medium"

                sites.append(
                    LiabilitySite(
                        motif=motif,
                        position=pos,
                        liability_type="isomerization",
                        in_cdr=in_cdr,
                        cdr_name=cdr_name,
                        severity=severity,
                    )
                )

        return sites

    def find_glycosylation_sites(self, sequence: str) -> list[LiabilitySite]:
        """Find N-glycosylation sites (N-X-S/T where X != P)."""
        sites = []
        sequence = sequence.upper()

        for match in re.finditer(GLYCOSYLATION_PATTERN, sequence):
            pos = match.start()
            motif = match.group()
            in_cdr, cdr_name = self._is_in_cdr(pos)
            # Glycosylation in CDRs is almost always problematic
            severity = "high" if in_cdr else "low"

            sites.append(
                LiabilitySite(
                    motif=motif,
                    position=pos,
                    liability_type="glycosylation",
                    in_cdr=in_cdr,
                    cdr_name=cdr_name,
                    severity=severity,
                )
            )

        return sites

    def find_oxidation_sites(self, sequence: str) -> list[LiabilitySite]:
        """Find potential oxidation sites (M, W residues)."""
        sites = []
        sequence = sequence.upper()

        for residue in OXIDATION_RESIDUES:
            for i, aa in enumerate(sequence):
                if aa == residue:
                    in_cdr, cdr_name = self._is_in_cdr(i)
                    # Oxidation is mainly a concern in exposed CDR residues
                    severity = "medium" if in_cdr else "low"

                    sites.append(
                        LiabilitySite(
                            motif=residue,
                            position=i,
                            liability_type="oxidation",
                            in_cdr=in_cdr,
                            cdr_name=cdr_name,
                            severity=severity,
                        )
                    )

        return sites

    def count_unpaired_cysteines(self, sequence: str) -> int:
        """Count potentially unpaired cysteines (odd count suggests unpaired)."""
        cys_count = sequence.upper().count("C")
        return cys_count % 2

    def scan(self, sequence: str) -> LiabilityReport:
        """Perform complete liability scan on a sequence.

        Args:
            sequence: Amino acid sequence to scan.

        Returns:
            LiabilityReport with all identified liabilities.

        Note:
            cdr_liabilities counts only hard liabilities (deamidation, isomerization,
            glycosylation) in CDRs. Oxidation in CDRs is tracked separately in
            cdr_oxidation_count as a soft filter per the README specification.
        """
        deamidation = self.find_deamidation_sites(sequence)
        isomerization = self.find_isomerization_sites(sequence)
        glycosylation = self.find_glycosylation_sites(sequence)
        oxidation = self.find_oxidation_sites(sequence)
        unpaired_cys = self.count_unpaired_cysteines(sequence)

        all_sites = deamidation + isomerization + glycosylation + oxidation

        # Hard liabilities in CDRs (excludes oxidation which is a soft filter)
        hard_liability_sites = deamidation + isomerization + glycosylation
        cdr_hard_sites = [s for s in hard_liability_sites if s.in_cdr]

        # Soft filter: oxidation in CDRs (flagged but not rejected)
        cdr_oxidation = [s for s in oxidation if s.in_cdr]

        return LiabilityReport(
            sequence=sequence,
            deamidation_sites=deamidation,
            isomerization_sites=isomerization,
            glycosylation_sites=glycosylation,
            oxidation_sites=oxidation,
            unpaired_cysteines=unpaired_cys,
            total_liabilities=len(all_sites) + (1 if unpaired_cys else 0),
            cdr_liabilities=len(cdr_hard_sites),
            cdr_oxidation_count=len(cdr_oxidation),
        )

    def scan_with_cdr_detection(
        self,
        sequence: str,
        chain_type: str = "H",
        numbering_scheme: str = "imgt",
    ) -> LiabilityReport:
        """Scan sequence with automatic CDR detection via ANARCI.

        Args:
            sequence: Amino acid sequence to scan.
            chain_type: 'H' for heavy chain, 'L' for light chain.
            numbering_scheme: Numbering scheme ('imgt', 'chothia', 'kabat').

        Returns:
            LiabilityReport with CDR-aware analysis.
        """
        # Import here to avoid circular dependency
        from src.analysis.numbering import get_cdr_positions

        cdr_positions = get_cdr_positions(sequence, chain_type, numbering_scheme)
        self.cdr_positions = cdr_positions

        return self.scan(sequence)


def scan_sequence(
    sequence: str,
    cdr_positions: Optional[dict[str, tuple[int, int]]] = None,
) -> LiabilityReport:
    """Convenience function to scan a sequence for liabilities.

    Args:
        sequence: Amino acid sequence to scan.
        cdr_positions: Optional CDR positions for region-specific analysis.

    Returns:
        LiabilityReport with all identified liabilities.
    """
    scanner = LiabilityScanner(cdr_positions)
    return scanner.scan(sequence)


def filter_by_liabilities(
    sequences: list[str],
    max_cdr_liabilities: int = 0,
    allow_unpaired_cys: bool = False,
    cdr_positions_list: Optional[list[dict]] = None,
) -> list[tuple[str, LiabilityReport]]:
    """Filter sequences by liability criteria.

    Args:
        sequences: List of sequences to filter.
        max_cdr_liabilities: Maximum allowed CDR liabilities.
        allow_unpaired_cys: Whether to allow unpaired cysteines.
        cdr_positions_list: CDR positions for each sequence.

    Returns:
        List of (sequence, report) tuples that pass filters.
    """
    results = []
    cdr_positions_list = cdr_positions_list or [None] * len(sequences)

    for seq, cdr_pos in zip(sequences, cdr_positions_list):
        scanner = LiabilityScanner(cdr_pos)
        report = scanner.scan(seq)

        if report.cdr_liabilities > max_cdr_liabilities:
            continue
        if not allow_unpaired_cys and report.unpaired_cysteines > 0:
            continue

        results.append((seq, report))

    return results
