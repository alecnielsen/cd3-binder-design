"""Binding interface analysis for antibody-antigen complexes.

Provides metrics for assessing binding quality and
comparing epitopes between different binders.
"""

from dataclasses import dataclass, field
from typing import Optional
from pathlib import Path


@dataclass
class InterfaceMetrics:
    """Metrics characterizing a binding interface."""

    # Basic counts
    num_contacts: int = 0
    num_interface_residues_binder: int = 0
    num_interface_residues_target: int = 0

    # Areas
    interface_area: float = 0.0  # Buried surface area (Å²)

    # Contact types
    num_hydrogen_bonds: int = 0
    num_salt_bridges: int = 0
    num_hydrophobic_contacts: int = 0

    # Shape complementarity
    shape_complementarity: float = 0.0

    # Residue lists
    interface_residues_binder: list[int] = field(default_factory=list)
    interface_residues_target: list[int] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "num_contacts": self.num_contacts,
            "num_interface_residues_binder": self.num_interface_residues_binder,
            "num_interface_residues_target": self.num_interface_residues_target,
            "interface_area": self.interface_area,
            "num_hydrogen_bonds": self.num_hydrogen_bonds,
            "num_salt_bridges": self.num_salt_bridges,
            "num_hydrophobic_contacts": self.num_hydrophobic_contacts,
            "shape_complementarity": self.shape_complementarity,
            "interface_residues_binder": self.interface_residues_binder,
            "interface_residues_target": self.interface_residues_target,
        }


@dataclass
class EpitopeComparison:
    """Comparison of epitopes between two binders."""

    binder_1_name: str
    binder_2_name: str

    # Overlap metrics
    target_residue_overlap: list[int] = field(default_factory=list)
    overlap_count: int = 0
    overlap_fraction: float = 0.0  # Jaccard similarity

    # Individual epitopes
    epitope_1: list[int] = field(default_factory=list)
    epitope_2: list[int] = field(default_factory=list)

    # Classification
    epitope_class: str = "unknown"  # "same", "overlapping", "distinct"

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "binder_1_name": self.binder_1_name,
            "binder_2_name": self.binder_2_name,
            "target_residue_overlap": self.target_residue_overlap,
            "overlap_count": self.overlap_count,
            "overlap_fraction": self.overlap_fraction,
            "epitope_1": self.epitope_1,
            "epitope_2": self.epitope_2,
            "epitope_class": self.epitope_class,
        }


class InterfaceAnalyzer:
    """Analyzer for protein-protein binding interfaces."""

    # Cached OKT3 epitope residues (extracted from 1SY6)
    _cached_okt3_epitope: list[int] = None

    @classmethod
    def get_okt3_epitope(cls, force_refresh: bool = False) -> list[int]:
        """Get OKT3 epitope residues, extracting from 1SY6 if needed.

        This method dynamically extracts the OKT3 epitope from the 1SY6
        crystal structure rather than using hardcoded values. The result
        is cached for efficiency.

        Args:
            force_refresh: If True, re-extract from 1SY6 even if cached.

        Returns:
            Sorted list of CD3ε residue numbers that form the OKT3 epitope.
        """
        if cls._cached_okt3_epitope is None or force_refresh:
            try:
                from src.structure.pdb_utils import get_okt3_epitope_from_1sy6
                cls._cached_okt3_epitope = get_okt3_epitope_from_1sy6()
            except Exception as e:
                import warnings
                warnings.warn(
                    f"Failed to extract OKT3 epitope from 1SY6: {e}. "
                    "Using fallback hardcoded values."
                )
                # Fallback to hardcoded values only if extraction fails
                cls._cached_okt3_epitope = [
                    23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 38, 39, 40, 41, 42, 45, 47
                ]
        return cls._cached_okt3_epitope

    def __init__(
        self,
        contact_distance: float = 5.0,
        okt3_epitope_residues: list[int] = None,
    ):
        """Initialize analyzer.

        Args:
            contact_distance: Distance cutoff for contacts (Å).
            okt3_epitope_residues: Custom OKT3 epitope residues. If None,
                dynamically extracts from 1SY6 structure.
        """
        self.contact_distance = contact_distance
        # Use provided residues or dynamically extract from 1SY6
        if okt3_epitope_residues is not None:
            self.okt3_epitope_residues = okt3_epitope_residues
        else:
            self.okt3_epitope_residues = self.get_okt3_epitope()

    def analyze_interface(
        self,
        pdb_string: str,
        binder_chain: str = "B",
        target_chain: str = "A",
    ) -> InterfaceMetrics:
        """Analyze binding interface between two chains.

        Args:
            pdb_string: PDB content as string.
            binder_chain: Chain ID for binder.
            target_chain: Chain ID for target.

        Returns:
            InterfaceMetrics object.
        """
        from src.structure.pdb_utils import (
            calculate_interface_residues,
            count_contacts,
            estimate_interface_area,
        )

        # Get interface residues
        binder_residues, target_residues = calculate_interface_residues(
            pdb_string, binder_chain, target_chain, self.contact_distance
        )

        # Count contacts
        num_contacts = count_contacts(
            pdb_string, binder_chain, target_chain, self.contact_distance
        )

        # Estimate interface area
        interface_area = estimate_interface_area(
            pdb_string, binder_chain, target_chain
        )

        return InterfaceMetrics(
            num_contacts=num_contacts,
            num_interface_residues_binder=len(binder_residues),
            num_interface_residues_target=len(target_residues),
            interface_area=interface_area,
            interface_residues_binder=binder_residues,
            interface_residues_target=target_residues,
        )

    def compare_epitopes(
        self,
        epitope_1: list[int],
        epitope_2: list[int],
        name_1: str = "binder_1",
        name_2: str = "binder_2",
    ) -> EpitopeComparison:
        """Compare two epitopes.

        Args:
            epitope_1: Residue positions for first epitope.
            epitope_2: Residue positions for second epitope.
            name_1: Name for first binder.
            name_2: Name for second binder.

        Returns:
            EpitopeComparison object.
        """
        set_1 = set(epitope_1)
        set_2 = set(epitope_2)

        overlap = set_1 & set_2
        union = set_1 | set_2

        jaccard = len(overlap) / len(union) if union else 0.0

        # Classify
        if jaccard > 0.8:
            epitope_class = "same"
        elif jaccard > 0.3:
            epitope_class = "overlapping"
        else:
            epitope_class = "distinct"

        return EpitopeComparison(
            binder_1_name=name_1,
            binder_2_name=name_2,
            target_residue_overlap=sorted(overlap),
            overlap_count=len(overlap),
            overlap_fraction=jaccard,
            epitope_1=sorted(epitope_1),
            epitope_2=sorted(epitope_2),
            epitope_class=epitope_class,
        )

    def compare_to_okt3(
        self,
        epitope: list[int],
        binder_name: str = "design",
    ) -> EpitopeComparison:
        """Compare an epitope to the OKT3 epitope.

        Args:
            epitope: Residue positions for the epitope.
            binder_name: Name for the binder.

        Returns:
            EpitopeComparison with OKT3.
        """
        return self.compare_epitopes(
            epitope_1=epitope,
            epitope_2=self.okt3_epitope_residues,
            name_1=binder_name,
            name_2="OKT3",
        )

    def annotate_epitope_class(
        self,
        epitope: list[int],
        overlap_threshold: float = 0.5,
    ) -> tuple[str, float]:
        """Annotate whether an epitope is OKT3-like or novel.

        Args:
            epitope: Residue positions for the epitope.
            overlap_threshold: Threshold for OKT3-like classification.

        Returns:
            Tuple of (class, overlap_fraction).
        """
        comparison = self.compare_to_okt3(epitope)

        if comparison.overlap_fraction >= overlap_threshold:
            return "OKT3-like", comparison.overlap_fraction
        else:
            return "novel_epitope", comparison.overlap_fraction


def analyze_complex_interface(
    pdb_path: str,
    binder_chain: str = "B",
    target_chain: str = "A",
    compare_to_okt3: bool = True,
) -> dict:
    """Convenience function for interface analysis.

    Args:
        pdb_path: Path to complex PDB file.
        binder_chain: Chain ID for binder.
        target_chain: Chain ID for target.
        compare_to_okt3: If True, compare epitope to OKT3.

    Returns:
        Dictionary with analysis results.
    """
    with open(pdb_path, "r") as f:
        pdb_string = f.read()

    analyzer = InterfaceAnalyzer()
    metrics = analyzer.analyze_interface(pdb_string, binder_chain, target_chain)

    result = {
        "interface_metrics": metrics.to_dict(),
    }

    if compare_to_okt3:
        epitope_class, overlap = analyzer.annotate_epitope_class(
            metrics.interface_residues_target
        )
        result["epitope_annotation"] = {
            "epitope_class": epitope_class,
            "okt3_overlap_fraction": overlap,
            "target_contact_residues": metrics.interface_residues_target,
            "okt3_epitope_residues": analyzer.okt3_epitope_residues,
        }

    return result


def batch_epitope_annotation(
    complex_pdbs: list[str],
    binder_names: list[str],
    binder_chain: str = "B",
    target_chain: str = "A",
    overlap_threshold: float = 0.5,
) -> list[dict]:
    """Annotate epitopes for multiple complexes.

    Args:
        complex_pdbs: List of paths to complex PDB files.
        binder_names: Names for each binder.
        binder_chain: Chain ID for binder.
        target_chain: Chain ID for target.
        overlap_threshold: Threshold for OKT3-like classification.

    Returns:
        List of annotation dictionaries.
    """
    analyzer = InterfaceAnalyzer()
    results = []

    for pdb_path, name in zip(complex_pdbs, binder_names):
        try:
            with open(pdb_path, "r") as f:
                pdb_string = f.read()

            metrics = analyzer.analyze_interface(pdb_string, binder_chain, target_chain)
            epitope_class, overlap = analyzer.annotate_epitope_class(
                metrics.interface_residues_target, overlap_threshold
            )

            results.append({
                "binder_name": name,
                "epitope_class": epitope_class,
                "okt3_overlap_fraction": overlap,
                "target_contact_residues": metrics.interface_residues_target,
                "num_contacts": metrics.num_contacts,
                "interface_area": metrics.interface_area,
            })

        except Exception as e:
            print(f"Warning: Failed to analyze {name}: {e}")
            results.append({
                "binder_name": name,
                "error": str(e),
            })

    return results
