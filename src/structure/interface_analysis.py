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
    """Analyzer for protein-protein binding interfaces.

    IMPORTANT: This analyzer uses canonical CD3ε numbering (matching 1XIW
    chain A) for OKT3 epitope comparisons. When comparing epitopes from
    predicted complexes, ensure the target chain uses the same numbering
    scheme, or provide the target sequence for alignment-based comparison.
    """

    # Cached OKT3 epitope residues (in canonical CD3ε numbering)
    _cached_okt3_epitope: list[int] = None
    # Cached canonical CD3ε sequence (from 1XIW chain A)
    _cached_cd3e_sequence: str = None

    @classmethod
    def get_okt3_epitope(cls, force_refresh: bool = False) -> list[int]:
        """Get OKT3 epitope residues in canonical CD3ε numbering.

        This method extracts the OKT3 epitope from the 1SY6 crystal structure
        and maps it to canonical CD3ε numbering (1XIW chain A) via sequence
        alignment. The result is cached for efficiency.

        Args:
            force_refresh: If True, re-extract from 1SY6 even if cached.

        Returns:
            Sorted list of CD3ε residue numbers (canonical 1XIW numbering)
            that form the OKT3 epitope.
        """
        if cls._cached_okt3_epitope is None or force_refresh:
            try:
                from src.structure.pdb_utils import get_okt3_epitope_from_1sy6
                # map_to_canonical=True ensures we get 1XIW-compatible numbering
                cls._cached_okt3_epitope = get_okt3_epitope_from_1sy6(
                    map_to_canonical=True
                )
            except Exception as e:
                import warnings
                warnings.warn(
                    f"Failed to extract OKT3 epitope from 1SY6: {e}. "
                    "Using fallback hardcoded values (canonical CD3ε numbering)."
                )
                # Fallback: canonical CD3ε numbering from literature
                cls._cached_okt3_epitope = [
                    23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 38, 39, 40, 41, 42, 45, 47
                ]
        return cls._cached_okt3_epitope

    @classmethod
    def get_canonical_cd3e_sequence(cls, force_refresh: bool = False) -> str:
        """Get canonical CD3ε sequence from 1XIW chain A.

        Args:
            force_refresh: If True, re-extract even if cached.

        Returns:
            CD3ε amino acid sequence.
        """
        if cls._cached_cd3e_sequence is None or force_refresh:
            try:
                from src.structure.pdb_utils import (
                    download_pdb,
                    extract_sequence_with_numbering,
                )
                import os

                cache_dir = "data/targets"
                cached_path = os.path.join(cache_dir, "1XIW.pdb")

                if os.path.exists(cached_path):
                    with open(cached_path, "r") as f:
                        pdb_content = f.read()
                else:
                    os.makedirs(cache_dir, exist_ok=True)
                    pdb_content = download_pdb("1XIW", cached_path)

                # Chain A (or E) is CD3ε in 1XIW; chain D is UCHT1 VH
                seq, _ = extract_sequence_with_numbering(pdb_content, "A")
                cls._cached_cd3e_sequence = seq
            except Exception as e:
                import warnings
                warnings.warn(f"Failed to get CD3ε sequence: {e}")
                cls._cached_cd3e_sequence = ""

        return cls._cached_cd3e_sequence

    def __init__(
        self,
        contact_distance: float = 5.0,
        okt3_epitope_residues: list[int] = None,
        canonical_cd3e_sequence: str = None,
    ):
        """Initialize analyzer.

        Args:
            contact_distance: Distance cutoff for contacts (Å).
            okt3_epitope_residues: Custom OKT3 epitope residues. If None,
                dynamically extracts from 1SY6 (in canonical numbering).
            canonical_cd3e_sequence: Canonical CD3ε sequence for alignment.
                If None, uses 1XIW chain A sequence.
        """
        self.contact_distance = contact_distance
        # Use provided residues or dynamically extract from 1SY6
        if okt3_epitope_residues is not None:
            self.okt3_epitope_residues = okt3_epitope_residues
        else:
            self.okt3_epitope_residues = self.get_okt3_epitope()

        # Store canonical sequence for alignment-based comparisons
        if canonical_cd3e_sequence is not None:
            self.canonical_cd3e_sequence = canonical_cd3e_sequence
        else:
            self.canonical_cd3e_sequence = self.get_canonical_cd3e_sequence()

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

        IMPORTANT: This comparison assumes both epitopes use the SAME residue
        numbering scheme. For comparing epitopes from different structures
        (e.g., 1XIW vs 1SY6), use compare_epitopes_aligned() instead.

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

    def compare_epitopes_aligned(
        self,
        epitope_1: list[int],
        epitope_2: list[int],
        sequence_1: str,
        sequence_2: str,
        name_1: str = "binder_1",
        name_2: str = "binder_2",
        seq1_start: int = 1,
        seq2_start: int = 1,
    ) -> EpitopeComparison:
        """Compare epitopes from different structures using sequence alignment.

        This method should be used when comparing epitopes from structures
        with different residue numbering (e.g., 1XIW vs 1SY6). It aligns
        the target sequences and maps residue numbers to a common reference.

        NOTE: This method assumes SEQUENTIAL residue numbering from seq1_start
        and seq2_start. It does not handle PDB files with numbering gaps or
        insertion codes. For predicted structures (which typically use 1-indexed
        sequential numbering), this is usually appropriate.

        Args:
            epitope_1: Residue positions for first epitope (in seq1 numbering).
                Must use sequential numbering from seq1_start.
            epitope_2: Residue positions for second epitope (in seq2 numbering).
                Must use sequential numbering from seq2_start.
            sequence_1: Target sequence for first epitope.
            sequence_2: Target sequence for second epitope.
            name_1: Name for first binder.
            name_2: Name for second binder.
            seq1_start: Starting residue number for sequence 1 (default 1).
            seq2_start: Starting residue number for sequence 2 (default 1).

        Returns:
            EpitopeComparison with epitopes mapped to common numbering.
        """
        from src.structure.pdb_utils import align_residue_numbers

        # Map epitope_2 residues to sequence_1 numbering
        mapped_epitope_2 = align_residue_numbers(
            query_sequence=sequence_1,
            reference_sequence=sequence_2,
            reference_residue_numbers=epitope_2,
            query_start=seq1_start,
            reference_start=seq2_start,
        )

        if not mapped_epitope_2:
            import warnings
            warnings.warn(
                f"Failed to align epitopes between {name_1} and {name_2}. "
                "Falling back to direct comparison (may be unreliable)."
            )
            mapped_epitope_2 = epitope_2

        return self.compare_epitopes(
            epitope_1=epitope_1,
            epitope_2=mapped_epitope_2,
            name_1=name_1,
            name_2=name_2,
        )

    def compare_to_okt3(
        self,
        epitope: list[int],
        binder_name: str = "design",
        target_sequence: str = None,
    ) -> EpitopeComparison:
        """Compare an epitope to the OKT3 epitope.

        IMPORTANT: The OKT3 epitope uses canonical CD3ε numbering (1XIW chain A).
        If the input epitope uses different numbering (e.g., from a predicted
        structure with 1-indexed residue numbers), provide target_sequence
        to enable alignment-based comparison.

        Args:
            epitope: Residue positions for the epitope.
            binder_name: Name for the binder.
            target_sequence: Target sequence for alignment. If provided,
                maps the input epitope to canonical numbering before comparison.

        Returns:
            EpitopeComparison with OKT3.
        """
        if target_sequence and self.canonical_cd3e_sequence:
            # Use alignment-based comparison
            return self.compare_epitopes_aligned(
                epitope_1=self.okt3_epitope_residues,
                epitope_2=epitope,
                sequence_1=self.canonical_cd3e_sequence,
                sequence_2=target_sequence,
                name_1="OKT3",
                name_2=binder_name,
                seq1_start=1,  # Canonical numbering typically starts at 1
                seq2_start=1,  # Predicted structures typically start at 1
            )
        else:
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
        target_sequence: str = None,
    ) -> tuple[str, float]:
        """Annotate whether an epitope is OKT3-like or novel.

        IMPORTANT: The OKT3 epitope uses canonical CD3ε numbering (1XIW chain A).
        If the input epitope uses different numbering (e.g., from a predicted
        structure), provide target_sequence to enable alignment-based comparison.

        Args:
            epitope: Residue positions for the epitope.
            overlap_threshold: Threshold for OKT3-like classification.
            target_sequence: Target sequence for alignment-based comparison.

        Returns:
            Tuple of (class, overlap_fraction).
        """
        comparison = self.compare_to_okt3(
            epitope,
            binder_name="design",
            target_sequence=target_sequence,
        )

        if comparison.overlap_fraction >= overlap_threshold:
            return "OKT3-like", comparison.overlap_fraction
        else:
            return "novel_epitope", comparison.overlap_fraction


def analyze_complex_interface(
    pdb_path: str,
    binder_chain: str = "B",
    target_chain: str = "A",
    compare_to_okt3: bool = True,
    target_sequence: str = None,
) -> dict:
    """Convenience function for interface analysis.

    Args:
        pdb_path: Path to complex PDB file.
        binder_chain: Chain ID for binder.
        target_chain: Chain ID for target.
        compare_to_okt3: If True, compare epitope to OKT3.
        target_sequence: Target sequence for alignment-based OKT3 comparison.
            If None, extracts from the PDB file. Provide this when the PDB
            uses non-canonical residue numbering (e.g., predicted structures).

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
        # Extract target sequence if not provided (for alignment-based comparison)
        if target_sequence is None:
            from src.structure.pdb_utils import extract_sequence_with_numbering
            target_sequence, _ = extract_sequence_with_numbering(
                pdb_string, target_chain
            )

        epitope_class, overlap = analyzer.annotate_epitope_class(
            metrics.interface_residues_target,
            target_sequence=target_sequence,
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

    Uses sequence alignment for OKT3 epitope comparison to handle
    predicted structures with different residue numbering.

    Args:
        complex_pdbs: List of paths to complex PDB files.
        binder_names: Names for each binder.
        binder_chain: Chain ID for binder.
        target_chain: Chain ID for target.
        overlap_threshold: Threshold for OKT3-like classification.

    Returns:
        List of annotation dictionaries.
    """
    from src.structure.pdb_utils import extract_sequence_with_numbering

    analyzer = InterfaceAnalyzer()
    results = []

    for pdb_path, name in zip(complex_pdbs, binder_names):
        try:
            with open(pdb_path, "r") as f:
                pdb_string = f.read()

            metrics = analyzer.analyze_interface(pdb_string, binder_chain, target_chain)

            # Extract target sequence for alignment-based comparison
            target_sequence, _ = extract_sequence_with_numbering(
                pdb_string, target_chain
            )

            epitope_class, overlap = analyzer.annotate_epitope_class(
                metrics.interface_residues_target,
                overlap_threshold,
                target_sequence=target_sequence,
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
