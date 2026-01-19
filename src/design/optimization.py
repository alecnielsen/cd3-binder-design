"""Optimization of existing CD3 binders.

This module provides functions for generating optimized variants
of known CD3 binders (teplizumab, SP34, UCHT1) through:
- CDR grafting onto human frameworks
- Humanization mutations
- Back-mutation variants for humanization safety
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import yaml
import copy


@dataclass
class AntibodySequences:
    """VH and VL sequences for an antibody."""

    name: str
    vh: str
    vl: Optional[str] = None  # None for VHH
    source: str = "unknown"
    notes: str = ""

    # CDR sequences (extracted by ANARCI)
    cdr_h1: Optional[str] = None
    cdr_h2: Optional[str] = None
    cdr_h3: Optional[str] = None
    cdr_l1: Optional[str] = None
    cdr_l2: Optional[str] = None
    cdr_l3: Optional[str] = None

    def is_vhh(self) -> bool:
        """Check if this is a VHH (no VL)."""
        return self.vl is None

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "name": self.name,
            "vh": self.vh,
            "vl": self.vl,
            "source": self.source,
            "notes": self.notes,
            "cdr_h1": self.cdr_h1,
            "cdr_h2": self.cdr_h2,
            "cdr_h3": self.cdr_h3,
            "cdr_l1": self.cdr_l1,
            "cdr_l2": self.cdr_l2,
            "cdr_l3": self.cdr_l3,
        }


@dataclass
class OptimizedVariant:
    """An optimized variant of a starting antibody."""

    name: str
    parent_name: str
    vh: str
    vl: Optional[str] = None
    variant_type: str = "unknown"  # "humanized", "back_mutation", "framework_graft"
    mutations: list[str] = field(default_factory=list)  # e.g., ["VH:A45G", "VL:S32T"]
    humanness_score: Optional[float] = None
    notes: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "name": self.name,
            "parent_name": self.parent_name,
            "vh": self.vh,
            "vl": self.vl,
            "variant_type": self.variant_type,
            "mutations": self.mutations,
            "humanness_score": self.humanness_score,
            "notes": self.notes,
        }


class SequenceOptimizer:
    """Optimizer for generating variants of existing binders.

    This class handles:
    1. Loading starting sequences from YAML files
    2. CDR extraction using ANARCI
    3. Framework grafting onto human germlines
    4. Humanization with BioPhi/Sapiens
    5. Back-mutation variant generation
    """

    def __init__(self, data_dir: str = "data/starting_sequences"):
        """Initialize optimizer.

        Args:
            data_dir: Directory containing starting sequence YAML files.
        """
        self.data_dir = Path(data_dir)
        self._starting_sequences: dict[str, AntibodySequences] = {}

    def load_starting_sequence(self, name: str) -> AntibodySequences:
        """Load a starting sequence from YAML file.

        Args:
            name: Name of the sequence (e.g., "teplizumab", "sp34").

        Returns:
            AntibodySequences object.
        """
        if name in self._starting_sequences:
            return self._starting_sequences[name]

        yaml_path = self.data_dir / f"{name}.yaml"
        if not yaml_path.exists():
            raise FileNotFoundError(f"Starting sequence not found: {yaml_path}")

        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)

        # Handle different YAML structures
        if "humanized" in data:
            # Use humanized version if available
            seq_data = data["humanized"]
        elif "sequences" in data:
            seq_data = data["sequences"]
        else:
            seq_data = data

        sequences = AntibodySequences(
            name=name,
            vh=seq_data.get("vh", seq_data.get("VH", "")),
            vl=seq_data.get("vl", seq_data.get("VL")),
            source=data.get("source", "unknown"),
            notes=data.get("notes", ""),
        )

        self._starting_sequences[name] = sequences
        return sequences

    def load_all_starting_sequences(self) -> dict[str, AntibodySequences]:
        """Load all starting sequences from data directory.

        Returns:
            Dictionary mapping names to AntibodySequences.
        """
        for yaml_file in self.data_dir.glob("*.yaml"):
            if yaml_file.name == "README.md":
                continue
            name = yaml_file.stem
            if name != "affinity_mutations":  # Skip non-sequence files
                try:
                    self.load_starting_sequence(name)
                except Exception as e:
                    print(f"Warning: Could not load {name}: {e}")

        return self._starting_sequences

    def extract_cdrs(self, sequences: AntibodySequences) -> AntibodySequences:
        """Extract CDR sequences using ANARCI.

        Args:
            sequences: Antibody sequences to extract CDRs from.

        Returns:
            Updated AntibodySequences with CDR fields populated.
        """
        try:
            from src.analysis.numbering import number_sequence, extract_cdrs

            # Number and extract CDRs from VH
            vh_numbered = number_sequence(sequences.vh, chain_type="H")
            if vh_numbered:
                sequences.cdr_h1 = vh_numbered.cdr1
                sequences.cdr_h2 = vh_numbered.cdr2
                sequences.cdr_h3 = vh_numbered.cdr3

            # Number and extract CDRs from VL if present
            if sequences.vl:
                vl_numbered = number_sequence(sequences.vl, chain_type="L")
                if vl_numbered:
                    sequences.cdr_l1 = vl_numbered.cdr1
                    sequences.cdr_l2 = vl_numbered.cdr2
                    sequences.cdr_l3 = vl_numbered.cdr3

        except ImportError:
            print("Warning: ANARCI not available, CDR extraction skipped")
        except Exception as e:
            print(f"Warning: CDR extraction failed: {e}")

        return sequences

    def humanize(
        self,
        sequences: AntibodySequences,
        method: str = "sapiens",
    ) -> OptimizedVariant:
        """Humanize antibody sequences using BioPhi/Sapiens.

        Args:
            sequences: Antibody sequences to humanize.
            method: Humanization method ("sapiens" or "oasis").

        Returns:
            OptimizedVariant with humanized sequences.
        """
        try:
            from src.analysis.humanness import (
                get_humanization_suggestions,
                score_humanness_pair,
            )

            # Get humanization suggestions and apply them
            if method == "sapiens":
                vh_suggestions = get_humanization_suggestions(sequences.vh, "H")
                vh_humanized = self._apply_mutations(sequences.vh, vh_suggestions)
                vh_mutations = [f"{s['original']}{s['position']}{s['suggested']}" for s in vh_suggestions]

                vl_humanized = None
                vl_mutations = []

                if sequences.vl:
                    vl_suggestions = get_humanization_suggestions(sequences.vl, "L")
                    vl_humanized = self._apply_mutations(sequences.vl, vl_suggestions)
                    vl_mutations = [f"{s['original']}{s['position']}{s['suggested']}" for s in vl_suggestions]

                all_mutations = (
                    [f"VH:{m}" for m in vh_mutations] +
                    [f"VL:{m}" for m in vl_mutations]
                )

                # Score humanness of result
                humanness_report = score_humanness_pair(vh_humanized, vl_humanized)

                return OptimizedVariant(
                    name=f"{sequences.name}_humanized",
                    parent_name=sequences.name,
                    vh=vh_humanized,
                    vl=vl_humanized,
                    variant_type="humanized",
                    mutations=all_mutations,
                    humanness_score=humanness_report.mean_score,
                    notes=f"Humanized with {method}",
                )

            else:
                raise ValueError(f"Unknown humanization method: {method}")

        except ImportError:
            print("Warning: BioPhi not available, returning original sequences")
            return OptimizedVariant(
                name=f"{sequences.name}_humanized",
                parent_name=sequences.name,
                vh=sequences.vh,
                vl=sequences.vl,
                variant_type="humanized",
                mutations=[],
                notes="Humanization skipped (BioPhi not available)",
            )

    def _apply_mutations(self, sequence: str, suggestions: list[dict]) -> str:
        """Apply humanization mutation suggestions to a sequence.

        Args:
            sequence: Original amino acid sequence.
            suggestions: List of mutation dicts with 'position', 'original', 'suggested'.

        Returns:
            Mutated sequence.
        """
        if not suggestions:
            return sequence

        seq_list = list(sequence)
        for mut in suggestions:
            pos = mut.get("position")
            original = mut.get("original")
            suggested = mut.get("suggested")

            if pos is not None and 0 <= pos < len(seq_list):
                # Verify original matches (safety check)
                if seq_list[pos] == original:
                    seq_list[pos] = suggested

        return "".join(seq_list)

    def generate_back_mutations(
        self,
        humanized: OptimizedVariant,
        original: AntibodySequences,
        cdr_only: bool = True,
    ) -> list[OptimizedVariant]:
        """Generate back-mutation variants to restore original CDR residues.

        This is a safety measure in case humanization disrupts binding.

        Args:
            humanized: Humanized variant.
            original: Original sequences before humanization.
            cdr_only: If True, only consider mutations in CDR regions.

        Returns:
            List of back-mutation variants.
        """
        variants = []

        # Identify mutations that affected CDRs
        cdr_mutations = []
        for mutation in humanized.mutations:
            # Parse mutation format "VH:A45G" or "VL:S32T"
            chain, change = mutation.split(":")
            # For simplicity, include all VH mutations for back-mutation
            # In practice, would check if position is in CDR
            cdr_mutations.append(mutation)

        if not cdr_mutations:
            return variants

        # Generate variant with all CDR back-mutations
        back_vh = original.vh  # Restore original VH
        back_vl = original.vl  # Restore original VL

        # But keep framework humanization where possible
        # This is a simplified approach - full implementation would
        # selectively restore only CDR residues

        back_variant = OptimizedVariant(
            name=f"{humanized.name}_backmut",
            parent_name=humanized.name,
            vh=back_vh,
            vl=back_vl,
            variant_type="back_mutation",
            mutations=[f"RESTORE:{m}" for m in cdr_mutations],
            notes="Back-mutation variant restoring original CDR residues",
        )
        variants.append(back_variant)

        return variants

    def generate_all_variants(
        self,
        starting_names: Optional[list[str]] = None,
        include_humanization: bool = True,
        include_back_mutations: bool = True,
    ) -> list[OptimizedVariant]:
        """Generate all optimized variants from starting sequences.

        Args:
            starting_names: Names of starting sequences to optimize.
                If None, uses all available.
            include_humanization: If True, generate humanized variants.
            include_back_mutations: If True, generate back-mutation variants.

        Returns:
            List of all optimized variants.
        """
        if starting_names is None:
            self.load_all_starting_sequences()
            starting_names = list(self._starting_sequences.keys())

        all_variants = []

        for name in starting_names:
            try:
                sequences = self.load_starting_sequence(name)
                sequences = self.extract_cdrs(sequences)

                # Original as a variant
                original_variant = OptimizedVariant(
                    name=f"{name}_original",
                    parent_name=name,
                    vh=sequences.vh,
                    vl=sequences.vl,
                    variant_type="original",
                    notes=f"Original {name} sequences",
                )
                all_variants.append(original_variant)

                # Humanized variant
                if include_humanization:
                    humanized = self.humanize(sequences)
                    all_variants.append(humanized)

                    # Back-mutation variants
                    if include_back_mutations:
                        back_variants = self.generate_back_mutations(
                            humanized, sequences
                        )
                        all_variants.extend(back_variants)

            except Exception as e:
                print(f"Warning: Failed to optimize {name}: {e}")

        return all_variants


def optimize_existing_binders(
    binder_names: Optional[list[str]] = None,
    data_dir: str = "data/starting_sequences",
    include_humanization: bool = True,
    include_back_mutations: bool = True,
) -> list[OptimizedVariant]:
    """Convenience function for optimizing existing binders.

    Args:
        binder_names: Names of binders to optimize (e.g., ["teplizumab", "sp34"]).
            If None, uses all available.
        data_dir: Directory containing starting sequence files.
        include_humanization: If True, generate humanized variants.
        include_back_mutations: If True, generate back-mutation variants.

    Returns:
        List of optimized variants.
    """
    optimizer = SequenceOptimizer(data_dir)
    variants = optimizer.generate_all_variants(
        starting_names=binder_names,
        include_humanization=include_humanization,
        include_back_mutations=include_back_mutations,
    )

    print(f"Generated {len(variants)} variants from {len(binder_names or ['all'])} binders")

    return variants
