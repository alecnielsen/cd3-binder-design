"""Affinity variant generation for CD3 binders.

This module generates affinity-attenuated variants of CD3 binders
using literature-known mutations that reduce binding affinity.

The goal is to generate a panel spanning different affinities
(WT, ~10x weaker, ~100x weaker) to test the hypothesis that
~50 nM Kd balances efficacy with reduced CRS.

IMPORTANT: Mutation positions in the YAML file use IMGT numbering.
This module uses ANARCI to map IMGT positions to sequence indices
before applying mutations. This ensures mutations are applied to
the correct residue regardless of sequence length variations.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import warnings
import yaml
import itertools


@dataclass
class AffinityMutation:
    """A mutation known to affect CD3 binding affinity."""

    position: int  # IMGT position number (NOT raw sequence index)
    chain: str  # "VH" or "VL"
    wild_type: str  # Original residue at this IMGT position
    mutant: str  # Mutated residue
    effect: str  # "weaker" or "stronger"
    fold_change: float  # Approximate fold change in Kd (>1 = weaker)
    source: str  # Literature reference
    notes: str = ""

    @property
    def mutation_str(self) -> str:
        """String representation of mutation (e.g., 'VH:Y32A')."""
        return f"{self.chain}:{self.wild_type}{self.position}{self.mutant}"

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "position": self.position,
            "chain": self.chain,
            "wild_type": self.wild_type,
            "mutant": self.mutant,
            "effect": self.effect,
            "fold_change": self.fold_change,
            "source": self.source,
            "notes": self.notes,
        }


@dataclass
class AffinityVariant:
    """An affinity variant of a binder."""

    name: str
    parent_name: str
    vh: str
    vl: Optional[str] = None
    mutations: list[AffinityMutation] = field(default_factory=list)
    predicted_fold_change: float = 1.0  # Combined fold change
    affinity_class: str = "wild_type"  # "wild_type", "10x_weaker", "100x_weaker"
    notes: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "name": self.name,
            "parent_name": self.parent_name,
            "vh": self.vh,
            "vl": self.vl,
            "mutations": [m.to_dict() for m in self.mutations],
            "mutation_strings": [m.mutation_str for m in self.mutations],
            "predicted_fold_change": self.predicted_fold_change,
            "affinity_class": self.affinity_class,
            "notes": self.notes,
        }


class AffinityMutationLibrary:
    """Library of known affinity-affecting mutations for CD3 binders."""

    def __init__(self, data_path: str = "data/starting_sequences/affinity_mutations.yaml"):
        """Initialize library.

        Args:
            data_path: Path to affinity mutations YAML file.
        """
        self.data_path = Path(data_path)
        self._mutations: list[AffinityMutation] = []
        self._loaded = False

    def load(self) -> list[AffinityMutation]:
        """Load mutations from YAML file.

        Supports two YAML formats:
        1. Flat format: top-level 'mutations' list with wild_type/mutant keys
        2. Nested format: variants under 'teplizumab_variants' with from/to keys
        """
        if self._loaded:
            return self._mutations

        if not self.data_path.exists():
            print(f"Warning: Affinity mutations file not found: {self.data_path}")
            return []

        with open(self.data_path, "r") as f:
            data = yaml.safe_load(f)

        # Try flat format first (top-level 'mutations' list)
        mutations_data = data.get("mutations", [])

        # If no flat mutations, parse nested teplizumab_variants structure
        if not mutations_data and "teplizumab_variants" in data:
            seen_mutations = set()  # Deduplicate mutations across variants
            for variant_key, variant_data in data["teplizumab_variants"].items():
                if not isinstance(variant_data, dict):
                    continue
                variant_mutations = variant_data.get("mutations", [])
                fold_change = variant_data.get("kd_fold_change", 10.0)
                variant_notes = variant_data.get("notes", "")

                for m in variant_mutations:
                    # Create a unique key for deduplication
                    mutation_key = (m.get("chain"), m.get("position"), m.get("to"))
                    if mutation_key in seen_mutations:
                        continue
                    seen_mutations.add(mutation_key)

                    # Map from/to to wild_type/mutant for compatibility
                    mutations_data.append({
                        "position": m.get("position"),
                        "chain": m.get("chain"),
                        "wild_type": m.get("from") or m.get("wild_type"),
                        "mutant": m.get("to") or m.get("mutant"),
                        "effect": m.get("effect", "weaker"),
                        "fold_change": m.get("fold_change", fold_change),
                        "source": m.get("source", "literature"),
                        "notes": m.get("notes", variant_notes),
                        "cdr": m.get("cdr", ""),
                    })

        for m in mutations_data:
            mutation = AffinityMutation(
                position=m["position"],
                chain=m["chain"],
                wild_type=m["wild_type"],
                mutant=m["mutant"],
                effect=m.get("effect", "weaker"),
                fold_change=m.get("fold_change", 10.0),
                source=m.get("source", "unknown"),
                notes=m.get("notes", ""),
            )
            self._mutations.append(mutation)

        self._loaded = True
        return self._mutations

    def get_mutations_for_effect(
        self,
        target_fold_change: float,
        tolerance: float = 0.5,
    ) -> list[AffinityMutation]:
        """Get mutations that achieve approximately the target fold change.

        Args:
            target_fold_change: Desired fold change in Kd.
            tolerance: Allowed deviation (log scale factor).

        Returns:
            List of mutations near target fold change.
        """
        self.load()

        import math
        log_target = math.log10(target_fold_change)

        matching = []
        for m in self._mutations:
            log_fc = math.log10(m.fold_change)
            if abs(log_fc - log_target) <= tolerance:
                matching.append(m)

        return matching

    def get_weaker_mutations(self, min_fold_change: float = 5.0) -> list[AffinityMutation]:
        """Get mutations that weaken binding.

        Args:
            min_fold_change: Minimum fold change to include.

        Returns:
            List of weakening mutations.
        """
        self.load()
        return [m for m in self._mutations if m.effect == "weaker" and m.fold_change >= min_fold_change]


class AffinityVariantGenerator:
    """Generator for affinity variant panels.

    Creates variants at different affinity levels:
    - Wild-type (no mutations)
    - ~10x weaker
    - ~100x weaker

    Uses ANARCI to map IMGT positions to sequence indices, ensuring
    mutations are applied to the correct residue regardless of
    sequence length variations between different antibodies.
    """

    def __init__(
        self,
        mutation_library: Optional[AffinityMutationLibrary] = None,
    ):
        """Initialize generator.

        Args:
            mutation_library: Library of known mutations.
        """
        self.library = mutation_library or AffinityMutationLibrary()
        # Cache for IMGT→sequence index mappings
        self._imgt_mapping_cache: dict[str, dict[int, int]] = {}

    def _get_imgt_mapping(self, sequence: str, chain_type: str) -> dict[int, int]:
        """Get or create IMGT position to sequence index mapping.

        Args:
            sequence: Amino acid sequence.
            chain_type: "H" for VH, "L" for VL.

        Returns:
            Dict mapping IMGT position (int) to 0-indexed sequence position.
        """
        cache_key = f"{chain_type}:{sequence}"
        if cache_key in self._imgt_mapping_cache:
            return self._imgt_mapping_cache[cache_key]

        try:
            from src.analysis.numbering import get_imgt_to_sequence_mapping
            mapping = get_imgt_to_sequence_mapping(sequence, chain_type)
        except ImportError as e:
            # ANARCI is REQUIRED for mutation application since mutations are
            # specified in IMGT numbering. Fallback to raw indexing would apply
            # mutations to incorrect residues.
            raise ImportError(
                "ANARCI is required for affinity mutation application. "
                "Mutation positions use IMGT numbering which requires ANARCI to map "
                "to sequence positions. Install with: pip install anarci\n"
                f"Original error: {e}"
            ) from e

        self._imgt_mapping_cache[cache_key] = mapping
        return mapping

    def apply_mutation(
        self,
        sequence: str,
        mutation: AffinityMutation,
        imgt_mapping: Optional[dict[int, int]] = None,
    ) -> tuple[str, bool, str]:
        """Apply a mutation to a sequence using IMGT position mapping.

        Args:
            sequence: Amino acid sequence.
            mutation: Mutation to apply (position is IMGT numbered).
            imgt_mapping: Pre-computed IMGT→sequence index mapping.
                If None, will be computed (slower).

        Returns:
            Tuple of (mutated_sequence, success, error_message).
            error_message is empty string on success.
        """
        # Get or compute IMGT mapping
        if imgt_mapping is None:
            chain_type = "H" if mutation.chain == "VH" else "L"
            imgt_mapping = self._get_imgt_mapping(sequence, chain_type)

        # Map IMGT position to sequence index
        imgt_pos = mutation.position
        seq_idx = imgt_mapping.get(imgt_pos)

        if seq_idx is None:
            return sequence, False, f"IMGT position {imgt_pos} not found in sequence"

        if seq_idx < 0 or seq_idx >= len(sequence):
            return sequence, False, f"Mapped index {seq_idx} out of bounds"

        actual_residue = sequence[seq_idx]
        if actual_residue != mutation.wild_type:
            return sequence, False, (
                f"Residue mismatch at IMGT {imgt_pos} (seq idx {seq_idx}): "
                f"expected {mutation.wild_type}, found {actual_residue}"
            )

        mutated = sequence[:seq_idx] + mutation.mutant + sequence[seq_idx + 1:]
        return mutated, True, ""

    def apply_mutation_legacy(
        self,
        sequence: str,
        mutation: AffinityMutation,
    ) -> tuple[str, bool]:
        """Apply a mutation using direct position indexing (legacy behavior).

        WARNING: This treats the mutation position as a raw 1-indexed sequence
        position, NOT IMGT numbering. Only use for sequences where positions
        were specified as sequence indices.

        Args:
            sequence: Amino acid sequence.
            mutation: Mutation to apply.

        Returns:
            Tuple of (mutated_sequence, success).
        """
        pos = mutation.position - 1  # Convert to 0-indexed

        if pos < 0 or pos >= len(sequence):
            return sequence, False

        if sequence[pos] != mutation.wild_type:
            return sequence, False

        mutated = sequence[:pos] + mutation.mutant + sequence[pos + 1:]
        return mutated, True

    def generate_single_mutants(
        self,
        parent_name: str,
        vh: str,
        vl: Optional[str] = None,
    ) -> list[AffinityVariant]:
        """Generate single-mutation variants.

        Args:
            parent_name: Name of parent sequence.
            vh: VH sequence.
            vl: VL sequence (optional).

        Returns:
            List of single-mutation variants.
        """
        variants = []
        mutations = self.library.load()

        # Pre-compute IMGT mappings for efficiency
        vh_mapping = self._get_imgt_mapping(vh, "H")
        vl_mapping = self._get_imgt_mapping(vl, "L") if vl else {}

        for mutation in mutations:
            if mutation.chain == "VH":
                mutated_vh, success, error = self.apply_mutation(vh, mutation, vh_mapping)
                if success:
                    variant = AffinityVariant(
                        name=f"{parent_name}_{mutation.mutation_str}",
                        parent_name=parent_name,
                        vh=mutated_vh,
                        vl=vl,
                        mutations=[mutation],
                        predicted_fold_change=mutation.fold_change,
                        affinity_class=self._classify_affinity(mutation.fold_change),
                        notes=f"Single mutation: {mutation.mutation_str}",
                    )
                    variants.append(variant)
                elif error:
                    warnings.warn(f"Skipping mutation {mutation.mutation_str}: {error}")

            elif mutation.chain == "VL" and vl:
                mutated_vl, success, error = self.apply_mutation(vl, mutation, vl_mapping)
                if success:
                    variant = AffinityVariant(
                        name=f"{parent_name}_{mutation.mutation_str}",
                        parent_name=parent_name,
                        vh=vh,
                        vl=mutated_vl,
                        mutations=[mutation],
                        predicted_fold_change=mutation.fold_change,
                        affinity_class=self._classify_affinity(mutation.fold_change),
                        notes=f"Single mutation: {mutation.mutation_str}",
                    )
                    variants.append(variant)
                elif error:
                    warnings.warn(f"Skipping mutation {mutation.mutation_str}: {error}")

        return variants

    def generate_combination_mutants(
        self,
        parent_name: str,
        vh: str,
        vl: Optional[str] = None,
        max_mutations: int = 2,
    ) -> list[AffinityVariant]:
        """Generate combinatorial mutation variants.

        Args:
            parent_name: Name of parent sequence.
            vh: VH sequence.
            vl: VL sequence (optional).
            max_mutations: Maximum number of mutations to combine.

        Returns:
            List of combinatorial variants.
        """
        variants = []
        mutations = self.library.load()

        # Pre-compute IMGT mapping for VH
        vh_mapping = self._get_imgt_mapping(vh, "H")

        # Group mutations by chain
        vh_mutations = [m for m in mutations if m.chain == "VH"]
        vl_mutations = [m for m in mutations if m.chain == "VL"]

        # Generate VH combinations
        for r in range(2, min(max_mutations + 1, len(vh_mutations) + 1)):
            for combo in itertools.combinations(vh_mutations, r):
                # Check mutations don't overlap
                positions = [m.position for m in combo]
                if len(positions) != len(set(positions)):
                    continue

                # Apply all mutations sequentially
                # Note: We re-compute mapping after each mutation since indices shift
                mutated_vh = vh
                current_mapping = vh_mapping.copy()
                all_success = True
                for mutation in combo:
                    mutated_vh, success, error = self.apply_mutation(
                        mutated_vh, mutation, current_mapping
                    )
                    if not success:
                        all_success = False
                        break
                    # Note: For same-length substitutions, mapping doesn't change

                if all_success:
                    # Estimate combined fold change (multiplicative assumption)
                    # WARNING: This assumes mutations are independent, which often fails
                    # due to epistasis. CDR mutations can be synergistic or antagonistic.
                    # Actual fold changes may differ significantly from predictions.
                    combined_fc = 1.0
                    for m in combo:
                        combined_fc *= m.fold_change

                    mutation_strs = [m.mutation_str for m in combo]
                    epistasis_warning = (
                        "EPISTASIS WARNING: Fold change is estimated assuming "
                        "mutation independence. Actual affinity may differ significantly "
                        "due to synergistic or antagonistic interactions between mutations."
                    )
                    variant = AffinityVariant(
                        name=f"{parent_name}_{'_'.join([m.mutant + str(m.position) for m in combo])}",
                        parent_name=parent_name,
                        vh=mutated_vh,
                        vl=vl,
                        mutations=list(combo),
                        predicted_fold_change=combined_fc,
                        affinity_class=self._classify_affinity(combined_fc),
                        notes=f"Combination: {', '.join(mutation_strs)}. {epistasis_warning}",
                    )
                    variants.append(variant)

        return variants

    def _classify_affinity(self, fold_change: float) -> str:
        """Classify affinity based on fold change.

        Args:
            fold_change: Fold change in Kd (>1 = weaker).

        Returns:
            Affinity class string.
        """
        if fold_change < 3:
            return "wild_type"
        elif fold_change < 30:
            return "10x_weaker"
        else:
            return "100x_weaker"

    def generate_affinity_panel(
        self,
        parent_name: str,
        vh: str,
        vl: Optional[str] = None,
        include_wild_type: bool = True,
        include_singles: bool = True,
        include_combinations: bool = True,
        target_classes: Optional[list[str]] = None,
    ) -> list[AffinityVariant]:
        """Generate a panel of affinity variants.

        Args:
            parent_name: Name of parent sequence.
            vh: VH sequence.
            vl: VL sequence (optional).
            include_wild_type: Include unmodified sequence.
            include_singles: Include single mutations.
            include_combinations: Include combination mutations.
            target_classes: If specified, only include these affinity classes.

        Returns:
            List of affinity variants covering different affinity levels.
        """
        panel = []

        # Wild-type
        if include_wild_type:
            wt = AffinityVariant(
                name=f"{parent_name}_WT",
                parent_name=parent_name,
                vh=vh,
                vl=vl,
                mutations=[],
                predicted_fold_change=1.0,
                affinity_class="wild_type",
                notes="Wild-type (no mutations)",
            )
            panel.append(wt)

        # Single mutants
        if include_singles:
            singles = self.generate_single_mutants(parent_name, vh, vl)
            panel.extend(singles)

        # Combination mutants
        if include_combinations:
            combos = self.generate_combination_mutants(parent_name, vh, vl)
            panel.extend(combos)

        # Filter by target classes if specified
        if target_classes:
            panel = [v for v in panel if v.affinity_class in target_classes]

        return panel


def generate_affinity_variants(
    parent_name: str,
    vh: str,
    vl: Optional[str] = None,
    mutations_file: str = "data/starting_sequences/affinity_mutations.yaml",
    target_classes: Optional[list[str]] = None,
) -> list[AffinityVariant]:
    """Convenience function for generating affinity variants.

    Args:
        parent_name: Name of parent sequence.
        vh: VH sequence.
        vl: VL sequence (optional).
        mutations_file: Path to mutations YAML file.
        target_classes: Target affinity classes to include.

    Returns:
        List of affinity variants.
    """
    library = AffinityMutationLibrary(mutations_file)
    generator = AffinityVariantGenerator(library)

    panel = generator.generate_affinity_panel(
        parent_name=parent_name,
        vh=vh,
        vl=vl,
        target_classes=target_classes,
    )

    # Summarize
    classes = {}
    for v in panel:
        classes[v.affinity_class] = classes.get(v.affinity_class, 0) + 1

    print(f"Generated {len(panel)} affinity variants:")
    for cls, count in sorted(classes.items()):
        print(f"  {cls}: {count}")

    return panel
