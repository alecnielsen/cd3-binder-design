"""Base classes and utilities for bispecific antibody formatting."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional
from pathlib import Path
import yaml


@dataclass
class AntibodyChain:
    """A single antibody chain (heavy or light)."""

    name: str
    sequence: str
    chain_type: str  # "heavy", "light", "vhh", "scfv"
    components: list[str] = field(default_factory=list)  # e.g., ["VH", "CH1", "hinge", "CH2", "CH3"]

    def __len__(self) -> int:
        return len(self.sequence)


@dataclass
class BispecificConstruct:
    """A complete bispecific antibody construct."""

    name: str
    format_type: str  # "crossmab", "fab_scfv", "fab_vhh", "igg_scfv", "igg_vhh"
    chains: list[AntibodyChain]
    target_1: str  # e.g., "HER2"
    target_2: str  # e.g., "CD3"
    notes: str = ""

    def get_chain(self, name: str) -> Optional[AntibodyChain]:
        """Get chain by name."""
        for chain in self.chains:
            if chain.name == name:
                return chain
        return None

    def to_fasta(self) -> str:
        """Convert to FASTA format."""
        lines = []
        for chain in self.chains:
            lines.append(f">{self.name}_{chain.name}")
            # Wrap sequence at 80 characters
            seq = chain.sequence
            for i in range(0, len(seq), 80):
                lines.append(seq[i:i+80])
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "format_type": self.format_type,
            "target_1": self.target_1,
            "target_2": self.target_2,
            "chains": [
                {
                    "name": c.name,
                    "sequence": c.sequence,
                    "chain_type": c.chain_type,
                    "components": c.components,
                    "length": len(c.sequence),
                }
                for c in self.chains
            ],
            "notes": self.notes,
        }


class SequenceLibrary:
    """Library of constant region and linker sequences."""

    def __init__(self, data_dir: Optional[Path] = None):
        """Initialize sequence library.

        Args:
            data_dir: Directory containing YAML sequence files.
        """
        self.data_dir = data_dir or Path("data/frameworks")
        self._constant_regions = {}
        self._linkers = {}
        self._knob_hole = {}
        self._loaded = False

    def _load(self):
        """Load sequences from YAML files."""
        if self._loaded:
            return

        # Load constant regions
        const_file = self.data_dir / "igg1_constant_regions.yaml"
        if const_file.exists():
            with open(const_file) as f:
                data = yaml.safe_load(f)
                self._constant_regions = data

        # Load linkers
        linker_file = self.data_dir / "linkers.yaml"
        if linker_file.exists():
            with open(linker_file) as f:
                data = yaml.safe_load(f)
                self._linkers = data

        # Load knob-hole mutations
        kh_file = self.data_dir / "knob_hole_mutations.yaml"
        if kh_file.exists():
            with open(kh_file) as f:
                data = yaml.safe_load(f)
                self._knob_hole = data

        self._loaded = True

    def get_ch1(self) -> str:
        """Get CH1 sequence."""
        self._load()
        return self._constant_regions.get("ch1", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_hinge(self) -> str:
        """Get hinge sequence."""
        self._load()
        return self._constant_regions.get("hinge", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_ch2(self) -> str:
        """Get CH2 sequence."""
        self._load()
        return self._constant_regions.get("ch2", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_ch3_knob(self) -> str:
        """Get CH3 with knob mutation."""
        self._load()
        return self._knob_hole.get("complete_sequences", {}).get("fc_knob", {}).get("ch3", "").replace("\n", "").replace(" ", "")

    def get_ch3_hole(self) -> str:
        """Get CH3 with hole mutations."""
        self._load()
        return self._knob_hole.get("complete_sequences", {}).get("fc_hole", {}).get("ch3", "").replace("\n", "").replace(" ", "")

    def get_ch3_standard(self) -> str:
        """Get standard CH3 (no knob-hole)."""
        self._load()
        return self._constant_regions.get("ch3", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_kappa_cl(self) -> str:
        """Get kappa light chain constant region."""
        self._load()
        return self._constant_regions.get("kappa_cl", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_lambda_cl(self) -> str:
        """Get lambda light chain constant region."""
        self._load()
        return self._constant_regions.get("lambda_cl", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_fc_lala(self) -> str:
        """Get Fc with LALA silencing mutations."""
        self._load()
        return self._constant_regions.get("fc_lala", {}).get("sequence", "").replace("\n", "").replace(" ", "")

    def get_scfv_linker(self) -> str:
        """Get standard scFv linker (G4S)3."""
        self._load()
        return self._linkers.get("scfv_linkers", {}).get("standard_15aa", {}).get("sequence", "GGGGSGGGGSGGGGS")

    def get_fc_fusion_linker(self) -> str:
        """Get Fc fusion linker (G4S)2."""
        self._load()
        return self._linkers.get("fc_fusion_linkers", {}).get("standard_10aa", {}).get("sequence", "GGGGSGGGGS")


class BispecificFormatter(ABC):
    """Abstract base class for bispecific antibody formatters."""

    def __init__(self, sequence_library: Optional[SequenceLibrary] = None):
        """Initialize formatter.

        Args:
            sequence_library: Library of constant region sequences.
        """
        self.library = sequence_library or SequenceLibrary()

    @property
    @abstractmethod
    def format_name(self) -> str:
        """Name of this bispecific format."""
        pass

    @abstractmethod
    def assemble(
        self,
        target_vh: str,
        target_vl: str,
        cd3_binder: str,  # VHH or scFv sequence
        cd3_binder_vl: Optional[str] = None,  # For scFv, the VL portion
        name: str = "bispecific",
        target_name: str = "HER2",
    ) -> BispecificConstruct:
        """Assemble bispecific construct.

        Args:
            target_vh: VH sequence for tumor target arm.
            target_vl: VL sequence for tumor target arm.
            cd3_binder: VHH sequence or VH portion of scFv for CD3 arm.
            cd3_binder_vl: VL portion for scFv CD3 arm (None for VHH).
            name: Name for the construct.
            target_name: Name of tumor target.

        Returns:
            BispecificConstruct with all chains.
        """
        pass

    def make_scfv(self, vh: str, vl: str, linker: Optional[str] = None) -> str:
        """Create scFv from VH and VL.

        Args:
            vh: VH sequence.
            vl: VL sequence.
            linker: Linker sequence (default: (G4S)3).

        Returns:
            scFv sequence (VH-linker-VL).
        """
        linker = linker or self.library.get_scfv_linker()
        return vh + linker + vl
