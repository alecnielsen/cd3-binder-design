"""BoltzGen runner for de novo binder design.

This module provides the interface for running BoltzGen on Modal
for GPU-accelerated binder design.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import json
import hashlib


@dataclass
class BoltzGenConfig:
    """Configuration for BoltzGen design runs."""

    # Design parameters
    binder_type: str = "vhh"  # "vhh" or "scfv"
    num_designs: int = 100

    # Target specification
    target_pdb_path: Optional[str] = None
    target_chain: str = "A"  # Chain ID for CD3ε

    # Hotspot residues (optional - guide design to specific epitope)
    hotspot_residues: list[int] = field(default_factory=list)

    # Sampling parameters
    seed: int = 42
    temperature: float = 1.0
    num_recycles: int = 3

    # Output
    output_dir: str = "data/outputs/denovo"

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "binder_type": self.binder_type,
            "num_designs": self.num_designs,
            "target_pdb_path": self.target_pdb_path,
            "target_chain": self.target_chain,
            "hotspot_residues": self.hotspot_residues,
            "seed": self.seed,
            "temperature": self.temperature,
            "num_recycles": self.num_recycles,
            "output_dir": self.output_dir,
        }

    def config_hash(self) -> str:
        """Generate hash of config for reproducibility tracking."""
        config_str = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.sha256(config_str.encode()).hexdigest()[:12]


@dataclass
class BoltzGenDesign:
    """Result from a single BoltzGen design."""

    sequence: str
    confidence: float
    design_id: str
    binder_type: str
    target_structure: str

    # Optional structural information
    plddt: Optional[float] = None
    ptm: Optional[float] = None

    # Metadata
    seed: Optional[int] = None
    temperature: Optional[float] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "sequence": self.sequence,
            "confidence": self.confidence,
            "design_id": self.design_id,
            "binder_type": self.binder_type,
            "target_structure": self.target_structure,
            "plddt": self.plddt,
            "ptm": self.ptm,
            "seed": self.seed,
            "temperature": self.temperature,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "BoltzGenDesign":
        """Create from dictionary."""
        return cls(**data)


class BoltzGenRunner:
    """Runner for BoltzGen de novo binder design.

    This class manages the interface to BoltzGen, which runs on Modal
    for GPU acceleration. It handles:
    - Configuration validation
    - Job submission to Modal
    - Result collection and parsing
    - Reproducibility tracking
    """

    def __init__(self, config: Optional[BoltzGenConfig] = None):
        """Initialize runner with optional config.

        Args:
            config: BoltzGen configuration. If None, uses defaults.
        """
        self.config = config or BoltzGenConfig()
        self._modal_available = None

    def check_modal_available(self) -> bool:
        """Check if Modal is available and configured."""
        if self._modal_available is not None:
            return self._modal_available

        try:
            import modal
            # Just check that modal is importable (auth is checked on first call)
            self._modal_available = True
        except ImportError:
            self._modal_available = False

        return self._modal_available

    def validate_config(self) -> list[str]:
        """Validate configuration and return list of errors."""
        errors = []

        if self.config.binder_type not in ["vhh", "scfv"]:
            errors.append(f"Invalid binder_type: {self.config.binder_type}. Must be 'vhh' or 'scfv'.")

        if self.config.num_designs < 1:
            errors.append(f"num_designs must be >= 1, got {self.config.num_designs}")

        if self.config.target_pdb_path:
            target_path = Path(self.config.target_pdb_path)
            if not target_path.exists():
                errors.append(f"Target PDB not found: {self.config.target_pdb_path}")
        else:
            errors.append("target_pdb_path is required")

        if self.config.temperature <= 0:
            errors.append(f"temperature must be > 0, got {self.config.temperature}")

        return errors

    def run_local_mock(self) -> list[BoltzGenDesign]:
        """Run a mock design for testing without GPU.

        Returns mock designs for testing the pipeline without
        actually running BoltzGen (which requires GPU).
        """
        import random

        random.seed(self.config.seed)

        designs = []
        for i in range(min(self.config.num_designs, 10)):  # Limit mock to 10
            if self.config.binder_type == "vhh":
                # Generate mock VHH sequence (~120 aa)
                seq_len = random.randint(115, 125)
            else:
                # Generate mock scFv sequence (~250 aa)
                seq_len = random.randint(240, 260)

            # Generate random sequence (not realistic, just for testing)
            aa = "ACDEFGHIKLMNPQRSTVWY"
            sequence = "".join(random.choices(aa, k=seq_len))

            design = BoltzGenDesign(
                sequence=sequence,
                confidence=random.uniform(0.5, 0.95),
                design_id=f"mock_{self.config.binder_type}_{i:04d}",
                binder_type=self.config.binder_type,
                target_structure=self.config.target_pdb_path or "unknown",
                plddt=random.uniform(70, 95),
                ptm=random.uniform(0.6, 0.9),
                seed=self.config.seed,
                temperature=self.config.temperature,
            )
            designs.append(design)

        return designs

    def run(self, use_modal: bool = True) -> list[BoltzGenDesign]:
        """Run BoltzGen design.

        Args:
            use_modal: If True, run on Modal. If False, run local mock.

        Returns:
            List of BoltzGenDesign results.

        Raises:
            ValueError: If configuration is invalid.
            RuntimeError: If Modal is not available when requested.
        """
        # Validate config
        errors = self.validate_config()
        if errors:
            raise ValueError(f"Invalid configuration: {'; '.join(errors)}")

        if not use_modal:
            return self.run_local_mock()

        if not self.check_modal_available():
            raise RuntimeError(
                "Modal is not available. Install with 'pip install modal' "
                "and run 'modal setup' to authenticate."
            )

        # Import and run Modal function
        return self._run_modal()

    def _run_modal(self) -> list[BoltzGenDesign]:
        """Run BoltzGen on Modal.

        This method calls the deployed Modal function.
        """
        # Import the Modal app
        try:
            import modal

            # Get the deployed function (Modal API v0.60+)
            run_boltzgen = modal.Function.from_name("boltzgen-cd3", "run_boltzgen")

            # Read target PDB
            with open(self.config.target_pdb_path, "r") as f:
                target_pdb_content = f.read()

            # Map binder_type to protocol and binder_length
            if self.config.binder_type == "vhh":
                protocol = "nanobody-anything"
                binder_length = 120  # VHH ~120 aa
            else:  # scfv
                protocol = "nanobody-anything"  # BoltzGen uses same protocol
                binder_length = 250  # scFv ~250 aa

            # Call Modal function
            results = run_boltzgen.remote(
                target_pdb_content=target_pdb_content,
                target_chain=self.config.target_chain,
                num_designs=self.config.num_designs,
                binder_length=binder_length,
                hotspot_residues=self.config.hotspot_residues,
                protocol=protocol,
                seed=self.config.seed,
            )

            # Parse results
            designs = []
            for i, result in enumerate(results):
                design = BoltzGenDesign(
                    sequence=result["sequence"],
                    confidence=result.get("confidence", 0.0),
                    design_id=f"boltzgen_{self.config.binder_type}_{i:04d}",
                    binder_type=self.config.binder_type,
                    target_structure=self.config.target_pdb_path,
                    plddt=result.get("plddt"),
                    ptm=result.get("ptm"),
                    seed=self.config.seed,
                    temperature=self.config.temperature,
                )
                designs.append(design)

            return designs

        except Exception as e:
            raise RuntimeError(f"Failed to run BoltzGen on Modal: {e}")

    def save_designs(
        self,
        designs: list[BoltzGenDesign],
        output_path: Optional[str] = None,
    ) -> str:
        """Save designs to JSON file.

        Args:
            designs: List of designs to save.
            output_path: Path to save to. If None, uses config output_dir.

        Returns:
            Path to saved file.
        """
        import datetime

        if output_path is None:
            output_dir = Path(self.config.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = output_dir / f"boltzgen_designs_{timestamp}.json"

        output_path = Path(output_path)

        output_data = {
            "config": self.config.to_dict(),
            "config_hash": self.config.config_hash(),
            "num_designs": len(designs),
            "designs": [d.to_dict() for d in designs],
        }

        with open(output_path, "w") as f:
            json.dump(output_data, f, indent=2)

        return str(output_path)


def run_boltzgen_batch(
    target_pdbs: list[str],
    binder_type: str = "vhh",
    num_designs_per_target: int = 100,
    seed: int = 42,
    use_modal: bool = True,
) -> dict[str, list[BoltzGenDesign]]:
    """Run BoltzGen on multiple target structures.

    Args:
        target_pdbs: List of paths to target PDB files.
        binder_type: Type of binder to design ("vhh" or "scfv").
        num_designs_per_target: Number of designs per target.
        seed: Random seed for reproducibility.
        use_modal: If True, run on Modal.

    Returns:
        Dictionary mapping target paths to lists of designs.
    """
    all_designs = {}

    for i, target_pdb in enumerate(target_pdbs):
        config = BoltzGenConfig(
            binder_type=binder_type,
            num_designs=num_designs_per_target,
            target_pdb_path=target_pdb,
            seed=seed + i,  # Different seed per target
        )

        runner = BoltzGenRunner(config)
        designs = runner.run(use_modal=use_modal)
        all_designs[target_pdb] = designs

    return all_designs
