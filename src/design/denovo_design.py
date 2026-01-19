"""De novo binder design orchestration.

This module provides high-level functions for running de novo
VHH and scFv design against CD3 targets.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import json
import datetime

from src.design.boltzgen_runner import (
    BoltzGenRunner,
    BoltzGenConfig,
    BoltzGenDesign,
)


@dataclass
class DeNovoDesignConfig:
    """Configuration for de novo design campaign."""

    # Target structures
    target_structures: list[str] = field(default_factory=list)

    # Design parameters
    num_vhh_designs: int = 200
    num_scfv_designs: int = 200

    # Sampling
    seed: int = 42
    temperature: float = 1.0

    # Hotspot guidance (optional)
    # These are CD3ε residues to guide binding (e.g., OKT3 epitope)
    hotspot_residues: list[int] = field(default_factory=list)

    # Output
    output_dir: str = "data/outputs/denovo"

    # Execution
    use_modal: bool = True


@dataclass
class DeNovoDesignResult:
    """Result from a de novo design campaign."""

    vhh_designs: list[BoltzGenDesign]
    scfv_designs: list[BoltzGenDesign]
    config: DeNovoDesignConfig
    timestamp: str
    target_structures: list[str]

    @property
    def total_designs(self) -> int:
        """Total number of designs generated."""
        return len(self.vhh_designs) + len(self.scfv_designs)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "vhh_designs": [d.to_dict() for d in self.vhh_designs],
            "scfv_designs": [d.to_dict() for d in self.scfv_designs],
            "config": {
                "target_structures": self.config.target_structures,
                "num_vhh_designs": self.config.num_vhh_designs,
                "num_scfv_designs": self.config.num_scfv_designs,
                "seed": self.config.seed,
                "temperature": self.config.temperature,
                "hotspot_residues": self.config.hotspot_residues,
            },
            "timestamp": self.timestamp,
            "target_structures": self.target_structures,
            "total_designs": self.total_designs,
        }

    def save(self, output_path: Optional[str] = None) -> str:
        """Save results to JSON file."""
        if output_path is None:
            output_dir = Path(self.config.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / f"denovo_results_{self.timestamp}.json"

        output_path = Path(output_path)

        with open(output_path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

        return str(output_path)


class DeNovoDesigner:
    """Orchestrator for de novo binder design campaigns.

    This class manages the full de novo design workflow:
    1. Configure design parameters
    2. Run BoltzGen for VHH and scFv designs
    3. Collect and organize results
    4. Save outputs with provenance tracking
    """

    def __init__(self, config: Optional[DeNovoDesignConfig] = None):
        """Initialize designer.

        Args:
            config: Design configuration. If None, uses defaults.
        """
        self.config = config or DeNovoDesignConfig()

    def validate_targets(self) -> list[str]:
        """Validate target structures exist.

        Returns:
            List of error messages (empty if valid).
        """
        errors = []

        if not self.config.target_structures:
            errors.append("No target structures specified")
            return errors

        for target in self.config.target_structures:
            if not Path(target).exists():
                errors.append(f"Target structure not found: {target}")

        return errors

    def run_vhh_design(self) -> list[BoltzGenDesign]:
        """Run VHH design on all targets.

        Returns:
            List of VHH designs from all targets.
        """
        all_designs = []
        num_targets = len(self.config.target_structures)

        if num_targets == 0:
            return all_designs

        # Distribute designs evenly with remainder going to first targets
        base_per_target = self.config.num_vhh_designs // num_targets
        remainder = self.config.num_vhh_designs % num_targets

        for i, target in enumerate(self.config.target_structures):
            # First 'remainder' targets get one extra design
            num_designs = base_per_target + (1 if i < remainder else 0)

            if num_designs == 0:
                continue

            config = BoltzGenConfig(
                binder_type="vhh",
                num_designs=num_designs,
                target_pdb_path=target,
                hotspot_residues=self.config.hotspot_residues,
                seed=self.config.seed + i,
                temperature=self.config.temperature,
                output_dir=self.config.output_dir,
            )

            runner = BoltzGenRunner(config)
            designs = runner.run(use_modal=self.config.use_modal)

            # Update design IDs to include target info
            target_name = Path(target).stem
            for j, design in enumerate(designs):
                design.design_id = f"vhh_{target_name}_{j:04d}"

            all_designs.extend(designs)

        return all_designs

    def run_scfv_design(self) -> list[BoltzGenDesign]:
        """Run scFv design on all targets.

        Returns:
            List of scFv designs from all targets.
        """
        all_designs = []
        num_targets = len(self.config.target_structures)

        if num_targets == 0:
            return all_designs

        # Distribute designs evenly with remainder going to first targets
        base_per_target = self.config.num_scfv_designs // num_targets
        remainder = self.config.num_scfv_designs % num_targets

        for i, target in enumerate(self.config.target_structures):
            # First 'remainder' targets get one extra design
            num_designs = base_per_target + (1 if i < remainder else 0)

            if num_designs == 0:
                continue

            config = BoltzGenConfig(
                binder_type="scfv",
                num_designs=num_designs,
                target_pdb_path=target,
                hotspot_residues=self.config.hotspot_residues,
                seed=self.config.seed + 1000 + i,  # Different seed space for scFv
                temperature=self.config.temperature,
                output_dir=self.config.output_dir,
            )

            runner = BoltzGenRunner(config)
            designs = runner.run(use_modal=self.config.use_modal)

            # Update design IDs to include target info
            target_name = Path(target).stem
            for j, design in enumerate(designs):
                design.design_id = f"scfv_{target_name}_{j:04d}"

            all_designs.extend(designs)

        return all_designs

    def run(self) -> DeNovoDesignResult:
        """Run full de novo design campaign.

        Returns:
            DeNovoDesignResult containing all designs.

        Raises:
            ValueError: If configuration is invalid.
        """
        # Validate
        errors = self.validate_targets()
        if errors:
            raise ValueError(f"Invalid configuration: {'; '.join(errors)}")

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        # Run VHH design
        print(f"Running VHH design ({self.config.num_vhh_designs} designs)...")
        vhh_designs = self.run_vhh_design()
        print(f"  Generated {len(vhh_designs)} VHH designs")

        # Run scFv design
        print(f"Running scFv design ({self.config.num_scfv_designs} designs)...")
        scfv_designs = self.run_scfv_design()
        print(f"  Generated {len(scfv_designs)} scFv designs")

        result = DeNovoDesignResult(
            vhh_designs=vhh_designs,
            scfv_designs=scfv_designs,
            config=self.config,
            timestamp=timestamp,
            target_structures=self.config.target_structures,
        )

        print(f"Total: {result.total_designs} designs")

        return result


def run_denovo_design(
    target_structures: list[str],
    num_vhh: int = 200,
    num_scfv: int = 200,
    seed: int = 42,
    hotspot_residues: Optional[list[int]] = None,
    output_dir: str = "data/outputs/denovo",
    use_modal: bool = True,
) -> DeNovoDesignResult:
    """Convenience function for running de novo design.

    Args:
        target_structures: Paths to target PDB files.
        num_vhh: Number of VHH designs to generate.
        num_scfv: Number of scFv designs to generate.
        seed: Random seed for reproducibility.
        hotspot_residues: Optional CD3ε residues to guide binding.
        output_dir: Directory for output files.
        use_modal: If True, run on Modal.

    Returns:
        DeNovoDesignResult containing all designs.
    """
    config = DeNovoDesignConfig(
        target_structures=target_structures,
        num_vhh_designs=num_vhh,
        num_scfv_designs=num_scfv,
        seed=seed,
        hotspot_residues=hotspot_residues or [],
        output_dir=output_dir,
        use_modal=use_modal,
    )

    designer = DeNovoDesigner(config)
    result = designer.run()

    # Save results
    output_path = result.save()
    print(f"Results saved to: {output_path}")

    return result
