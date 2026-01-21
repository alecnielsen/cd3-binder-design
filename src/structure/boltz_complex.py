"""Boltz-2 wrapper for complex structure prediction.

Boltz-2 predicts protein complex structures and provides
confidence scores (pDockQ, pTM) for binding assessment.

IMPORTANT: pDockQ is a STRUCTURAL CONFIDENCE score, not an
affinity predictor. High pDockQ means the model is confident
in the predicted structure, not that binding is strong.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import json
import tempfile


@dataclass
class ComplexPredictionResult:
    """Result from Boltz-2 complex prediction."""

    pdb_string: str
    binder_sequence: str
    target_sequence: str

    # Confidence scores
    pdockq: float  # Predicted DockQ - structural confidence
    ptm: float  # Predicted TM-score
    plddt_mean: float  # Mean pLDDT
    ipae: float  # Interface PAE (lower = better)

    # Interface analysis
    interface_residues_binder: list[int] = field(default_factory=list)
    interface_residues_target: list[int] = field(default_factory=list)
    num_contacts: int = 0
    interface_area: float = 0.0  # Buried surface area in Å²

    # Metadata
    model_version: str = "boltz2"
    seed: Optional[int] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "binder_sequence": self.binder_sequence,
            "target_sequence": self.target_sequence,
            "pdockq": self.pdockq,
            "pdockq_note": "Structural confidence, NOT affinity predictor",
            "ptm": self.ptm,
            "plddt_mean": self.plddt_mean,
            "ipae": self.ipae,
            "interface_residues_binder": self.interface_residues_binder,
            "interface_residues_target": self.interface_residues_target,
            "num_contacts": self.num_contacts,
            "interface_area": self.interface_area,
            "model_version": self.model_version,
            "seed": self.seed,
        }

    def save_pdb(self, output_path: str) -> str:
        """Save structure to PDB file."""
        with open(output_path, "w") as f:
            f.write(self.pdb_string)
        return output_path

    def passes_threshold(
        self,
        min_pdockq: float = 0.5,
        min_interface_area: float = 800.0,
        min_contacts: int = 10,
    ) -> tuple[bool, list[str]]:
        """Check if result passes quality thresholds.

        Args:
            min_pdockq: Minimum pDockQ score.
            min_interface_area: Minimum interface area in Å².
            min_contacts: Minimum number of contacts.

        Returns:
            Tuple of (passes, list of failure reasons).
        """
        failures = []

        if self.pdockq < min_pdockq:
            failures.append(f"pDockQ {self.pdockq:.3f} < {min_pdockq}")

        if self.interface_area < min_interface_area:
            failures.append(f"Interface area {self.interface_area:.1f} < {min_interface_area}")

        if self.num_contacts < min_contacts:
            failures.append(f"Contacts {self.num_contacts} < {min_contacts}")

        return len(failures) == 0, failures


class Boltz2Predictor:
    """Wrapper for Boltz-2 complex structure prediction.

    Boltz-2 predicts protein-protein complex structures and
    provides confidence metrics for assessing binding.

    NOTE: This wrapper can run locally (if Boltz-2 is installed)
    or via Modal for GPU acceleration.
    """

    def __init__(self, use_modal: bool = True):
        """Initialize predictor.

        Args:
            use_modal: If True, use Modal for GPU compute.
        """
        self.use_modal = use_modal
        self._model = None
        self._available = None

    def is_available(self) -> bool:
        """Check if Boltz-2 is available locally."""
        if self._available is not None:
            return self._available

        try:
            import boltz
            self._available = True
        except ImportError:
            self._available = False

        return self._available

    def predict_complex(
        self,
        binder_sequence: str,
        target_pdb_path: str,
        target_chain: str = "A",
        seed: int = 42,
    ) -> ComplexPredictionResult:
        """Predict binder-target complex structure.

        Args:
            binder_sequence: Binder amino acid sequence.
            target_pdb_path: Path to target structure PDB.
            target_chain: Chain ID in target PDB.
            seed: Random seed for reproducibility.

        Returns:
            ComplexPredictionResult with predicted complex.
        """
        if self.use_modal:
            return self._predict_modal(binder_sequence, target_pdb_path, target_chain, seed)
        else:
            return self._predict_local(binder_sequence, target_pdb_path, target_chain, seed)

    def _predict_local(
        self,
        binder_sequence: str,
        target_pdb_path: str,
        target_chain: str,
        seed: int,
    ) -> ComplexPredictionResult:
        """Run prediction locally."""
        if not self.is_available():
            raise RuntimeError(
                "Boltz-2 not available locally. Install with: pip install boltz "
                "or use use_modal=True for GPU compute."
            )

        # Extract target sequence from PDB
        from src.structure.pdb_utils import extract_sequence_from_pdb
        target_sequence = extract_sequence_from_pdb(target_pdb_path, target_chain)

        # Local Boltz-2 prediction (simplified - actual API may differ)
        import boltz

        result = boltz.predict_complex(
            sequences=[binder_sequence, target_sequence],
            seed=seed,
        )

        return self._parse_boltz_result(result, binder_sequence, target_sequence, seed)

    def _predict_modal(
        self,
        binder_sequence: str,
        target_pdb_path: str,
        target_chain: str,
        seed: int,
    ) -> ComplexPredictionResult:
        """Run prediction on Modal."""
        try:
            from modal import Function

            # Get the deployed function
            predict_complex = Function.lookup("boltz2-cd3", "predict_complex")

            # Read target PDB
            with open(target_pdb_path, "r") as f:
                target_pdb_content = f.read()

            # Call Modal function
            result = predict_complex.remote(
                binder_sequence=binder_sequence,
                target_pdb_content=target_pdb_content,
                target_chain=target_chain,
                seed=seed,
            )

            return ComplexPredictionResult(
                pdb_string=result["pdb_string"],
                binder_sequence=binder_sequence,
                target_sequence=result["target_sequence"],
                pdockq=result["pdockq"],
                ptm=result["ptm"],
                plddt_mean=result["plddt_mean"],
                ipae=result.get("ipae", 0.0),
                interface_residues_binder=result.get("interface_residues_binder", []),
                interface_residues_target=result.get("interface_residues_target", []),
                num_contacts=result.get("num_contacts", 0),
                interface_area=result.get("interface_area", 0.0),
                model_version="boltz2",
                seed=seed,
            )

        except Exception as e:
            raise RuntimeError(f"Failed to run Boltz-2 on Modal: {e}")

    def _parse_boltz_result(
        self,
        result,
        binder_sequence: str,
        target_sequence: str,
        seed: int,
    ) -> ComplexPredictionResult:
        """Parse Boltz-2 output into result object."""
        # This is a placeholder - actual parsing depends on Boltz-2 API
        return ComplexPredictionResult(
            pdb_string=result.get_pdb_string(),
            binder_sequence=binder_sequence,
            target_sequence=target_sequence,
            pdockq=result.pdockq,
            ptm=result.ptm,
            plddt_mean=result.plddt_mean,
            ipae=getattr(result, "ipae", 0.0),
            interface_residues_binder=getattr(result, "interface_binder", []),
            interface_residues_target=getattr(result, "interface_target", []),
            num_contacts=getattr(result, "num_contacts", 0),
            interface_area=getattr(result, "interface_area", 0.0),
            model_version="boltz2",
            seed=seed,
        )

    def predict_batch(
        self,
        binder_sequences: list[str],
        target_pdb_path: str,
        target_chain: str = "A",
        seed: int = 42,
    ) -> list[ComplexPredictionResult]:
        """Predict complexes for multiple binders.

        Args:
            binder_sequences: List of binder sequences.
            target_pdb_path: Path to target structure.
            target_chain: Target chain ID.
            seed: Base random seed.

        Returns:
            List of ComplexPredictionResult objects.
        """
        results = []
        for i, seq in enumerate(binder_sequences):
            try:
                result = self.predict_complex(
                    binder_sequence=seq,
                    target_pdb_path=target_pdb_path,
                    target_chain=target_chain,
                    seed=seed + i,
                )
                results.append(result)
            except Exception as e:
                print(f"Warning: Failed to predict complex for sequence {i}: {e}")
                results.append(None)

        return results


def predict_binder_complex(
    binder_sequence: str,
    target_pdb_path: str,
    target_chain: str = "A",
    output_path: Optional[str] = None,
    use_modal: bool = True,
    seed: int = 42,
) -> ComplexPredictionResult:
    """Convenience function for complex prediction.

    Args:
        binder_sequence: Binder amino acid sequence.
        target_pdb_path: Path to target structure PDB.
        target_chain: Chain ID in target PDB.
        output_path: Path to save complex PDB (optional).
        use_modal: If True, use Modal for GPU compute.
        seed: Random seed.

    Returns:
        ComplexPredictionResult with predicted complex.
    """
    predictor = Boltz2Predictor(use_modal=use_modal)
    result = predictor.predict_complex(
        binder_sequence=binder_sequence,
        target_pdb_path=target_pdb_path,
        target_chain=target_chain,
        seed=seed,
    )

    if output_path:
        result.save_pdb(output_path)
        print(f"Complex structure saved to: {output_path}")

    return result


def run_calibration(
    known_binder_sequences: list[str],
    target_pdb_path: str,
    target_chain: str = "A",
    use_modal: bool = True,
    pdockq_margin: float = 0.05,
    interface_area_margin: float = 100.0,
    contacts_margin: int = 2,
) -> dict:
    """Run calibration to establish filter thresholds.

    Uses known binders to determine appropriate thresholds
    for pDockQ, interface area, and contact counts.

    Args:
        known_binder_sequences: List of known binder sequences.
        target_pdb_path: Path to target structure.
        target_chain: Target chain ID.
        use_modal: If True, use Modal.
        pdockq_margin: Margin to subtract from min pDockQ.
        interface_area_margin: Margin to subtract from min interface area.
        contacts_margin: Margin to subtract from min contacts.

    Returns:
        Dictionary with calibrated thresholds.
    """
    predictor = Boltz2Predictor(use_modal=use_modal)

    results = []
    for seq in known_binder_sequences:
        try:
            result = predictor.predict_complex(seq, target_pdb_path, target_chain)
            results.append(result)
        except Exception as e:
            print(f"Warning: Calibration failed for sequence: {e}")

    if not results:
        raise RuntimeError("Calibration failed - no successful predictions")

    # Calculate thresholds (minimum of known binders minus margin)
    pdockq_values = [r.pdockq for r in results]
    area_values = [r.interface_area for r in results]
    contact_values = [r.num_contacts for r in results]

    calibration = {
        "known_binder_results": [r.to_dict() for r in results],
        "calibrated_thresholds": {
            "min_pdockq": max(0.0, min(pdockq_values) - pdockq_margin),
            "min_interface_area": max(0.0, min(area_values) - interface_area_margin),
            "min_contacts": max(0, min(contact_values) - contacts_margin),
        },
        "margins_used": {
            "pdockq_margin": pdockq_margin,
            "interface_area_margin": interface_area_margin,
            "contacts_margin": contacts_margin,
        },
        "known_binder_stats": {
            "pdockq": {"min": min(pdockq_values), "max": max(pdockq_values), "mean": sum(pdockq_values) / len(pdockq_values)},
            "interface_area": {"min": min(area_values), "max": max(area_values), "mean": sum(area_values) / len(area_values)},
            "contacts": {"min": min(contact_values), "max": max(contact_values), "mean": sum(contact_values) / len(contact_values)},
        },
    }

    return calibration
