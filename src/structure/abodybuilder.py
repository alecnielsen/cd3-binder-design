"""ABodyBuilder2 wrapper for antibody structure prediction.

ABodyBuilder2 (ImmuneBuilder) predicts antibody Fv structures
with high accuracy, especially for CDR loops.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import tempfile


@dataclass
class ABodyBuilderResult:
    """Result from ABodyBuilder2 structure prediction."""

    pdb_string: str
    sequence_vh: str
    sequence_vl: Optional[str]
    plddt_mean: float
    plddt_per_residue: list[float]
    model_confidence: float

    def save_pdb(self, output_path: str) -> str:
        """Save structure to PDB file."""
        with open(output_path, "w") as f:
            f.write(self.pdb_string)
        return output_path


class ABodyBuilder:
    """Wrapper for ABodyBuilder2 antibody structure prediction.

    ABodyBuilder2 is part of ImmuneBuilder and uses deep learning
    to predict antibody Fv structures.
    """

    def __init__(self):
        """Initialize ABodyBuilder wrapper."""
        self._model = None
        self._available = None

    def is_available(self) -> bool:
        """Check if ABodyBuilder2 is available."""
        if self._available is not None:
            return self._available

        try:
            from ImmuneBuilder import ABodyBuilder2
            self._available = True
        except ImportError:
            self._available = False

        return self._available

    def _load_model(self):
        """Load the ABodyBuilder2 model."""
        if self._model is not None:
            return self._model

        if not self.is_available():
            raise RuntimeError(
                "ABodyBuilder2 not available. Install with: "
                "pip install ImmuneBuilder"
            )

        from ImmuneBuilder import ABodyBuilder2
        self._model = ABodyBuilder2()
        return self._model

    def predict_structure(
        self,
        vh_sequence: str,
        vl_sequence: Optional[str] = None,
        num_recycles: int = 3,
    ) -> ABodyBuilderResult:
        """Predict antibody Fv structure.

        Args:
            vh_sequence: VH (or VHH) sequence.
            vl_sequence: VL sequence (optional, None for VHH).
            num_recycles: Number of recycling iterations.

        Returns:
            ABodyBuilderResult with predicted structure.
        """
        model = self._load_model()

        # Create sequence dict for ABodyBuilder2
        if vl_sequence:
            sequences = {
                "H": vh_sequence,
                "L": vl_sequence,
            }
        else:
            # For VHH (nanobody), only provide heavy chain
            sequences = {
                "H": vh_sequence,
            }

        # Run prediction
        output = model.predict(sequences)

        # Get structure as PDB string
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
            temp_path = f.name

        output.save(temp_path)
        with open(temp_path, "r") as f:
            pdb_string = f.read()

        Path(temp_path).unlink()  # Clean up temp file

        # Extract confidence metrics
        plddt_per_residue = output.plddt.tolist() if hasattr(output, "plddt") else []
        plddt_mean = sum(plddt_per_residue) / len(plddt_per_residue) if plddt_per_residue else 0.0

        return ABodyBuilderResult(
            pdb_string=pdb_string,
            sequence_vh=vh_sequence,
            sequence_vl=vl_sequence,
            plddt_mean=plddt_mean,
            plddt_per_residue=plddt_per_residue,
            model_confidence=plddt_mean / 100.0,  # Normalize to 0-1
        )

    def predict_batch(
        self,
        sequences: list[tuple[str, Optional[str]]],
        num_recycles: int = 3,
    ) -> list[ABodyBuilderResult]:
        """Predict structures for multiple sequences.

        Args:
            sequences: List of (vh, vl) tuples. vl can be None for VHH.
            num_recycles: Number of recycling iterations.

        Returns:
            List of ABodyBuilderResult objects.
        """
        results = []
        for vh, vl in sequences:
            try:
                result = self.predict_structure(vh, vl, num_recycles)
                results.append(result)
            except Exception as e:
                print(f"Warning: Failed to predict structure: {e}")
                results.append(None)

        return results


def predict_antibody_structure(
    vh_sequence: str,
    vl_sequence: Optional[str] = None,
    output_path: Optional[str] = None,
) -> ABodyBuilderResult:
    """Convenience function for structure prediction.

    Args:
        vh_sequence: VH (or VHH) sequence.
        vl_sequence: VL sequence (optional).
        output_path: Path to save PDB file (optional).

    Returns:
        ABodyBuilderResult with predicted structure.
    """
    predictor = ABodyBuilder()
    result = predictor.predict_structure(vh_sequence, vl_sequence)

    if output_path:
        result.save_pdb(output_path)
        print(f"Structure saved to: {output_path}")

    return result
