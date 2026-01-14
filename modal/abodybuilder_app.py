"""ABodyBuilder2 Modal deployment for antibody structure prediction.

This module provides a Modal app for running ABodyBuilder2
on GPU infrastructure for antibody Fv structure prediction.

Deploy with:
    modal deploy modal/abodybuilder_app.py

Run locally:
    modal run modal/abodybuilder_app.py::predict_structure --vh-seq "EVQL..." --vl-seq "DIQM..."
"""

import modal

# Create Modal app
app = modal.App("abodybuilder2-cd3")

# Define container image with ABodyBuilder2 dependencies
abodybuilder_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "torch>=2.0.0",
        "numpy",
        "biopython",
        "einops",
        "scipy",
    )
    .pip_install("ImmuneBuilder")
)


@app.function(
    image=abodybuilder_image,
    gpu="T4",  # T4 is sufficient for antibody structure prediction
    timeout=600,  # 10 min timeout
    retries=2,
)
def predict_structure(
    vh_sequence: str,
    vl_sequence: str = None,
    num_recycles: int = 3,
) -> dict:
    """Predict antibody Fv structure.

    Args:
        vh_sequence: VH (or VHH) sequence.
        vl_sequence: VL sequence (optional, None for VHH).
        num_recycles: Number of recycling iterations.

    Returns:
        Dictionary with prediction results.
    """
    import tempfile
    import os
    from pathlib import Path

    try:
        from ImmuneBuilder import ABodyBuilder2

        # Load model
        model = ABodyBuilder2()

        # Prepare sequences
        if vl_sequence:
            sequences = {"H": vh_sequence, "L": vl_sequence}
        else:
            sequences = {"H": vh_sequence}

        # Run prediction
        output = model.predict(sequences)

        # Save to temp file and read back
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
            temp_path = f.name

        output.save(temp_path)
        with open(temp_path, "r") as f:
            pdb_string = f.read()

        Path(temp_path).unlink()

        # Extract confidence metrics
        plddt_per_residue = output.plddt.tolist() if hasattr(output, "plddt") else []
        plddt_mean = sum(plddt_per_residue) / len(plddt_per_residue) if plddt_per_residue else 0.0

        return {
            "pdb_string": pdb_string,
            "vh_sequence": vh_sequence,
            "vl_sequence": vl_sequence,
            "plddt_mean": plddt_mean,
            "plddt_per_residue": plddt_per_residue,
            "model_confidence": plddt_mean / 100.0,
            "is_vhh": vl_sequence is None,
        }

    except Exception as e:
        raise RuntimeError(f"ABodyBuilder2 prediction failed: {e}")


@app.function(
    image=abodybuilder_image,
    gpu="T4",
    timeout=3600,  # 1 hour for batch
)
def predict_structure_batch(
    sequences: list[dict],  # [{"vh": "...", "vl": "..."}, ...]
    num_recycles: int = 3,
) -> list[dict]:
    """Predict structures for multiple antibodies.

    Args:
        sequences: List of dicts with "vh" and optional "vl" keys.
        num_recycles: Number of recycling iterations.

    Returns:
        List of prediction results.
    """
    results = []

    for seq_dict in sequences:
        try:
            result = predict_structure.local(
                vh_sequence=seq_dict["vh"],
                vl_sequence=seq_dict.get("vl"),
                num_recycles=num_recycles,
            )
            results.append(result)
        except Exception as e:
            results.append({
                "error": str(e),
                "vh_sequence": seq_dict["vh"],
                "vl_sequence": seq_dict.get("vl"),
            })

    return results


@app.local_entrypoint()
def main(
    vh_seq: str = None,
    vl_seq: str = None,
    output_pdb: str = None,
):
    """Local entrypoint for testing.

    Args:
        vh_seq: VH (or VHH) sequence.
        vl_seq: VL sequence (optional).
        output_pdb: Path to save output PDB (optional).
    """
    if vh_seq is None:
        print("Usage: modal run modal/abodybuilder_app.py --vh-seq 'EVQL...' [--vl-seq 'DIQM...']")
        return

    print(f"Running ABodyBuilder2 structure prediction")
    print(f"  VH length: {len(vh_seq)} aa")
    if vl_seq:
        print(f"  VL length: {len(vl_seq)} aa")
    else:
        print("  Type: VHH (no VL)")

    # Run on Modal
    result = predict_structure.remote(
        vh_sequence=vh_seq,
        vl_sequence=vl_seq,
    )

    print(f"\nPrediction results:")
    print(f"  pLDDT mean: {result['plddt_mean']:.1f}")
    print(f"  Model confidence: {result['model_confidence']:.3f}")
    print(f"  Is VHH: {result['is_vhh']}")

    if output_pdb:
        with open(output_pdb, "w") as f:
            f.write(result["pdb_string"])
        print(f"  Structure saved to: {output_pdb}")
