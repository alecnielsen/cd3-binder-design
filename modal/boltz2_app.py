"""Boltz-2 Modal deployment for complex structure prediction.

This module provides a Modal app for running Boltz-2
on GPU infrastructure for protein complex prediction.

Deploy with:
    modal deploy modal/boltz2_app.py

Run locally:
    modal run modal/boltz2_app.py::predict_complex --binder-seq "EVQL..." --target-pdb data/targets/cd3.pdb
"""

import modal

# Create Modal app
app = modal.App("boltz2-cd3")

# Define container image with Boltz-2 dependencies
boltz2_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "torch>=2.0.0",
        "numpy",
        "biopython",
        "einops",
        "scipy",
    )
    # Note: Boltz installation may require additional steps
    .pip_install("boltz")
)


def parse_pdb_sequence(pdb_content: str, chain_id: str = "A") -> str:
    """Extract amino acid sequence from PDB content."""
    aa_3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    sequence = []
    seen_residues = set()

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            chain = line[21]
            if chain != chain_id:
                continue

            res_name = line[17:20].strip()
            res_num = int(line[22:26])
            res_key = (res_num, res_name)

            if res_key not in seen_residues and res_name in aa_3to1:
                seen_residues.add(res_key)
                sequence.append((res_num, aa_3to1[res_name]))

    sequence.sort(key=lambda x: x[0])
    return "".join([aa for _, aa in sequence])


def calculate_interface_metrics(
    pdb_string: str,
    binder_chain: str = "B",
    target_chain: str = "A",
    distance_cutoff: float = 5.0,
) -> dict:
    """Calculate interface metrics from predicted complex."""
    # Parse atom coordinates
    binder_atoms = []
    target_atoms = []
    binder_residues = set()
    target_residues = set()

    for line in pdb_string.split("\n"):
        if not line.startswith("ATOM"):
            continue

        chain = line[21]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        res_num = int(line[22:26])

        if chain == binder_chain:
            binder_atoms.append((res_num, x, y, z))
        elif chain == target_chain:
            target_atoms.append((res_num, x, y, z))

    # Find contacts - count unique residue pairs
    contact_pairs = set()
    cutoff_sq = distance_cutoff ** 2

    for res_b, xb, yb, zb in binder_atoms:
        for res_t, xt, yt, zt in target_atoms:
            dist_sq = (xb - xt) ** 2 + (yb - yt) ** 2 + (zb - zt) ** 2
            if dist_sq <= cutoff_sq:
                contact_pairs.add((res_b, res_t))
                binder_residues.add(res_b)
                target_residues.add(res_t)

    # Estimate interface area (rough approximation)
    interface_area = (len(binder_residues) + len(target_residues)) * 50.0

    return {
        "num_contacts": len(contact_pairs),
        "interface_residues_binder": sorted(binder_residues),
        "interface_residues_target": sorted(target_residues),
        "interface_area": interface_area,
    }


@app.function(
    image=boltz2_image,
    gpu="A100",
    timeout=1800,  # 30 min timeout
    retries=2,
)
def predict_complex(
    binder_sequence: str,
    target_pdb_content: str,
    target_chain: str = "A",
    seed: int = 42,
) -> dict:
    """Predict binder-target complex structure.

    Args:
        binder_sequence: Binder amino acid sequence.
        target_pdb_content: Target PDB file content.
        target_chain: Target chain ID.
        seed: Random seed.

    Returns:
        Dictionary with prediction results.
    """
    import tempfile
    import os

    # Extract target sequence
    target_sequence = parse_pdb_sequence(target_pdb_content, target_chain)

    try:
        # Import Boltz
        import boltz

        # Predict complex structure
        result = boltz.predict_complex(
            sequences=[binder_sequence, target_sequence],
            seed=seed,
        )

        # Get PDB string
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            temp_path = f.name

        result.save(temp_path)
        with open(temp_path, "r") as f:
            pdb_string = f.read()

        os.unlink(temp_path)

        # Calculate interface metrics
        interface = calculate_interface_metrics(pdb_string)

        return {
            "pdb_string": pdb_string,
            "binder_sequence": binder_sequence,
            "target_sequence": target_sequence,
            "pdockq": getattr(result, "pdockq", 0.0),
            "ptm": getattr(result, "ptm", 0.0),
            "plddt_mean": getattr(result, "plddt_mean", 0.0),
            "ipae": getattr(result, "ipae", 0.0),
            "num_contacts": interface["num_contacts"],
            "interface_residues_binder": interface["interface_residues_binder"],
            "interface_residues_target": interface["interface_residues_target"],
            "interface_area": interface["interface_area"],
            "seed": seed,
        }

    except Exception as e:
        raise RuntimeError(f"Boltz-2 prediction failed: {e}")


@app.function(
    image=boltz2_image,
    gpu="A100",
    timeout=7200,  # 2 hours for batch
)
def predict_complex_batch(
    binder_sequences: list[str],
    target_pdb_content: str,
    target_chain: str = "A",
    seed: int = 42,
) -> list[dict]:
    """Predict complexes for multiple binders.

    Args:
        binder_sequences: List of binder sequences.
        target_pdb_content: Target PDB content.
        target_chain: Target chain ID.
        seed: Base random seed.

    Returns:
        List of prediction results.
    """
    results = []

    for i, seq in enumerate(binder_sequences):
        try:
            result = predict_complex.local(
                binder_sequence=seq,
                target_pdb_content=target_pdb_content,
                target_chain=target_chain,
                seed=seed + i,
            )
            results.append(result)
        except Exception as e:
            results.append({
                "error": str(e),
                "binder_sequence": seq,
                "seed": seed + i,
            })

    return results


@app.function(
    image=boltz2_image,
    gpu="A100",
    timeout=3600,
)
def run_calibration(
    known_binder_sequences: list[str],
    target_pdb_content: str,
    target_chain: str = "A",
    seed: int = 42,
) -> dict:
    """Run calibration using known binders.

    Uses known binders to establish filter thresholds.

    Args:
        known_binder_sequences: List of known binder sequences.
        target_pdb_content: Target PDB content.
        target_chain: Target chain ID.
        seed: Random seed.

    Returns:
        Calibration results with recommended thresholds.
    """
    results = []

    for i, seq in enumerate(known_binder_sequences):
        try:
            result = predict_complex.local(
                binder_sequence=seq,
                target_pdb_content=target_pdb_content,
                target_chain=target_chain,
                seed=seed + i,
            )
            results.append(result)
        except Exception as e:
            print(f"Warning: Calibration failed for sequence {i}: {e}")

    if not results:
        raise RuntimeError("Calibration failed - no successful predictions")

    # Calculate thresholds
    pdockq_values = [r["pdockq"] for r in results]
    area_values = [r["interface_area"] for r in results]
    contact_values = [r["num_contacts"] for r in results]

    return {
        "known_binder_results": results,
        "calibrated_thresholds": {
            "min_pdockq": min(pdockq_values) - 0.05,
            "min_interface_area": min(area_values) - 100,
            "min_contacts": max(0, min(contact_values) - 2),
        },
        "known_binder_stats": {
            "pdockq": {
                "min": min(pdockq_values),
                "max": max(pdockq_values),
                "mean": sum(pdockq_values) / len(pdockq_values),
            },
            "interface_area": {
                "min": min(area_values),
                "max": max(area_values),
                "mean": sum(area_values) / len(area_values),
            },
            "contacts": {
                "min": min(contact_values),
                "max": max(contact_values),
                "mean": sum(contact_values) / len(contact_values),
            },
        },
    }


@app.local_entrypoint()
def main(
    binder_seq: str = None,
    target_pdb: str = None,
    seed: int = 42,
):
    """Local entrypoint for testing.

    Args:
        binder_seq: Binder amino acid sequence.
        target_pdb: Path to target PDB file.
        seed: Random seed.
    """
    if binder_seq is None or target_pdb is None:
        print("Usage: modal run modal/boltz2_app.py --binder-seq 'EVQL...' --target-pdb <path>")
        return

    # Read target PDB
    with open(target_pdb, "r") as f:
        target_content = f.read()

    print(f"Running Boltz-2 complex prediction")
    print(f"  Binder length: {len(binder_seq)} aa")
    print(f"  Target: {target_pdb}")
    print(f"  Seed: {seed}")

    # Run on Modal
    result = predict_complex.remote(
        binder_sequence=binder_seq,
        target_pdb_content=target_content,
        seed=seed,
    )

    print(f"\nPrediction results:")
    print(f"  pDockQ: {result['pdockq']:.3f}")
    print(f"  pTM: {result['ptm']:.3f}")
    print(f"  pLDDT mean: {result['plddt_mean']:.1f}")
    print(f"  Interface contacts: {result['num_contacts']}")
    print(f"  Interface area: {result['interface_area']:.1f} A^2")
    print(f"  Target contact residues: {result['interface_residues_target'][:10]}...")
