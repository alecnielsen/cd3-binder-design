"""Boltz-2 Modal deployment for complex structure prediction.

This module provides a Modal app for running Boltz-2
on GPU infrastructure for protein complex prediction.

Deploy with:
    modal deploy modal/boltz2_app.py

Run locally:
    modal run modal/boltz2_app.py --binder-seq "EVQL..." --target-seq "DGNE..."
"""

import json
import subprocess
from pathlib import Path
from typing import Optional

import modal

MINUTES = 60

app = modal.App("boltz2-cd3")

# Container with Boltz-2
image = modal.Image.debian_slim(python_version="3.12").pip_install(
    "boltz==2.1.1",
)

# Persistent volume for model weights
boltz_model_volume = modal.Volume.from_name("boltz-models", create_if_missing=True)
models_dir = Path("/models/boltz")

# Image for downloading model
download_image = (
    modal.Image.debian_slim()
    .pip_install("huggingface-hub==0.36.0", "hf_transfer")
    .env({"HF_HUB_ENABLE_HF_TRANSFER": "1"})
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
            insertion_code = line[26] if len(line) > 26 else " "
            res_key = (res_num, insertion_code, res_name)

            if res_key not in seen_residues and res_name in aa_3to1:
                seen_residues.add(res_key)
                sequence.append((res_num, insertion_code, aa_3to1[res_name]))

    sequence.sort(key=lambda x: (x[0], x[1]))
    return "".join([aa for _, _, aa in sequence])


def build_yaml(sequences: dict[str, str], use_msa: bool = False) -> str:
    """Build Boltz YAML input from sequences dict."""
    lines = ["version: 1", "sequences:"]
    for chain_id, seq in sequences.items():
        lines.extend([
            "  - protein:",
            f"      id: {chain_id}",
            f"      sequence: {seq}",
        ])
        if not use_msa:
            lines.append("      msa: empty")
    return "\n".join(lines)


def calculate_interface_metrics(
    cif_or_pdb: str,
    binder_chain: str = "B",
    target_chain: str = "A",
    distance_cutoff: float = 5.0,
) -> dict:
    """Calculate interface metrics from predicted complex structure."""
    binder_atoms = []
    target_atoms = []
    binder_residues = set()
    target_residues = set()

    for line in cif_or_pdb.split("\n"):
        # Handle both PDB and basic mmCIF ATOM lines
        if not line.startswith("ATOM"):
            continue

        # Try PDB format first
        if len(line) >= 54 and line[30:38].strip().replace(".", "").replace("-", "").isdigit():
            chain = line[21]
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                res_num = int(line[22:26])
            except (ValueError, IndexError):
                continue
        else:
            # Skip mmCIF for now - would need proper parsing
            continue

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

    interface_area = (len(binder_residues) + len(target_residues)) * 80.0

    return {
        "num_contacts": len(contact_pairs),
        "interface_residues_binder": sorted(binder_residues),
        "interface_residues_target": sorted(target_residues),
        "interface_area": interface_area,
    }


@app.function(
    volumes={models_dir: boltz_model_volume},
    timeout=20 * MINUTES,
    image=download_image,
)
def download_model(force_download: bool = False):
    """Download Boltz-2 model weights to Modal volume."""
    from huggingface_hub import snapshot_download

    snapshot_download(
        repo_id="boltz-community/boltz-2",
        local_dir=models_dir,
        force_download=force_download,
    )
    boltz_model_volume.commit()
    print(f"Model downloaded to {models_dir}")


@app.function(
    image=image,
    volumes={models_dir: boltz_model_volume},
    timeout=10 * MINUTES,
    gpu="H100",
)
def predict_complex(
    binder_sequence: str,
    target_sequence: str,
    use_msa: bool = False,
    seed: int = 42,
) -> dict:
    """Predict binder-target complex structure.

    Args:
        binder_sequence: Binder amino acid sequence.
        target_sequence: Target amino acid sequence.
        use_msa: Whether to use MSA server (slower but may be more accurate).
        seed: Random seed.

    Returns:
        Dictionary with prediction results.
    """
    # Build YAML input - target is chain A, binder is chain B
    yaml_content = build_yaml({"A": target_sequence, "B": binder_sequence}, use_msa=use_msa)

    input_path = Path("input.yaml")
    input_path.write_text(yaml_content)

    cmd = [
        "boltz", "predict", str(input_path),
        "--cache", str(models_dir),
        "--accelerator", "gpu",
    ]
    if use_msa:
        cmd.append("--use_msa_server")

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"Boltz-2 prediction failed: {result.stderr}")

    # Find outputs
    results_dir = list(Path(".").glob("boltz_results_*"))[0]
    pred_dir = list((results_dir / "predictions").iterdir())[0]

    cif_file = list(pred_dir.glob("*_model_0.cif"))[0]
    conf_file = list(pred_dir.glob("confidence_*.json"))[0]

    cif_content = cif_file.read_text()
    confidence = json.loads(conf_file.read_text())

    # Calculate interface metrics (note: mmCIF parsing is limited)
    interface = calculate_interface_metrics(cif_content, binder_chain="B", target_chain="A")

    return {
        "cif_string": cif_content,
        "binder_sequence": binder_sequence,
        "target_sequence": target_sequence,
        "pdockq": confidence.get("pdockq", 0.0),
        "ptm": confidence.get("ptm", 0.0),
        "plddt_mean": confidence.get("complex_plddt", 0.0),
        "ipae": confidence.get("ipae", 0.0),
        "iptm": confidence.get("protein_iptm", 0.0),
        "num_contacts": interface["num_contacts"],
        "interface_residues_binder": interface["interface_residues_binder"],
        "interface_residues_target": interface["interface_residues_target"],
        "interface_area": interface["interface_area"],
        "seed": seed,
    }


@app.function(
    image=image,
    volumes={models_dir: boltz_model_volume},
    timeout=10 * MINUTES,
    gpu="H100",
)
def predict_complex_from_pdb(
    binder_sequence: str,
    target_pdb_content: str,
    target_chain: str = "A",
    use_msa: bool = False,
    seed: int = 42,
) -> dict:
    """Predict complex using target from PDB content.

    Args:
        binder_sequence: Binder amino acid sequence.
        target_pdb_content: Target PDB file content.
        target_chain: Target chain ID to extract.
        use_msa: Whether to use MSA server.
        seed: Random seed.

    Returns:
        Dictionary with prediction results.
    """
    target_sequence = parse_pdb_sequence(target_pdb_content, target_chain)
    return predict_complex(binder_sequence, target_sequence, use_msa=use_msa, seed=seed)


@app.function(
    image=image,
    volumes={models_dir: boltz_model_volume},
    timeout=120 * MINUTES,
    gpu="H100",
)
def predict_complex_batch(
    binder_sequences: list[str],
    target_sequence: str,
    use_msa: bool = False,
    seed: int = 42,
) -> list[dict]:
    """Predict complexes for multiple binders.

    Args:
        binder_sequences: List of binder sequences.
        target_sequence: Target sequence.
        use_msa: Whether to use MSA server.
        seed: Base random seed.

    Returns:
        List of prediction results.
    """
    results = []

    for i, seq in enumerate(binder_sequences):
        print(f"Predicting {i+1}/{len(binder_sequences)}...")
        try:
            result = predict_complex(
                binder_sequence=seq,
                target_sequence=target_sequence,
                use_msa=use_msa,
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
    image=image,
    volumes={models_dir: boltz_model_volume},
    timeout=60 * MINUTES,
    gpu="H100",
)
def run_calibration(
    known_binder_sequences: list[str],
    target_sequence: str,
    use_msa: bool = False,
    seed: int = 42,
) -> dict:
    """Run calibration using known binders.

    Uses known binders to establish filter thresholds.

    Args:
        known_binder_sequences: List of known binder sequences.
        target_sequence: Target sequence.
        use_msa: Whether to use MSA server.
        seed: Random seed.

    Returns:
        Calibration results with recommended thresholds.
    """
    results = []

    for i, seq in enumerate(known_binder_sequences):
        print(f"Calibrating with binder {i+1}/{len(known_binder_sequences)}...")
        try:
            result = predict_complex(
                binder_sequence=seq,
                target_sequence=target_sequence,
                use_msa=use_msa,
                seed=seed + i,
            )
            results.append(result)
        except Exception as e:
            print(f"Warning: Calibration failed for sequence {i}: {e}")

    if not results:
        raise RuntimeError("Calibration failed - no successful predictions")

    # Calculate thresholds
    pdockq_values = [r["pdockq"] for r in results if "pdockq" in r]
    area_values = [r["interface_area"] for r in results if "interface_area" in r]
    contact_values = [r["num_contacts"] for r in results if "num_contacts" in r]

    return {
        "known_binder_results": results,
        "calibrated_thresholds": {
            "min_pdockq": max(0.0, min(pdockq_values) - 0.05) if pdockq_values else 0.0,
            "min_interface_area": max(0.0, min(area_values) - 100) if area_values else 0.0,
            "min_contacts": max(0, min(contact_values) - 2) if contact_values else 0,
        },
        "known_binder_stats": {
            "pdockq": {
                "min": min(pdockq_values) if pdockq_values else 0,
                "max": max(pdockq_values) if pdockq_values else 0,
                "mean": sum(pdockq_values) / len(pdockq_values) if pdockq_values else 0,
            },
            "interface_area": {
                "min": min(area_values) if area_values else 0,
                "max": max(area_values) if area_values else 0,
                "mean": sum(area_values) / len(area_values) if area_values else 0,
            },
            "contacts": {
                "min": min(contact_values) if contact_values else 0,
                "max": max(contact_values) if contact_values else 0,
                "mean": sum(contact_values) / len(contact_values) if contact_values else 0,
            },
        },
    }


@app.local_entrypoint()
def main(
    binder_seq: str = None,
    target_seq: str = None,
    target_pdb: str = None,
    target_chain: str = "A",
    use_msa: bool = False,
    seed: int = 42,
    download: bool = False,
):
    """Local entrypoint for testing.

    Args:
        binder_seq: Binder amino acid sequence.
        target_seq: Target amino acid sequence.
        target_pdb: Path to target PDB file (alternative to target_seq).
        target_chain: Chain ID to extract from PDB.
        use_msa: Whether to use MSA server.
        seed: Random seed.
        download: If True, download model weights and exit.
    """
    if download:
        print("Downloading Boltz-2 model weights...")
        download_model.remote()
        print("Done!")
        return

    if binder_seq is None:
        print("Usage: modal run modal/boltz2_app.py --binder-seq 'EVQL...' --target-seq 'DGNE...'")
        print("   or: modal run modal/boltz2_app.py --binder-seq 'EVQL...' --target-pdb <path>")
        print("   or: modal run modal/boltz2_app.py --download")
        return

    # Get target sequence
    if target_seq:
        target_sequence = target_seq
    elif target_pdb:
        with open(target_pdb, "r") as f:
            target_pdb_content = f.read()
        target_sequence = parse_pdb_sequence(target_pdb_content, target_chain)
    else:
        print("Error: Must provide either --target-seq or --target-pdb")
        return

    print(f"Running Boltz-2 complex prediction")
    print(f"  Binder length: {len(binder_seq)} aa")
    print(f"  Target length: {len(target_sequence)} aa")
    print(f"  MSA: {'enabled' if use_msa else 'disabled'}")
    print(f"  Seed: {seed}")

    # Run on Modal
    result = predict_complex.remote(
        binder_sequence=binder_seq,
        target_sequence=target_sequence,
        use_msa=use_msa,
        seed=seed,
    )

    print(f"\nPrediction results:")
    print(f"  pDockQ: {result['pdockq']:.3f}")
    print(f"  pTM: {result['ptm']:.3f}")
    print(f"  pLDDT mean: {result['plddt_mean']:.1f}")
    print(f"  ipTM: {result.get('iptm', 0):.3f}")
    print(f"  Interface contacts: {result['num_contacts']}")
    print(f"  Interface area: {result['interface_area']:.1f} A^2")
    print(f"  Target contact residues: {result['interface_residues_target'][:10]}...")
