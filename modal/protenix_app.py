"""Protenix structure prediction on Modal.

Protenix (Apache 2.0, https://github.com/bytedance/protenix) is ByteDance's
open-source AlphaFold3-style structure prediction. Used as cross-validation
for top candidates after Boltz-2 filtering.

Deploy with:
    modal deploy modal/protenix_app.py

Warm up (pre-download weights):
    modal run modal/protenix_app.py --warmup
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import modal

MINUTES = 60

app = modal.App("protenix-cd3")

# Protenix auto-downloads model weights on first prediction.
# We cache them in a persistent volume so subsequent runs are fast.
protenix_cache_volume = modal.Volume.from_name("protenix-cache", create_if_missing=True)
cache_dir = Path("/cache/protenix")

protenix_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.1-devel-ubuntu22.04",
        add_python="3.11",
    )
    .apt_install("kalign", "hmmer")
    .pip_install("protenix")
    .env({
        "CUDA_HOME": "/usr/local/cuda",
    })
)


@app.function(
    image=protenix_image,
    gpu="H100",
    volumes={str(cache_dir): protenix_cache_volume},
    timeout=30 * MINUTES,
)
def predict_complex(
    binder_sequence: str,
    target_sequence: str,
    binder_type: str = "vhh",
    binder_sequence_vl: str | None = None,
    seed: int = 101,
    use_msa: bool = False,
    model_name: str = "protenix_base_default_v1.0.0",
) -> dict:
    """Predict complex structure using Protenix.

    Args:
        binder_sequence: VHH or VH sequence.
        target_sequence: CD3 target sequence.
        binder_type: "vhh" or "scfv"/"fab".
        binder_sequence_vl: VL sequence (for Fab/scFv).
        seed: Random seed.
        use_msa: Whether to use MSA search (slower, not needed for cross-validation).
        model_name: Protenix model name.

    Returns:
        Dict with iptm, ptm, plddt_mean, ranking_score, chain_pair_iptm, cif_string.
    """
    # Point Protenix cache to persistent volume
    os.environ["PROTENIX_ROOT_DIR"] = str(cache_dir)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Build Protenix input JSON
        chains = []

        if binder_type == "vhh":
            chains.append({"proteinChain": {"sequence": binder_sequence, "count": 1}})
        else:
            # Fab/scFv: VH and VL as separate chains
            chains.append({"proteinChain": {"sequence": binder_sequence, "count": 1}})
            if binder_sequence_vl:
                chains.append({"proteinChain": {"sequence": binder_sequence_vl, "count": 1}})

        # Target chain
        chains.append({"proteinChain": {"sequence": target_sequence, "count": 1}})

        input_data = [
            {
                "name": "complex",
                "modelSeeds": [seed],
                "sequences": chains,
            }
        ]

        input_path = tmpdir / "input.json"
        with open(input_path, "w") as f:
            json.dump(input_data, f, indent=2)

        output_dir = tmpdir / "output"
        output_dir.mkdir()

        # Run Protenix prediction via CLI
        cmd = [
            "protenix", "pred",
            "-i", str(input_path),
            "-o", str(output_dir),
            "-s", str(seed),
            "-n", model_name,
            "--use_default_params", "true",
        ]

        if use_msa:
            cmd.extend(["--use_msa", "true"])

        print(f"Running Protenix: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=25 * 60)

        # Commit volume to persist any newly-downloaded weights
        protenix_cache_volume.commit()

        if result.returncode != 0:
            error_msg = result.stderr[-1000:] if result.stderr else "unknown error"
            print(f"Protenix stderr: {error_msg}")
            return {
                "error": f"Protenix failed: {error_msg[-500:]}",
                "iptm": None,
                "ptm": None,
                "plddt_mean": None,
                "ranking_score": None,
                "chain_pair_iptm": None,
                "cif_string": None,
            }

        # Parse output
        return _parse_protenix_output(output_dir, seed)


def _parse_protenix_output(output_dir: Path, seed: int) -> dict:
    """Parse Protenix output files for metrics and structure."""
    result = {
        "iptm": None,
        "ptm": None,
        "plddt_mean": None,
        "ranking_score": None,
        "chain_pair_iptm": None,
        "cif_string": None,
        "error": None,
    }

    # Find output files - Protenix outputs to output_dir/complex/seed_N/
    complex_dir = output_dir / "complex"
    if not complex_dir.exists():
        # Try flat structure
        complex_dir = output_dir

    # Look for confidence JSON
    confidence_files = list(complex_dir.rglob("*confidence*.json")) + list(complex_dir.rglob("*summary*.json"))
    for conf_file in confidence_files:
        try:
            with open(conf_file) as f:
                conf = json.load(f)

            # Protenix outputs these fields
            if isinstance(conf, dict):
                result["iptm"] = conf.get("iptm", conf.get("interface_ptm"))
                result["ptm"] = conf.get("ptm")
                result["plddt_mean"] = conf.get("plddt_mean", conf.get("mean_plddt"))
                result["ranking_score"] = conf.get("ranking_score")
                result["chain_pair_iptm"] = conf.get("chain_pair_iptm")
            elif isinstance(conf, list) and conf:
                c = conf[0]
                result["iptm"] = c.get("iptm", c.get("interface_ptm"))
                result["ptm"] = c.get("ptm")
                result["plddt_mean"] = c.get("plddt_mean", c.get("mean_plddt"))
                result["ranking_score"] = c.get("ranking_score")
                result["chain_pair_iptm"] = c.get("chain_pair_iptm")
            break
        except (json.JSONDecodeError, KeyError):
            continue

    # Find CIF output
    cif_files = list(complex_dir.rglob("*.cif"))
    if cif_files:
        try:
            result["cif_string"] = cif_files[0].read_text()
        except Exception:
            pass

    return result


@app.function(
    image=protenix_image,
    gpu="H100",
    volumes={str(cache_dir): protenix_cache_volume},
    timeout=60 * MINUTES,
)
def batch_predict(
    candidates: list[dict],
    target_sequence: str,
    seed: int = 101,
    use_msa: bool = False,
    model_name: str = "protenix_base_default_v1.0.0",
) -> list[dict]:
    """Batch predict structures for multiple candidates.

    Each candidate dict should have:
        - design_id: str
        - sequence: str (VH or VHH)
        - sequence_vl: Optional[str]
        - binder_type: str

    Args:
        candidates: List of candidate dicts.
        target_sequence: CD3 target sequence.
        seed: Random seed.
        use_msa: Whether to use MSA.
        model_name: Protenix model name.

    Returns:
        List of result dicts (one per candidate), each with design_id and metrics.
    """
    results = []
    for cand in candidates:
        print(f"Predicting {cand['design_id']}...")
        try:
            pred = predict_complex.local(
                binder_sequence=cand["sequence"],
                target_sequence=target_sequence,
                binder_type=cand.get("binder_type", "vhh"),
                binder_sequence_vl=cand.get("sequence_vl"),
                seed=seed,
                use_msa=use_msa,
                model_name=model_name,
            )
            pred["design_id"] = cand["design_id"]
            results.append(pred)
        except Exception as e:
            results.append({
                "design_id": cand["design_id"],
                "error": str(e),
                "iptm": None,
                "ptm": None,
                "plddt_mean": None,
                "ranking_score": None,
                "chain_pair_iptm": None,
                "cif_string": None,
            })
    return results


@app.function(
    image=protenix_image,
    gpu="H100",
    volumes={str(cache_dir): protenix_cache_volume},
    timeout=30 * MINUTES,
)
def warmup():
    """Pre-download Protenix model weights by running a minimal prediction.

    This ensures weights are cached in the persistent volume for fast subsequent runs.
    """
    os.environ["PROTENIX_ROOT_DIR"] = str(cache_dir)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Minimal input: two short protein chains
        input_data = [
            {
                "name": "warmup",
                "modelSeeds": [42],
                "sequences": [
                    {"proteinChain": {"sequence": "EVQLVESGGGLVQPGGSLRL", "count": 1}},
                    {"proteinChain": {"sequence": "QDGNEEMGGITQTP", "count": 1}},
                ],
            }
        ]

        input_path = tmpdir / "input.json"
        with open(input_path, "w") as f:
            json.dump(input_data, f, indent=2)

        output_dir = tmpdir / "output"
        output_dir.mkdir()

        cmd = [
            "protenix", "pred",
            "-i", str(input_path),
            "-o", str(output_dir),
            "-s", "42",
            "-n", "protenix_base_default_v1.0.0",
            "--use_default_params", "true",
            "-e", "1",  # Single sample for speed
        ]

        print("Running warmup prediction to cache model weights...")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=25 * 60)

        protenix_cache_volume.commit()

        if result.returncode == 0:
            print("Warmup complete! Model weights cached.")
        else:
            print(f"Warmup prediction failed (weights may still be cached): {result.stderr[-500:]}")


@app.local_entrypoint()
def main(warmup_flag: bool = False):
    if warmup_flag:
        print("Running warmup to pre-download Protenix weights...")
        warmup.remote()
        print("Done!")
    else:
        print("Usage:")
        print("  modal run modal/protenix_app.py --warmup-flag  # Pre-download weights")
        print("  modal deploy modal/protenix_app.py              # Deploy app")
        print()
        print("Call predict_complex() or batch_predict() from pipeline scripts.")
