"""HuDiff Modal deployment for post-hoc antibody humanization.

HuDiff (Tencent AI4S, Nature Machine Intelligence 2025) uses diffusion-based
framework regeneration to humanize antibody sequences while preserving CDRs.
Supports both antibodies (VH+VL) and nanobodies (VHH).

License: AFL-3.0 (Academic Free License) — permissive, compatible with project.

Deploy with:
    modal deploy modal/hudiff_app.py

Download model weights:
    modal run modal/hudiff_app.py --download

Test:
    modal run modal/hudiff_app.py --heavy-seq "EVQL..." --light-seq "DIQM..."
    modal run modal/hudiff_app.py --nanobody-seq "EVQL..."
"""

from __future__ import annotations

import json
import os
import subprocess
import tarfile
from pathlib import Path
from typing import Optional

import modal

MINUTES = 60

app = modal.App("hudiff-cd3")

# HuDiff needs PyTorch + CUDA, abnumber (ANARCI), and the HuDiff repo itself.
# PyTorch 1.13 is specified by HuDiff but newer versions work for inference.
hudiff_image = (
    modal.Image.debian_slim(python_version="3.9")
    .apt_install("git", "hmmer")
    .pip_install(
        "torch==2.1.0",
        "numpy<2",
        "pandas",
        "scipy",
        "scikit-learn",
        "tqdm",
        "easydict",
        "pyyaml",
        "einops",
        "biopython",
        "abnumber",
        "anarci",
    )
    .pip_install("sequence-models==1.6.0")
    .run_commands("git clone https://github.com/TencentAI4S/HuDiff.git /opt/hudiff")
    .env({"PYTHONPATH": "/opt/hudiff"})
)

# Persistent volume for model checkpoints (~2 GB)
hudiff_model_volume = modal.Volume.from_name("hudiff-models", create_if_missing=True)
models_dir = Path("/models/hudiff")

# Image for downloading model checkpoints
download_image = (
    modal.Image.debian_slim()
    .pip_install("huggingface-hub==0.36.0", "hf_transfer")
    .env({"HF_HUB_ENABLE_HF_TRANSFER": "1"})
)


@app.function(
    volumes={models_dir: hudiff_model_volume},
    timeout=30 * MINUTES,
    image=download_image,
)
def download_model(force_download: bool = False):
    """Download HuDiff model checkpoints from HuggingFace.

    Downloads release_data_dir.tar.gz from cloud77/HuDiff and extracts
    the checkpoint files needed for inference.
    """
    from huggingface_hub import hf_hub_download

    ckpt_ab = models_dir / "checkpoints" / "antibody" / "hudiffab.pt"
    ckpt_nb = models_dir / "checkpoints" / "nanobody" / "hudiffnb.pt"

    if ckpt_ab.exists() and ckpt_nb.exists() and not force_download:
        print("Checkpoints already exist. Use force_download=True to re-download.")
        return

    print("Downloading HuDiff checkpoints from cloud77/HuDiff...")
    tar_path = hf_hub_download(
        repo_id="cloud77/HuDiff",
        filename="release_data_dir.tar.gz",
        local_dir=str(models_dir / "download"),
    )

    print(f"Extracting checkpoints from {tar_path}...")
    with tarfile.open(tar_path, "r:gz") as tar:
        # Extract only checkpoint files (skip large LMDB training data)
        for member in tar.getmembers():
            if member.name.endswith(".pt") or member.name.endswith(".ckpt"):
                tar.extract(member, path=str(models_dir))
                print(f"  Extracted: {member.name}")

    hudiff_model_volume.commit()

    # Verify
    for ckpt, name in [(ckpt_ab, "antibody"), (ckpt_nb, "nanobody")]:
        if ckpt.exists():
            size_mb = ckpt.stat().st_size / (1024 * 1024)
            print(f"  {name} checkpoint: {size_mb:.1f} MB")
        else:
            # Try alternate paths from tar extraction
            alt_paths = list(models_dir.rglob(f"*{ckpt.name}"))
            if alt_paths:
                alt_paths[0].parent.mkdir(parents=True, exist_ok=True)
                os.rename(str(alt_paths[0]), str(ckpt))
                print(f"  {name} checkpoint moved to {ckpt}")
            else:
                print(f"  WARNING: {name} checkpoint not found!")

    print("Download complete!")


def _load_antibody_model(ckpt_path: str, device: str):
    """Load HuDiff antibody model from checkpoint."""
    import sys
    sys.path.insert(0, "/opt/hudiff")

    import torch
    from utils.train_utils import model_selected

    ckpt = torch.load(ckpt_path, map_location="cpu")
    pretrain_config = ckpt["pretrain_config"]

    model = model_selected(pretrain_config).to(device)
    model.load_state_dict(ckpt["model"])
    model.eval()

    return model, ckpt


def _load_nanobody_model(ckpt_path: str, device: str):
    """Load HuDiff nanobody model from checkpoint."""
    import sys
    sys.path.insert(0, "/opt/hudiff")

    import torch
    from utils.train_utils import model_selected
    from utils.tokenizer import Tokenizer

    ckpt = torch.load(ckpt_path, map_location="cpu")
    config = ckpt["config"]

    # Extract sub-model states
    from nanobody_scripts.nanosample import get_multi_model_state
    abnativ_state, _, infilling_state = get_multi_model_state(ckpt)

    # Load AbNatiV humanness scorer
    from model.nanoencoder.abnativ_model import AbNatiV_Model
    hparams = ckpt["abnativ_params"]
    abnativ_model = AbNatiV_Model(hparams)
    abnativ_model.load_state_dict(abnativ_state)
    abnativ_model.to(device)

    # Load infilling model
    from model.nanoencoder.model import NanoAntiTFNet
    infilling_params = ckpt["infilling_params"]
    infilling_model = NanoAntiTFNet(**infilling_params)
    infilling_model.load_state_dict(infilling_state)
    infilling_model.to(device)

    # Configure framework model
    config.model["equal_weight"] = True
    config.model["vhh_nativeness"] = False
    config.model["human_threshold"] = None
    config.model["human_all_seq"] = False
    config.model["temperature"] = False

    model_dict = {
        "abnativ": abnativ_model,
        "infilling": infilling_model,
        "target_infilling": infilling_model,
    }
    framework_model = model_selected(config, pretrained_model=model_dict, tokenizer=Tokenizer())
    model = framework_model.infilling_pretrain
    model.eval()

    return model, ckpt, config


@app.function(
    image=hudiff_image,
    volumes={models_dir: hudiff_model_volume},
    timeout=30 * MINUTES,
    gpu="A100",
)
def humanize_antibody(
    heavy_seq: str,
    light_seq: str,
    num_samples: int = 5,
    batch_size: int = 10,
    seed: int = 42,
) -> list:
    """Humanize an antibody VH+VL pair using HuDiff-Ab.

    Regenerates framework regions via diffusion while preserving CDR sequences.

    Args:
        heavy_seq: VH amino acid sequence.
        light_seq: VL amino acid sequence.
        num_samples: Number of humanized variants to generate.
        batch_size: Batch size for parallel sampling.
        seed: Random seed.

    Returns:
        List of dicts with 'vh', 'vl' keys for each humanized variant.
    """
    import sys
    sys.path.insert(0, "/opt/hudiff")

    import torch
    import numpy as np
    from abnumber import Chain
    from utils.misc import seed_all
    from antibody_scripts.sample import batch_input_element

    device = "cuda" if torch.cuda.is_available() else "cpu"
    seed_all(seed)

    ckpt_path = str(models_dir / "checkpoints" / "antibody" / "hudiffab.pt")
    model, ckpt = _load_antibody_model(ckpt_path, device)

    # IMGT-number the input sequences
    try:
        mouse_aa_h = Chain(heavy_seq, scheme="imgt").seq
        mouse_aa_l = Chain(light_seq, scheme="imgt").seq
    except Exception as e:
        print(f"IMGT numbering failed: {e}")
        print("Falling back to raw sequences")
        mouse_aa_h = heavy_seq
        mouse_aa_l = light_seq

    print(f"Input VH: {mouse_aa_h[:40]}... ({len(mouse_aa_h)} aa)")
    print(f"Input VL: {mouse_aa_l[:40]}... ({len(mouse_aa_l)} aa)")

    # Prepare input tensors
    pad_region = 0
    finetune = True
    (h_l_pad_seq_sample, h_l_pad_seq_region,
     chain_type, h_l_ms_batch, h_l_loc, ms_tokenizer) = batch_input_element(
        mouse_aa_h, mouse_aa_l, batch_size, pad_region, finetune=finetune
    )

    results = []
    duplicated_set = set()
    remaining = num_samples
    max_attempts = num_samples * 5  # Avoid infinite loops

    attempt = 0
    while remaining > 0 and attempt < max_attempts:
        attempt += 1
        # Shuffle sampling order for diversity
        np.random.shuffle(h_l_loc)

        # Re-initialize masked positions for each attempt
        (h_l_pad_seq_sample, h_l_pad_seq_region,
         chain_type, h_l_ms_batch, h_l_loc, ms_tokenizer) = batch_input_element(
            mouse_aa_h, mouse_aa_l, batch_size, pad_region, finetune=finetune
        )
        np.random.shuffle(h_l_loc)

        all_token = ms_tokenizer.toks
        with torch.no_grad():
            for i in h_l_loc:
                h_l_prediction = model(
                    h_l_pad_seq_sample.to(device),
                    h_l_pad_seq_region.to(device),
                    chain_type.to(device),
                )
                h_l_pred = h_l_prediction[:, i, :len(all_token) - 1]
                h_l_soft = torch.nn.functional.softmax(h_l_pred, dim=1)
                h_l_sample = torch.multinomial(h_l_soft, num_samples=1)
                h_l_pad_seq_sample[:, i] = h_l_sample.squeeze()

        # Split heavy (first 152 positions) and light (remaining)
        h_pad_seq_sample = h_l_pad_seq_sample[:, :152]
        l_pad_seq_sample = h_l_pad_seq_sample[:, 152:]
        h_untokenized = [ms_tokenizer.idx2seq(s) for s in h_pad_seq_sample]
        l_untokenized = [ms_tokenizer.idx2seq(s) for s in l_pad_seq_sample]

        for g_h, g_l in zip(h_untokenized, l_untokenized):
            if remaining == 0:
                break
            key = (g_h, g_l)
            if key not in duplicated_set:
                duplicated_set.add(key)
                # Remove gap characters
                vh_clean = g_h.replace("-", "")
                vl_clean = g_l.replace("-", "")
                if vh_clean and vl_clean:
                    results.append({"vh": vh_clean, "vl": vl_clean})
                    remaining -= 1
                    print(f"  Generated variant {len(results)}/{num_samples}")

    print(f"Generated {len(results)} unique antibody variants")
    return results


@app.function(
    image=hudiff_image,
    volumes={models_dir: hudiff_model_volume},
    timeout=30 * MINUTES,
    gpu="A100",
)
def humanize_nanobody(
    sequence: str,
    num_samples: int = 5,
    batch_size: int = 10,
    seed: int = 42,
) -> list:
    """Humanize a nanobody/VHH sequence using HuDiff-Nb.

    Regenerates framework regions via diffusion while preserving CDR sequences.

    Args:
        sequence: VHH amino acid sequence.
        num_samples: Number of humanized variants to generate.
        batch_size: Batch size for parallel sampling.
        seed: Random seed.

    Returns:
        List of dicts with 'sequence' key for each humanized variant.
    """
    import sys
    sys.path.insert(0, "/opt/hudiff")

    import torch
    import numpy as np
    from abnumber import Chain
    from utils.misc import seed_all
    from nanobody_scripts.nanosample import batch_input_element

    device = "cuda" if torch.cuda.is_available() else "cpu"
    seed_all(seed)

    ckpt_path = str(models_dir / "checkpoints" / "nanobody" / "hudiffnb.pt")
    model, ckpt, config = _load_nanobody_model(ckpt_path, device)

    # IMGT-number the input sequence
    try:
        vhh_seq = Chain(sequence, scheme="imgt").seq
    except Exception as e:
        print(f"IMGT numbering failed: {e}")
        print("Falling back to raw sequence")
        vhh_seq = sequence

    print(f"Input VHH: {vhh_seq[:40]}... ({len(vhh_seq)} aa)")

    results = []
    duplicated_set = set()
    remaining = num_samples
    max_attempts = num_samples * 5

    attempt = 0
    while remaining > 0 and attempt < max_attempts:
        attempt += 1

        (nano_pad_token, nano_pad_region,
         nano_loc, ms_tokenizer) = batch_input_element(
            vhh_seq, inpaint_sample=True, batch_size=batch_size
        )
        np.random.shuffle(nano_loc)

        all_token = ms_tokenizer.toks
        with torch.no_grad():
            for i in nano_loc:
                nano_prediction = model(
                    nano_pad_token.to(device),
                    nano_pad_region.to(device),
                    H_chn_type=None,
                )
                nano_pred = nano_prediction[:, i, :len(all_token) - 1]
                nano_soft = torch.nn.functional.softmax(nano_pred, dim=1)
                nano_sample = torch.multinomial(nano_soft, num_samples=1)
                nano_pad_token[:, i] = nano_sample.squeeze()

        nano_untokenized = [ms_tokenizer.idx2seq(s) for s in nano_pad_token]

        for g_h in nano_untokenized:
            if remaining == 0:
                break
            if g_h not in duplicated_set:
                # Validate via IMGT numbering
                seq_clean = g_h.replace("-", "")
                if seq_clean:
                    try:
                        Chain(seq_clean, scheme="imgt")
                        duplicated_set.add(g_h)
                        results.append({"sequence": seq_clean})
                        remaining -= 1
                        print(f"  Generated variant {len(results)}/{num_samples}")
                    except Exception:
                        pass  # Invalid sequence, skip

    print(f"Generated {len(results)} unique nanobody variants")
    return results


@app.function(
    image=hudiff_image,
    volumes={models_dir: hudiff_model_volume},
    timeout=60 * MINUTES,
    gpu="A100",
)
def batch_humanize(
    candidates: list,
) -> list:
    """Batch humanize multiple candidates.

    Each candidate dict should have:
        - design_id: str
        - binder_type: "vhh" or "fab"
        - sequence: VH or VHH sequence
        - sequence_vl: VL sequence (for fab/scfv)
        - num_samples: int (optional, default 5)

    Args:
        candidates: List of candidate dicts.

    Returns:
        List of result dicts, each with design_id and variants.
    """
    results = []
    for cand in candidates:
        design_id = cand["design_id"]
        binder_type = cand.get("binder_type", "vhh")
        num_samples = cand.get("num_samples", 5)

        print(f"\nHumanizing {design_id} ({binder_type})...")

        try:
            if binder_type in ("fab", "scfv") and cand.get("sequence_vl"):
                variants = humanize_antibody.local(
                    heavy_seq=cand["sequence"],
                    light_seq=cand["sequence_vl"],
                    num_samples=num_samples,
                    seed=hash(design_id) % 2**31,
                )
            else:
                variants = humanize_nanobody.local(
                    sequence=cand["sequence"],
                    num_samples=num_samples,
                    seed=hash(design_id) % 2**31,
                )

            results.append({
                "design_id": design_id,
                "binder_type": binder_type,
                "variants": variants,
                "num_generated": len(variants),
            })
        except Exception as e:
            print(f"  Error humanizing {design_id}: {e}")
            import traceback
            traceback.print_exc()
            results.append({
                "design_id": design_id,
                "binder_type": binder_type,
                "variants": [],
                "num_generated": 0,
                "error": str(e),
            })

    return results


@app.local_entrypoint()
def main(
    heavy_seq: str = None,
    light_seq: str = None,
    nanobody_seq: str = None,
    num_samples: int = 5,
    seed: int = 42,
    download: bool = False,
):
    """Local entrypoint for testing.

    Args:
        heavy_seq: VH sequence for antibody humanization.
        light_seq: VL sequence for antibody humanization.
        nanobody_seq: VHH sequence for nanobody humanization.
        num_samples: Number of variants to generate.
        seed: Random seed.
        download: If True, download model weights and exit.
    """
    if download:
        print("Downloading HuDiff model weights...")
        download_model.remote()
        print("Done!")
        return

    if heavy_seq and light_seq:
        print(f"Humanizing antibody (VH+VL)...")
        print(f"  VH: {heavy_seq[:40]}... ({len(heavy_seq)} aa)")
        print(f"  VL: {light_seq[:40]}... ({len(light_seq)} aa)")

        variants = humanize_antibody.remote(
            heavy_seq=heavy_seq,
            light_seq=light_seq,
            num_samples=num_samples,
            seed=seed,
        )

        print(f"\nGenerated {len(variants)} humanized variants:")
        for i, v in enumerate(variants):
            print(f"  {i+1}. VH: {v['vh'][:40]}... ({len(v['vh'])} aa)")
            print(f"      VL: {v['vl'][:40]}... ({len(v['vl'])} aa)")

    elif nanobody_seq:
        print(f"Humanizing nanobody (VHH)...")
        print(f"  VHH: {nanobody_seq[:40]}... ({len(nanobody_seq)} aa)")

        variants = humanize_nanobody.remote(
            sequence=nanobody_seq,
            num_samples=num_samples,
            seed=seed,
        )

        print(f"\nGenerated {len(variants)} humanized variants:")
        for i, v in enumerate(variants):
            print(f"  {i+1}. VHH: {v['sequence'][:40]}... ({len(v['sequence'])} aa)")

    else:
        print("Usage:")
        print("  # Download model weights (one-time)")
        print("  modal run modal/hudiff_app.py --download")
        print("")
        print("  # Humanize antibody (VH+VL)")
        print("  modal run modal/hudiff_app.py --heavy-seq 'EVQL...' --light-seq 'DIQM...'")
        print("")
        print("  # Humanize nanobody (VHH)")
        print("  modal run modal/hudiff_app.py --nanobody-seq 'EVQL...'")
        print("")
        print("  # Deploy for pipeline use")
        print("  modal deploy modal/hudiff_app.py")
