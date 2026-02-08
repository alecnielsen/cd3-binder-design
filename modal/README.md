# Modal GPU Deployments

This directory contains Modal apps for GPU-intensive operations in the CD3 binder design pipeline.

## Overview

Modal (https://modal.com) provides on-demand GPU compute for:
- **BoltzGen**: De novo binder design
- **Boltz-2**: Protein complex structure prediction
- **Protenix**: Cross-validation structure prediction (step 05b)
- **ABodyBuilder2**: Antibody structure prediction (optional)

## Why Modal?

1. **CUDA Required**: BoltzGen and Boltz-2 require NVIDIA GPUs with CUDA support
2. **No Local Setup**: Avoids complex local GPU configuration
3. **Pay-per-use**: ~$0.001-0.01 per design, ~$0.05 per complex prediction
4. **Scalable**: Can parallelize batch processing across multiple GPUs

**Note**: Apple MPS (Metal Performance Shaders) is NOT compatible with BoltzGen.

## Apps

### boltzgen_app.py

De novo binder design using BoltzGen. Supports two design modes:

| Mode | Protocol | Output | Use Case |
|------|----------|--------|----------|
| **VHH** | `nanobody-anything` | Single-domain ~120 aa | De novo nanobody |
| **Fab** | `antibody-anything` | VH ~120 aa + VL ~107 aa | CDR redesign on human scaffolds |

```bash
# Download model weights (first time only)
modal run modal/boltzgen_app.py --download

# Deploy
modal deploy modal/boltzgen_app.py

# VHH design (default)
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --design-type vhh --num-designs 10

# Fab CDR redesign (requires scaffold files)
python scripts/setup_fab_scaffolds.py  # Run once to download scaffolds
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --design-type fab --scaffolds adalimumab,belimumab --num-designs 10
```

#### Critical: Target Chain Extraction

**BoltzGen automatically extracts only the specified target chain before design.**

This is critical because crystal structures often contain bound antibodies:
- **1XIW**: Contains CD3εδ + UCHT1 Fab (8 chains total)
- **1SY6**: Contains CD3γ/ε fusion + OKT3 Fab (3 chains total)

If the full multi-chain PDB were passed to BoltzGen, designs could:
1. Target the wrong protein (e.g., UCHT1 instead of CD3)
2. Target an epitope only accessible when no antibody is bound
3. Generate antibody-like sequences by "learning" from the bound Fab

The `extract_chain_from_pdb_content()` function ensures only the target chain
(e.g., chain A = CD3ε) is passed to the model.

#### BoltzGen YAML Format

BoltzGen uses a specific YAML format for design specifications:

```yaml
entities:
  - protein:
      id: B
      sequence: 110..130  # Range notation for binder length
  - file:
      path: target.pdb    # File reference for target
      include:
        - chain:
            id: A

# Optional: specify binding site residues
binding_types:
  - chain:
      id: A
      binding: 23,24,25,50,51
```

**Important**: Do NOT use `sequence: XXX...` placeholders or `design: true` keys - these are not valid BoltzGen format.

**Functions:**
- `run_boltzgen()`: VHH design for a single target (extracts chain first)
- `run_boltzgen_fab()`: Fab CDR redesign using human antibody scaffolds
- `run_boltzgen_batch()`: Generate designs for multiple targets
- `extract_chain_from_pdb_content()`: Extract single chain from multi-chain PDB
- `build_design_spec_yaml()`: Generate correct YAML format for VHH
- `build_fab_target_yaml()`: Generate YAML format for Fab CDR redesign

**Parameters for `run_boltzgen()` (VHH):**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `target_pdb_content` | str | required | PDB file content (can be multi-chain) |
| `target_chain` | str | "A" | Chain ID to extract for design |
| `num_designs` | int | 100 | Number of designs to generate |
| `binder_length` | int | 120 | Binder length in aa (110-130 range used) |
| `hotspot_residues` | list[int] | None | Optional target residues for binding |
| `protocol` | str | "nanobody-anything" | BoltzGen protocol |
| `seed` | int | 42 | Random seed |

**Parameters for `run_boltzgen_fab()` (Fab CDR redesign):**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `target_pdb_content` | str | required | PDB file content (can be multi-chain) |
| `target_chain` | str | "A" | Chain ID to extract for design |
| `scaffold_names` | list[str] | required | Scaffold names (e.g., ["adalimumab"]) |
| `scaffold_yaml_contents` | dict | required | {scaffold_name: yaml_content} |
| `scaffold_cif_contents` | dict | required | {scaffold_name: cif_content} |
| `num_designs` | int | 100 | Number of designs to generate |
| `seed` | int | 42 | Random seed |

**Fab CDR Redesign Output:**
| Field | Description |
|-------|-------------|
| `vh_sequence` | VH chain sequence (~120 aa) |
| `vl_sequence` | VL chain sequence (~107 aa) |
| `pTM` | Predicted TM-score for Fab structure |
| `ipTM` | Interface predicted TM-score (binding quality) |

**Output Metrics:**
| Metric | Description |
|--------|-------------|
| `ipTM` | Interface predicted TM-score (0-1, higher = better binding) |
| `pTM` | Predicted TM-score for binder structure |
| `pae_min` | Minimum predicted aligned error at interface |
| `rmsd` | RMSD between designed and refolded structure |

**GPU Configuration:**
- GPU: A100 (40GB)
- Timeout: 1 hour (single), 2 hours (batch)

### boltz2_app.py

Protein complex structure prediction using Boltz-2.

```bash
# Download model weights (first time only)
modal run modal/boltz2_app.py --download

# Deploy
modal deploy modal/boltz2_app.py

# Run prediction (sequences only - no PDB needed)
modal run modal/boltz2_app.py --binder-seq "EVQLVESGGGLVQ..." --target-seq "QTPYKVSISGT..."
```

**Functions:**
- `predict_complex()`: Predict single binder-target complex
- `predict_complex_batch()`: Predict complexes for multiple binders
- `run_calibration()`: Calibrate thresholds using known binders

**Output Metrics:**
| Metric | Description |
|--------|-------------|
| `pdockq` | Structural confidence score (0-1, higher = more confident) |
| `ptm` | Predicted TM-score |
| `plddt_mean` | Mean pLDDT across complex |
| `ipae` | Interface predicted aligned error |
| `num_contacts` | Unique residue pairs within 5A |
| `interface_area` | Estimated buried surface area |
| `interface_residues_target` | CD3 residues at interface |

**Critical**: pDockQ is a structural confidence metric, NOT an affinity predictor. High pDockQ does not mean high binding affinity.

**PDB Parsing Notes:**
- Handles standard amino acids (20 canonical)
- Properly handles insertion codes (e.g., residue 52A comes after 52)
- Filters by chain ID
- Returns sequences in residue number order

**GPU Configuration:**
- GPU: A100 (40GB)
- Timeout: 30 min (single), 2 hours (batch), 1 hour (calibration)
- Retries: 2

### protenix_app.py

Cross-validation structure prediction using Protenix (Apache 2.0, ByteDance). Used in step 05b to re-predict top candidates after filtering, providing an orthogonal check against Boltz-2 predictions.

```bash
# Deploy (builds CUDA devel image — takes ~5 min first time)
modal deploy modal/protenix_app.py

# Pre-cache model weights via warmup prediction (first time only)
modal run modal/protenix_app.py --warmup-flag
```

**Important**: Protenix requires `nvidia/cuda:12.4.1-devel-ubuntu22.04` base image (not `debian_slim`). Protenix JIT-compiles custom CUDA kernels (e.g., `fast_layer_norm_cuda_v2`) which require the `nvcc` compiler. The `CUDA_HOME` env var is set to `/usr/local/cuda`.

**Functions:**
- `predict_complex()`: Predict single binder-target complex structure
- `batch_predict()`: Predict complexes for multiple candidates sequentially (60 min timeout)
- `warmup()`: Run minimal prediction to trigger auto-download and cache model weights

**Note**: There is no separate `download_weights` command. Protenix auto-downloads weights on first prediction. The `warmup()` function handles this by running a minimal 2-chain prediction.

**Input:** Binder sequence(s) + target sequence (pure sequence, no PDB needed)
- VHH: 2 chains (binder + target)
- Fab/scFv: 3 chains (VH + VL + target)

**Output Metrics:**
| Metric | Description |
|--------|-------------|
| `iptm` | Interface predicted TM-score |
| `ptm` | Predicted TM-score |
| `plddt_mean` | Mean pLDDT |
| `ranking_score` | Protenix ranking score |
| `chain_pair_iptm` | Per-chain-pair ipTM matrix |
| `cif_string` | Predicted structure in CIF format |

**GPU Configuration:**
- GPU: H100
- Image: `nvidia/cuda:12.4.1-devel-ubuntu22.04` with Python 3.11
- Dependencies: `protenix` (pip), `kalign` + `hmmer` (apt)
- Timeout: 30 min per prediction, 60 min for batch
- MSA: Disabled by default (`use_msa=False`) for speed
- Volume: `protenix-cache` for persistent model weight storage

**Cross-validation:** Candidates where Boltz-2 and Protenix ipTM disagree by >0.1 are flagged in the validation report.

### abodybuilder_app.py

Antibody structure prediction using ABodyBuilder2/ImmuneBuilder.

```bash
# Deploy
modal deploy modal/abodybuilder_app.py
```

**Note**: ABodyBuilder2 can also run locally without GPU. Modal deployment is optional for this tool.

## Deployment

### Initial Setup

```bash
# Install Modal CLI
pip install modal

# Authenticate
modal setup

# Deploy all apps
modal deploy modal/boltzgen_app.py
modal deploy modal/boltz2_app.py
modal deploy modal/protenix_app.py
modal deploy modal/abodybuilder_app.py
```

### Container Images

Each app defines its own container image with required dependencies:

```python
boltzgen_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "torch>=2.0.0",
        "numpy",
        "biopython",
        "einops",
        "scipy",
        "boltzgen",
    )
)
```

**Important**: Pin package versions in production for reproducibility.

## Usage from Pipeline

The pipeline scripts call Modal functions via the src wrappers:

```python
# From src/structure/boltz_complex.py
from modal import lookup

def predict_complex(binder_seq: str, target_pdb: str, use_modal: bool = True):
    if use_modal:
        predict_fn = lookup("boltz2-cd3", "predict_complex")
        return predict_fn.remote(
            binder_sequence=binder_seq,
            target_pdb_content=open(target_pdb).read(),
        )
    else:
        # Local fallback (requires GPU)
        ...
```

## Cost Estimation

| Operation | GPU Time | Estimated Cost |
|-----------|----------|----------------|
| BoltzGen VHH (100 designs) | ~20 min | $3-5 |
| BoltzGen Fab CDR redesign (100 designs) | ~30 min | $5-8 |
| Boltz-2 (1 complex) | ~2 min | $0.05 |
| Boltz-2 (100 complexes) | ~3 hours | $15-20 |
| Protenix (10 candidates) | ~30-60 min | $5-10 |
| Calibration (3 binders) | ~10 min | $1-2 |

**Full pipeline estimate**: $25-60 for 400 designs (200 VHH + 200 Fab) + structure predictions + validation.

## Reproducibility

For reproducible results:

1. **Random seeds**: All functions accept a `seed` parameter
2. **Version pinning**: Pin boltzgen, boltz, torch versions
3. **Container hashes**: Record image hashes in config

```yaml
# config.yaml
modal:
  boltzgen_image_hash: "sha256:abc123..."
  boltz2_image_hash: "sha256:def456..."
```

## Troubleshooting

### GPU Not Available

```
Error: No GPU available
```

Ensure your Modal account has GPU quota. Free tier has limits.

### Timeout Errors

```
Error: Function timed out
```

- Increase timeout in `@app.function()` decorator
- Reduce batch size
- Check for memory issues (OOM)

### Import Errors

```
ModuleNotFoundError: No module named 'boltzgen'
```

The package may not be available on PyPI yet. Check:
- BoltzGen installation instructions
- May require `pip install git+https://github.com/...`

## License

All tools used have permissive licenses:
- BoltzGen: MIT
- Boltz-2: MIT
- Protenix: Apache 2.0
- ABodyBuilder2/ImmuneBuilder: BSD 3-Clause

Modal itself requires an account but has no license restrictions on the tools run within it.
