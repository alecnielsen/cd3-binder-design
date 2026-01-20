# Modal GPU Deployments

This directory contains Modal apps for GPU-intensive operations in the CD3 binder design pipeline.

## Overview

Modal (https://modal.com) provides on-demand GPU compute for:
- **BoltzGen**: De novo binder design
- **Boltz-2**: Protein complex structure prediction
- **ABodyBuilder2**: Antibody structure prediction (optional)

## Why Modal?

1. **CUDA Required**: BoltzGen and Boltz-2 require NVIDIA GPUs with CUDA support
2. **No Local Setup**: Avoids complex local GPU configuration
3. **Pay-per-use**: ~$0.001-0.01 per design, ~$0.05 per complex prediction
4. **Scalable**: Can parallelize batch processing across multiple GPUs

**Note**: Apple MPS (Metal Performance Shaders) is NOT compatible with BoltzGen.

## Apps

### boltzgen_app.py

De novo binder design using BoltzGen.

```bash
# Deploy
modal deploy modal/boltzgen_app.py

# Test locally
modal run modal/boltzgen_app.py --target-pdb data/targets/cd3.pdb --num-designs 10
```

**Functions:**
- `run_boltzgen()`: Generate designs for a single target
- `run_boltzgen_batch()`: Generate designs for multiple targets

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `target_pdb_content` | str | required | PDB file content as string |
| `target_chain` | str | "A" | Target chain ID |
| `binder_type` | str | "vhh" | "vhh" or "scfv" |
| `num_designs` | int | 100 | Number of designs to generate |
| `hotspot_residues` | list[int] | None | Optional residues to target |
| `seed` | int | 42 | Random seed |
| `temperature` | float | 1.0 | Sampling temperature |
| `num_recycles` | int | 3 | Structure prediction recycles |

**GPU Configuration:**
- GPU: A100 (40GB)
- Timeout: 1 hour (single), 2 hours (batch)
- Retries: 2

### boltz2_app.py

Protein complex structure prediction using Boltz-2.

```bash
# Deploy
modal deploy modal/boltz2_app.py

# Test locally
modal run modal/boltz2_app.py --binder-seq "EVQL..." --target-pdb data/targets/cd3.pdb
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

**GPU Configuration:**
- GPU: A100 (40GB)
- Timeout: 30 min (single), 2 hours (batch), 1 hour (calibration)
- Retries: 2

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
| BoltzGen (100 VHH) | ~20 min | $3-5 |
| BoltzGen (100 scFv) | ~30 min | $5-8 |
| Boltz-2 (1 complex) | ~2 min | $0.05 |
| Boltz-2 (100 complexes) | ~3 hours | $15-20 |
| Calibration (3 binders) | ~10 min | $1-2 |

**Full pipeline estimate**: $20-50 for 400 designs + structure predictions.

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
- ABodyBuilder2/ImmuneBuilder: BSD 3-Clause

Modal itself requires an account but has no license restrictions on the tools run within it.
