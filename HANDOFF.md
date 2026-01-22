# Handoff Notes - BoltzGen Integration

## Current Issue: `modal/boltzgen_app.py` Has Incorrect API

The existing `modal/boltzgen_app.py` assumes a BoltzGen API that doesn't exist:

```python
# WRONG - This API doesn't exist:
from boltzgen import BoltzGen
model = BoltzGen.load()
designs = model.design(target_structure=..., binder_type="scfv", ...)
```

### Actual BoltzGen API

BoltzGen uses a **CLI-based workflow** with YAML design specifications:

```bash
boltzgen run design_spec.yaml \
    --protocol protein-anything \
    --output output_dir \
    --num_designs 100
```

Available protocols:
- `protein-anything` - Generic protein binder
- `nanobody-anything` - VHH nanobodies (~120 aa)
- `antibody-anything` - Antibody CDR design

### Key Finding: No scFv Linker in BoltzGen

**BoltzGen does NOT automatically add (G4S)3 linkers to outputs.**

Evidence:
- Searched entire BoltzGen package for "linker" and "GGGGS" - nothing found
- Design spec YAML has no linker configuration option
- BoltzGen produces single-domain binders, not VH-linker-VL scFvs

The (G4S)3 linker in `config.yaml:73` is for **downstream bispecific assembly**, not BoltzGen.

## Required Changes

### Option A: Fix `modal/boltzgen_app.py` to use actual BoltzGen CLI

```python
# Use subprocess to call boltzgen CLI
import subprocess
import yaml

def run_boltzgen(target_pdb_content, num_designs, ...):
    # 1. Write design spec YAML
    design_spec = {
        "entities": [
            {"protein": {"sequence": target_sequence, "design": False}},
            {"protein": {"sequence": "X" * 120, "design": True, "binding": ["A"]}}
        ]
    }

    # 2. Call boltzgen CLI
    subprocess.run(["boltzgen", "run", spec_path, "--protocol", "nanobody-anything", ...])

    # 3. Parse output CIF/FASTA files
```

### Option B: Remove `binder_type="scfv"` Support

Since BoltzGen doesn't natively produce scFvs:
- Only support `binder_type="vhh"` (nanobody)
- Add scFv linker in post-processing if needed

## Design Spec YAML Format

Valid keys: `entities`, `protein`, `id`, `sequence`, `design`, `binding`, `constraints`, `total_len`, `min`, `max`

Example:
```yaml
entities:
  - protein:
      id: A
      sequence: DGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVC
  - protein:
      id: B
      sequence: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      design: true
      binding: [A]
```

## Next Steps

1. Decide: Fix boltzgen_app.py to use CLI, or simplify to VHH-only?
2. If scFv needed: Add post-processing to concatenate VH + linker + VL
3. Update tests to match actual BoltzGen behavior

## References

- [BoltzGen GitHub](https://github.com/HannesStark/boltzgen)
- [BoltzGen CLI help](run `boltzgen run --help` in Modal container)
- `config.yaml:72-74` - downstream linker definitions
