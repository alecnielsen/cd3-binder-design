# Handoff Notes - Modal Integration

## Current Status (2026-01-25)

### Modal Apps: WORKING

Both Modal apps work when called directly:

```bash
# BoltzGen - generates antibody-like sequences
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --target-chain A --num-designs 5

# Boltz-2 - predicts complex structures
modal run modal/boltz2_app.py --binder-seq "QVQLVE..." --target-seq "QTPYK..."
```

### Pipeline Integration: PARTIALLY WORKING

Calibration runs but the full pipeline has issues (see "Remaining Issues" below).

```bash
# This works
PYTHONPATH=/Users/alec/kernel/cd3-binder-design python3 scripts/00_run_calibration.py

# This has issues
PYTHONPATH=/Users/alec/kernel/cd3-binder-design python3 scripts/run_full_pipeline.py --config config.yaml
```

## Session Progress (2026-01-25)

### Fixed Issues

1. **BoltzGen YAML format** - Complete rewrite to use correct BoltzGen specification
   - Use `sequence: 110..130` (range notation) instead of `sequence: XXX...`
   - Reference target via `file: { path: target.pdb }` instead of inline sequence
   - Remove invalid keys (`design: true`, `binding: [A]`)

2. **Modal API updates** - Updated from deprecated `Function.lookup()` to `Function.from_name()`

3. **Module shadowing** - Removed `modal/__init__.py` that was shadowing the `modal` pip package

4. **Config YAML corruption** - Changed `yaml.dump()` to `yaml.safe_dump()` to avoid `!!python/tuple` tags

5. **Boltz-2 inter-function calls** - Added `.local()` for calling between Modal functions

6. **Modal function parameter mismatches** - Updated `boltzgen_runner.py` to pass correct parameters

### Remaining Issues

1. **Interface metrics return 0** - `calculate_interface_metrics()` in boltz2_app.py only parses PDB format, but Boltz-2 outputs mmCIF. pTM/ipTM/pLDDT work fine.

2. **Some scripts don't accept --config** - `01_setup_targets.py` doesn't take `--config` but `run_full_pipeline.py` passes it. Use `--skip-setup` as workaround.

3. **Parameter mismatches in src/ code** - The `src/design/boltzgen_runner.py` may still have mismatches with Modal function signatures. Partially fixed but needs verification.

4. **Calibration shows identical results** - All 3 calibration binders show same pTM/pLDDT (possible caching issue or Boltz-2 not seeing different sequences).

5. **Long execution time** - 200 designs takes very long. Reduced to 5 in config for testing.

## Files Changed This Session

```
modal/boltzgen_app.py      - Fixed YAML format, CIF parsing, CSV metrics
modal/boltz2_app.py        - Added .local() calls for inter-function calls
modal/__init__.py          - DELETED (was shadowing modal package)
src/structure/boltz_complex.py - Updated Modal API, use predict_complex directly
src/design/boltzgen_runner.py  - Fixed Modal check, parameter mapping
src/pipeline/config.py     - Use yaml.safe_dump()
config.yaml                - Fixed target paths, reduced design count for testing
data/outputs/calibration.json - Calibration results (interface metrics are 0)
```

## Quick Test Commands

```bash
# Set Python path
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

# Test BoltzGen directly (WORKS)
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --target-chain A --num-designs 5

# Test Boltz-2 directly (WORKS)
modal run modal/boltz2_app.py --binder-seq "QVQLVESGGGLVQAGGSLRLSEAASGFTFSSYGMGWVRQAPGKGREWVAA" --target-seq "QTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCENCME"

# Run calibration (WORKS but metrics are 0)
python3 scripts/00_run_calibration.py

# Run de novo design step (NEEDS TESTING)
python3 scripts/02_run_denovo_design.py --config config.yaml
```

## Next Steps to Complete Pipeline

1. **Verify de novo design step** - Run with small count (5 designs) and verify output
2. **Fix interface metrics** - Update `calculate_interface_metrics()` to parse mmCIF
3. **Debug calibration** - Investigate why all binders show same scores
4. **Test remaining pipeline steps** - Steps 03-07 haven't been tested
5. **Increase design count** - Once working, increase to 100-200 designs

## Technical Notes

### Modal Apps Deployed
- `boltzgen-cd3` - De novo binder design
- `boltz2-cd3` - Complex structure prediction

### Model Weights (in Modal volumes)
- `boltzgen-models` - BoltzGen weights (~10GB)
- `boltz-models` - Boltz-2 weights (~2GB)

### Key Configuration
- `config.yaml` - Main pipeline config
- `data/targets/1XIW.pdb` - CD3εδ structure
- `data/targets/1SY6.pdb` - CD3εγ + OKT3 structure
