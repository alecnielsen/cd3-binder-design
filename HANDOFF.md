# Handoff Notes - Modal Integration

## Status: Modal Apps Rewritten

Both `modal/boltz2_app.py` and `modal/boltzgen_app.py` have been rewritten to use the actual CLI interfaces (via subprocess), matching the working pattern from `protein-predict`.

### What Changed

**boltz2_app.py:**
- Now uses `boltz predict input.yaml` via subprocess
- Builds YAML input with sequences
- Parses output CIF and confidence JSON
- Added model download function with persistent Modal volume
- Uses H100 GPU (same as protein-predict)

**boltzgen_app.py:**
- Now uses `boltzgen run spec.yaml --protocol nanobody-anything` via subprocess
- Builds design spec YAML with target sequence and designable binder positions
- Parses output FASTA/CIF files for designed sequences
- Added model download function with persistent Modal volume
- Uses A100 GPU

### Key Implementation Details

Both apps follow this pattern:
```python
# 1. Build YAML input
yaml_content = build_yaml(sequences)
Path("input.yaml").write_text(yaml_content)

# 2. Call CLI via subprocess
cmd = ["boltz", "predict", "input.yaml", "--accelerator", "gpu"]
result = subprocess.run(cmd, capture_output=True, text=True)

# 3. Parse output files
cif_file = list(Path(".").glob("boltz_results_*/predictions/*/*.cif"))[0]
```

### Next Steps

1. **Download model weights** (first time only):
   ```bash
   modal run modal/boltz2_app.py --download
   modal run modal/boltzgen_app.py --download
   ```

2. **Test the apps**:
   ```bash
   # Test Boltz-2
   modal run modal/boltz2_app.py --binder-seq "EVQLVES..." --target-seq "DGNE..."

   # Test BoltzGen
   modal run modal/boltzgen_app.py --target-pdb data/targets/cd3.pdb --target-chain A
   ```

3. **Run calibration** once Modal apps are verified working:
   ```bash
   python scripts/00_run_calibration.py
   ```

4. **Run full pipeline**:
   ```bash
   python scripts/run_full_pipeline.py --config config.yaml
   ```

### Potential Issues to Watch

1. **BoltzGen output format** - The output parsing assumes FASTA files; actual output format may differ. Check `boltzgen run --help` for exact output structure.

2. **Model weight paths** - The `--cache` flag paths may need adjustment based on actual BoltzGen/Boltz CLI expectations.

3. **HuggingFace repo IDs** - The download functions assume `boltz-community/boltz-2` and `boltz-community/boltzgen`. Verify these are correct.

4. **mmCIF parsing** - The interface metrics calculation uses simple PDB parsing; mmCIF output from Boltz-2 may need proper parsing for accurate contact calculations.

### References

- Working Boltz implementation: `/Users/alec/kernel/protein-predict/boltz_modal.py`
- BoltzGen GitHub: https://github.com/HannesStark/boltzgen
- Boltz GitHub: https://github.com/jwohlwend/boltz
