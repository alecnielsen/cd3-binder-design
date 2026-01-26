# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-01-26)

### Pipeline: FULLY WORKING

All pipeline steps are tested and working:

```bash
# Set Python path
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

# Full pipeline execution
python3 scripts/00_run_calibration.py              # ✅ Working - sets filter thresholds
python3 scripts/02_run_denovo_design.py --config config.yaml  # ✅ Working - BoltzGen designs
python3 scripts/03_run_optimization.py --config config.yaml   # ✅ Working - humanization/variants
python3 scripts/04_predict_structures.py --config config.yaml # ✅ Working - Boltz-2 predictions
python3 scripts/05_filter_candidates.py --config config.yaml  # ✅ Working - filtering cascade
python3 scripts/06_format_bispecifics.py --config config.yaml # ✅ Working - bispecific formats
python3 scripts/07_generate_report.py --config config.yaml    # ✅ Working - HTML/JSON reports
```

### Modal Apps: DEPLOYED AND WORKING

```bash
# BoltzGen - de novo binder design
python3 -m modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --target-chain A --num-designs 5

# Boltz-2 - complex structure prediction
python3 -m modal run modal/boltz2_app.py --binder-seq "QVQLVE..." --target-seq "QTPYK..."
```

## Session Progress (2026-01-26)

### Fixed Issues

1. **mmCIF parsing** - Added `parse_mmcif_atoms()` function to properly parse Boltz-2 mmCIF output
   - Interface metrics (contacts, area) now calculate correctly

2. **Identical calibration results** - Fixed two issues:
   - Added cleanup of old `boltz_results_*` directories before each prediction
   - Added `--seed` parameter to boltz command
   - Results now show different pTM/pLDDT for different binders

3. **Python 3.9 compatibility** - Fixed type hint syntax in `src/utils/constants.py`
   - Changed `tuple[int, int] | None` to `Optional[Tuple[int, int]]`

### Calibration Results (3 known binders)

| Binder | pTM | pLDDT | Contacts | Interface Area |
|--------|-----|-------|----------|----------------|
| Teplizumab | 0.944 | 0.964 | 42 | 2560 Å² |
| SP34 | 0.753 | 0.917 | 30 | 2160 Å² |
| UCHT1 | 0.784 | 0.941 | 32 | 2240 Å² |

Calibrated thresholds: min_contacts=28, min_interface_area=2060 Å²

## Known Limitations

1. **pDockQ always 0** - Boltz-2 v2.1.1 doesn't compute pDockQ by default. Use ipTM instead.

2. **ANARCI/BioPhi not installed** - Humanness scoring skipped (soft-fail). Install for full functionality:
   ```bash
   pip install anarci biophi
   ```

3. **BoltzGen metrics not captured** - ipTM/pTM from BoltzGen not being passed through to final designs (confidence=0). The sequences are valid.

## Output Files

```
data/outputs/
├── calibration.json          # Calibration results with thresholds
├── denovo/                   # BoltzGen designs
├── optimized/                # Humanized/affinity variants
├── structures/               # Boltz-2 predictions
├── filtered/                 # Filtered candidates
├── formatted/                # Bispecific constructs
│   ├── crossmab_constructs.json
│   ├── fab_scfv_constructs.json
│   └── ...
└── reports/                  # Final reports
    ├── report_*.html
    ├── report_*.json
    └── scorecards/
```

## Configuration

Key settings in `config.yaml`:
- `design.num_vhh_designs: 5` - Increase for production (100-200)
- `design.num_scfv_designs: 5` - Increase for production (100-200)
- `filtering.use_calibrated: true` - Uses calibrated thresholds
- `output.num_final_candidates: 10` - Final candidate count

## Next Steps

1. **Increase design count** - Set `num_vhh_designs` and `num_scfv_designs` to 100-200 for production
2. **Install optional packages** - `pip install anarci biophi` for humanness scoring
3. **Capture BoltzGen metrics** - Update boltzgen_runner.py to extract ipTM/pTM from output
4. **Review candidates** - Check generated reports in `data/outputs/reports/`

## Future Enhancements: Affinity Prediction

Currently, the pipeline uses **structural confidence metrics** (pTM, ipTM, contacts) which do NOT predict binding affinity. The following open-source tools could be integrated:

### Recommended Tools (All Permissively Licensed)

| Tool | License | Input | Predicts | Priority |
|------|---------|-------|----------|----------|
| **Boltz-2 IC50** | MIT | Structure | log10(IC50) | High - already available |
| **AttABseq** | MIT | Sequences | ΔΔG (mutation effects) | Medium |
| **Graphinity** | BSD-3 | Structures | ΔΔG (mutation effects) | Medium |
| **ANTIPASTI** | CC BY 4.0 | Structure | Binding affinity | Low |

### Implementation Notes

1. **Boltz-2 IC50** (lowest effort):
   - Add `--sampling_steps_affinity 200` to boltz command
   - Add `properties: { affinity: { ligand: B } }` to YAML
   - Parse `affinity_output.json` from results
   - NOT validated for antibodies - use with caution

2. **AttABseq** (sequence-only ΔΔG):
   - GitHub: https://github.com/ruofanjin/AttABseq
   - Useful for comparing affinity variants (WT vs 10x weaker)
   - No structure needed
   - `pip install` not available - requires cloning repo

3. **Graphinity** (structure-based ΔΔG):
   - GitHub: https://github.com/oxpig/Graphinity
   - From Oxford Protein Informatics Group
   - Requires both WT and mutant structures
   - Pearson r=0.87 on benchmarks

### Critical Caveat

**None of these reliably predict absolute Kd.** They predict:
- Structural confidence (current pipeline)
- Relative affinity change (ΔΔG) from mutations
- IC50 estimates (assay-dependent, not intrinsic Kd)

**Experimental validation with SPR/BLI remains essential.**
