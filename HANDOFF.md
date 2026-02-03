# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-02-03)

### Pipeline: VALIDATION TEST PASSED (VHH Track)

VHH de novo design validated. Fab track has known issues (scaffold chain IDs).

```bash
conda activate cd3-binder
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

python3 scripts/setup_fab_scaffolds.py             # ✅ Done - scaffolds downloaded
python3 scripts/00_run_calibration.py              # ✅ Working
python3 scripts/02_run_denovo_design.py --config config.yaml  # ✅ VHH working, Fab blocked
python3 scripts/03_run_optimization.py --config config.yaml   # ✅ Working (reformats known antibodies)
python3 scripts/04_predict_structures.py --config config.yaml # ✅ Working
python3 scripts/05_filter_candidates.py --config config.yaml  # ✅ Working
python3 scripts/06_format_bispecifics.py --config config.yaml # ✅ Working
python3 scripts/07_generate_report.py --config config.yaml    # ✅ Working
```

---

## Design Types

| Type | Method | Output | Status |
|------|--------|--------|--------|
| **VHH** | BoltzGen `nanobody-anything` | Single-domain ~120 aa | ✅ Validated |
| **Fab** | BoltzGen `antibody-anything` CDR redesign | VH ~120 aa + VL ~107 aa | ⚠️ Blocked (scaffold chain IDs) |
| **scFv from known Ab** | Optimization track | VH-linker-VL | ✅ Working |

**IMPORTANT**: The optimization track does NOT generate new designs. It reformats known antibody sequences (teplizumab, SP34, UCHT1) from `data/starting_sequences/*.yaml` as scFv for structure prediction. The only tracks producing **novel** binders are VHH and Fab de novo design.

---

## Validation Test Results (2026-02-03)

### Pipeline Run Summary

| Step | Status | Output |
|------|--------|--------|
| De novo design | ✅ | 10 VHH designs generated |
| Structure prediction | ✅ | 22 candidates predicted (10 new VHH + 2 old + 10 optimized) |
| Filtering | ✅ | 10 candidates passed (via fallback relaxation) |
| Bispecific formatting | ✅ | 29 constructs across 5 formats |
| Report generation | ✅ | 12 report files |

### Top Candidates

| Rank | Candidate | Score | Type |
|------|-----------|-------|------|
| 1 | sp34_original | 0.382 | Known Ab (scFv) |
| 2 | sp34_original_WT | 0.382 | Known Ab (scFv) |
| 3 | sp34_humanized | 0.382 | Known Ab (scFv) |
| 4 | vhh_1XIW_0002 | 0.354 | **De novo VHH** |
| 5 | vhh_1SY6_0003 | 0.353 | **De novo VHH** |

### Key Observations

1. **VHH design reliably produces 10 designs** - Confirmed across 4 consecutive runs
2. **pDockQ shows 0.000** - This is expected; pDockQ is NOT a native Boltz-2 metric (it's AlphaFold-Multimer specific). Use pTM/ipTM instead.
3. **Fallback filtering works** - 0 candidates passed strict thresholds, 10 passed after relaxation
4. **De novo VHH designs score competitively** - Within range of known antibody scFvs

---

## Known Issues

### 1. Fab Scaffold Chain ID Mismatch - BLOCKING

**Root cause**: The `setup_fab_scaffolds.py` defines chain IDs (H/L) that don't match the actual chain labels in downloaded CIF files from RCSB PDB.

**Example**: dupilumab YAML expects chains H/L, but 6WGB.cif has chains A/B/C/D.

**Error**: `ValueError: Specified chain id H not in file /root/fab_scaffolds/dupilumab.6wgb.cif`

**Affected scaffolds**: Most scaffolds except adalimumab (which uses B/A correctly).

**Fix needed**: Update chain mappings in `scripts/setup_fab_scaffolds.py` by checking actual chain IDs in each PDB structure.

### 2. Modal Connection Timeouts - INTERMITTENT

**Symptom**: `grpclib.exceptions.StreamTerminatedError: Connection lost` during long Fab jobs.

**Cause**: Fab CDR redesign has 90-minute timeout; Modal client connection drops before completion.

**Workaround**: Run Fab design separately with shorter batches or increased client timeout.

### 3. BoltzGen Modal Path Bug - FIXED ✅

**Root cause**: YAML path was being doubled (`fab_scaffolds/fab_scaffolds/...`).

**Fix**: Updated `modal/boltzgen_app.py` line 625 to keep relative path unchanged.

---

## Fixes Applied

### Session 2026-02-03

1. **Modal package installation** - Added to conda environment (was missing)
2. **BoltzGen Modal deployment** - Deployed app + downloaded model weights
3. **Path doubling bug** - Fixed in `modal/boltzgen_app.py`
4. **VHH-only validation** - Set `num_fab_designs: 0` to bypass Fab issues

### Session 2026-02-02

1. **Fab Scaffold Files** - Ran `setup_fab_scaffolds.py` to download CIF + YAML files
2. **Empty Sequence Bug** - Added handler for `vhh_designs`/`fab_designs` keys
3. **Conda Environment** - Created `cd3-binder` with ANARCI + Sapiens

---

## Next Steps for Production

### Before Production Run

1. **Fix Fab scaffold chain IDs** - Check each PDB and update `setup_fab_scaffolds.py`
2. **Re-enable Fab design** - Set `num_fab_designs: 10+` in config
3. **Test Fab end-to-end** - Verify at least adalimumab scaffold works

### Production Config

```yaml
design:
  num_vhh_designs: 50
  num_fab_designs: 50  # After scaffold fix
  fab_scaffolds:
    - adalimumab  # Only this one has correct chain IDs currently
```

### Experimental Validation

After computational pipeline:
1. SPR/BLI binding assays for top candidates
2. T-cell activation (CD69/CD25)
3. Cytokine release panel (IL-2, IFN-γ, TNF-α)
4. Expression titer assessment

---

## Output Files

```
data/outputs/
├── calibration.json          # Calibration thresholds
├── denovo/
│   ├── denovo_results_20260126_090617.json  # Old trial (2 VHH)
│   └── denovo_results_20260203_063803.json  # Validation (10 VHH)
├── optimized/                # Known antibody scFvs
├── structures/               # Boltz-2 predictions (22 candidates)
├── filtered/                 # 10 filtered candidates
├── formatted/                # 29 bispecific constructs
└── reports/                  # Scorecards and HTML report
```

---

## Key Documentation

- `docs/reference/calibration-methodology.md` - Explains scFv constructs and metrics vs affinity
- `CLAUDE.md` - Implementation notes for Claude sessions
- `README.md` - Full project documentation

---

## Calibration Results (Reference)

From known CD3 binder scFvs:

| Binder | pTM | pLDDT | Contacts | Interface Area |
|--------|-----|-------|----------|----------------|
| Teplizumab-scFv | 0.944 | 96.4 | 42 | 2560 Å² |
| SP34-scFv | 0.753 | 91.7 | 30 | 2160 Å² |
| UCHT1-scFv | 0.784 | 94.1 | 32 | 2240 Å² |

**Important**: These are VH-linker-VL scFvs, NOT full IgG. pTM/pLDDT are structural confidence scores, NOT binding affinity predictors.
