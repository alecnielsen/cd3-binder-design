# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-01-26)

### Pipeline: WORKING (with caveats)

All pipeline steps run successfully, but there's a **critical issue with de novo scFv design** (see below).

```bash
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

python3 scripts/00_run_calibration.py              # ✅ Working
python3 scripts/02_run_denovo_design.py --config config.yaml  # ⚠️ VHH only (scFv broken)
python3 scripts/03_run_optimization.py --config config.yaml   # ✅ Working
python3 scripts/04_predict_structures.py --config config.yaml # ✅ Working
python3 scripts/05_filter_candidates.py --config config.yaml  # ✅ Working
python3 scripts/06_format_bispecifics.py --config config.yaml # ✅ Working
python3 scripts/07_generate_report.py --config config.yaml    # ✅ Working
```

---

## CRITICAL: De Novo scFv Design is Broken

### The Problem

The current code sets `binder_length=250` with `nanobody-anything` protocol for scFv:

```python
# src/design/boltzgen_runner.py lines 236-241
if self.config.binder_type == "vhh":
    protocol = "nanobody-anything"
    binder_length = 120
else:  # scfv
    protocol = "nanobody-anything"  # WRONG - same protocol
    binder_length = 250  # Just makes longer single-domain blob
```

**This does NOT produce a proper scFv (VH-linker-VL).** It produces a ~250 aa single-domain protein.

### The Solution: BoltzGen Fab CDR Redesign

BoltzGen v0.2.0 supports **Fab (VH+VL) CDR redesign** using proven human scaffolds:

```yaml
# Example: adalimumab Fab scaffold
path: adalimumab.6cr1.cif
include:
  - chain:
      id: B  # VH (residues 1-121)
      res_index: 1..121
  - chain:
      id: A  # VL (residues 1-107)
      res_index: 1..107

design:  # Only CDRs are redesigned
  - chain:
      id: B
      res_index: 26..32,52..57,99..110  # H-CDR1, H-CDR2, H-CDR3
  - chain:
      id: A
      res_index: 24..34,50..56,89..97   # L-CDR1, L-CDR2, L-CDR3
```

**14 proven human Fab scaffolds available:**
adalimumab, belimumab, dupilumab, golimumab, guselkumab, nirsevimab, sarilumab, secukinumab, tezepelumab, tralokinumab, ustekinumab, mab1, necitumumab, crenezumab

### Recommended Fix

1. **VHH de novo**: Keep using `nanobody-anything` (works correctly)
2. **Fab/scFv design**: Implement BoltzGen Fab CDR redesign with human scaffolds
3. **Remove broken scFv option**: Don't pretend `binder_length=250` makes scFv

---

## Session Progress (2026-01-26)

### Fixed Issues

1. **mmCIF parsing** - Added `parse_mmcif_atoms()` for Boltz-2 output
2. **Identical calibration results** - Fixed directory cleanup + seed parameter
3. **Python 3.9 compatibility** - Fixed type hints in constants.py

### Calibration Results (scFv-derived constructs from known antibodies)

| Binder | pTM | pLDDT | Contacts | Interface Area |
|--------|-----|-------|----------|----------------|
| Teplizumab-scFv | 0.944 | 96.4 | 42 | 2560 Å² |
| SP34-scFv | 0.753 | 91.7 | 30 | 2160 Å² |
| UCHT1-scFv | 0.784 | 94.1 | 32 | 2240 Å² |

**Important**: These are VH-linker-VL scFvs derived from the parental antibodies, NOT full IgG. See `docs/reference/calibration-methodology.md`.

### Key Findings

1. **Structural metrics ≠ affinity** - pTM/pLDDT/contacts are structural confidence, NOT binding affinity
2. **pDockQ not available** - It's an AlphaFold-Multimer metric, not native to Boltz-2
3. **Boltz-2 has IC50 prediction** - Can be enabled but not validated for antibodies

---

## What Actually Works

| Design Type | Method | Status |
|-------------|--------|--------|
| **VHH de novo** | BoltzGen `nanobody-anything` | ✅ Works |
| **scFv de novo** | BoltzGen `binder_length=250` | ❌ Broken (makes single-domain blob) |
| **Fab CDR redesign** | BoltzGen Fab scaffolds | 🔧 Not implemented (but supported) |
| **scFv from known Ab** | Optimization track | ✅ Works |

---

## Next Steps (Priority Order)

### 1. Fix Fab/scFv Design (HIGH)
- Implement BoltzGen Fab CDR redesign using human scaffolds
- See `example/fab_scaffolds/` and `example/fab_targets/` in BoltzGen repo
- Target CD3ε with proven Fab scaffolds

### 2. Affinity Prediction (MEDIUM)
- Enable Boltz-2 IC50 prediction (already available)
- Consider adding AttABseq (MIT) for sequence-based ΔΔG
- Consider Graphinity (BSD-3) for structure-based ΔΔG

### 3. Production Run (AFTER FIXES)
- Increase design counts to 100-200
- Run full pipeline
- Review candidates

---

## Output Files

```
data/outputs/
├── calibration.json          # Calibration thresholds
├── denovo/                   # BoltzGen designs (VHH only reliable)
├── optimized/                # Humanized variants (scFv works here)
├── structures/               # Boltz-2 predictions
├── filtered/                 # Filtered candidates
├── formatted/                # Bispecific constructs
└── reports/                  # Final reports
```

---

## Key Documentation

- `docs/reference/calibration-methodology.md` - Explains scFv constructs and metrics vs affinity
- `CLAUDE.md` - Implementation notes for Claude sessions
- `README.md` - Full project documentation

---

## BoltzGen Fab Scaffold Reference

From `github.com/HannesStark/boltzgen/blob/main/example/fab_scaffolds/`:

```yaml
# Structure of a Fab scaffold file
path: <antibody>.cif
include:
  - chain: { id: B, res_index: 1..121 }  # VH
  - chain: { id: A, res_index: 1..107 }  # VL
design:  # CDR regions only
  - chain: { id: B, res_index: <H-CDR1>,<H-CDR2>,<H-CDR3> }
  - chain: { id: A, res_index: <L-CDR1>,<L-CDR2>,<L-CDR3> }
design_insertions:  # Variable CDR lengths
  - insertion: { id: B, res_index: 99, num_residues: 3..21 }  # H-CDR3
```

This approach:
- Uses proven human frameworks (lower immunogenicity)
- Redesigns only CDRs (proper antibody structure)
- Supports variable CDR lengths
- Maintains VH/VL pairing
