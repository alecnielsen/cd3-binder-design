# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-01-26)

### Pipeline: WORKING

All pipeline steps run successfully with both VHH and Fab design tracks.

```bash
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

python3 scripts/setup_fab_scaffolds.py             # Download Fab scaffolds (once)
python3 scripts/00_run_calibration.py              # ✅ Working
python3 scripts/02_run_denovo_design.py --config config.yaml  # ✅ Working (VHH + Fab)
python3 scripts/03_run_optimization.py --config config.yaml   # ✅ Working
python3 scripts/04_predict_structures.py --config config.yaml # ✅ Working
python3 scripts/05_filter_candidates.py --config config.yaml  # ✅ Working
python3 scripts/06_format_bispecifics.py --config config.yaml # ✅ Working
python3 scripts/07_generate_report.py --config config.yaml    # ✅ Working
```

---

## Design Types

| Type | Method | Output | Status |
|------|--------|--------|--------|
| **VHH** | BoltzGen `nanobody-anything` | Single-domain ~120 aa | ✅ Working |
| **Fab** | BoltzGen `antibody-anything` CDR redesign | VH ~120 aa + VL ~107 aa | ✅ Working |
| **scFv from known Ab** | Optimization track | VH-linker-VL | ✅ Working |

### Fab CDR Redesign

The Fab design uses BoltzGen's `antibody-anything` protocol with human antibody scaffolds. Only the CDR loops are redesigned while the human framework regions are preserved (lower immunogenicity).

**14 proven human Fab scaffolds available:**
adalimumab, belimumab, crenezumab, dupilumab, golimumab, guselkumab, mab1, necitumumab, nirsevimab, sarilumab, secukinumab, tezepelumab, tralokinumab, ustekinumab

Configure scaffolds in `config.yaml`:
```yaml
design:
  num_fab_designs: 5
  fab_scaffolds:
    - adalimumab
    - belimumab
    - dupilumab
```

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

## Next Steps

### 1. Production Run
- Increase design counts to 100-200
- Run full pipeline
- Review candidates

### 2. Affinity Prediction (Future Enhancement)
- Enable Boltz-2 IC50 prediction (already available)
- Consider adding AttABseq (MIT) for sequence-based ΔΔG
- Consider Graphinity (BSD-3) for structure-based ΔΔG

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

## BoltzGen Fab Implementation Details

### Scaffold Files

Located in `data/fab_scaffolds/`:
- `*.cif` - Structure files downloaded from RCSB PDB
- `*.yaml` - BoltzGen scaffold specs with CDR definitions

Generate with:
```bash
python scripts/setup_fab_scaffolds.py
```

### Scaffold YAML Format

```yaml
# Example: adalimumab.6cr1.yaml
path: adalimumab.6cr1.cif
include:
  - chain: { id: B, res_index: 1..121 }  # VH
  - chain: { id: A, res_index: 1..107 }  # VL
design:  # CDR regions only
  - chain: { id: B, res_index: 26..32,52..57,99..110 }  # H-CDRs
  - chain: { id: A, res_index: 24..34,50..56,89..97 }   # L-CDRs
design_insertions:  # Variable CDR lengths
  - insertion: { id: B, res_index: 99, num_residues: 3..21 }  # H-CDR3
  - insertion: { id: A, res_index: 89, num_residues: 3..12 }  # L-CDR3
```

This approach:
- Uses proven human frameworks (lower immunogenicity)
- Redesigns only CDRs (proper antibody structure)
- Supports variable CDR lengths
- Maintains VH/VL pairing
- Outputs separate VH and VL sequences (downstream formatting joins with linker)
