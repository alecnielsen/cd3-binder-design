# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-02-02)

### Pipeline: NEEDS FIXES BEFORE PRODUCTION

Trial run completed but revealed issues that must be resolved before scaling up.

```bash
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

python3 scripts/setup_fab_scaffolds.py             # Download Fab scaffolds (once)
python3 scripts/00_run_calibration.py              # ✅ Working
python3 scripts/02_run_denovo_design.py --config config.yaml  # ⚠️ VHH works, Fab produced 0 output
python3 scripts/03_run_optimization.py --config config.yaml   # ✅ Working (but NOT new designs)
python3 scripts/04_predict_structures.py --config config.yaml # ✅ Working
python3 scripts/05_filter_candidates.py --config config.yaml  # ⚠️ Bug: empty sequence candidate
python3 scripts/06_format_bispecifics.py --config config.yaml # ✅ Working
python3 scripts/07_generate_report.py --config config.yaml    # ✅ Working
```

---

## Design Types

| Type | Method | Output | Trial Run Status |
|------|--------|--------|------------------|
| **VHH** | BoltzGen `nanobody-anything` | Single-domain ~120 aa | ⚠️ 2 designs (expected 5) |
| **Fab** | BoltzGen `antibody-anything` CDR redesign | VH ~120 aa + VL ~107 aa | ❌ 0 designs |
| **scFv from known Ab** | Optimization track | VH-linker-VL | ✅ 9 variants |

**IMPORTANT**: The optimization track does NOT generate new designs. It reformats known antibody sequences (teplizumab, SP34, UCHT1) from `data/starting_sequences/*.yaml` as scFv for structure prediction. The only tracks producing novel binders are VHH and Fab de novo design.

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

## Trial Run Analysis (2026-02-02)

### What Worked

| Component | Evidence |
|-----------|----------|
| Calibration | All 3 known binders predicted with good pTM (0.75-0.94), contacts (30-42) |
| Structure prediction | Boltz-2 producing valid mmCIF outputs |
| Filtering cascade | 11 input → 7 first pass → 10 with fallback |
| Bispecific formatting | All 5 formats generated |
| Report generation | Scorecards created |

### Issues Found

| Issue | Severity | Details |
|-------|----------|---------|
| **Fab design: no output** | 🔴 High | Config had `num_fab_designs: 5` but 0 Fab designs were generated |
| **Empty sequence bug** | 🔴 High | Rank 1 candidate has `sequence: ""` - data flow issue |
| **Humanness scores null** | ⚠️ Medium | All `oasis_score_vh/vl` are null - BioPhi not returning data |
| **CDR-H3 length null** | ⚠️ Medium | ANARCI numbering may not be running |
| **VHH undercount** | ⚠️ Low | Got 2 VHH designs when config requested 5 |

---

## Recommendations Before Production Run

### 1. Fix Fab Design Track (BLOCKING)

The main de novo design approach (Fab CDR redesign) produced zero output. Must diagnose:

```bash
# Test Fab design directly on Modal
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb \
  --design-type fab --scaffolds adalimumab --num-designs 3
```

Check:
- Are scaffold files being loaded correctly?
- Is `run_boltzgen_fab()` being called in `02_run_denovo_design.py`?
- Are there errors in Modal logs?

### 2. Fix Empty Sequence Bug (BLOCKING)

The rank 1 candidate `"unknown"` has an empty sequence. Trace data flow in:
- `02_run_denovo_design.py` - is design output being parsed correctly?
- `05_filter_candidates.py` - is there a fallback creating empty candidates?

### 3. Fix Humanness Scoring (Important)

All OASis scores are null. Check:
- Is BioPhi installed and importable?
- `src/analysis/humanness.py` - is `score_humanness_pair()` being called?

### 4. Run Medium-Scale Test

Before 100+ designs, run with 10-20 designs to validate fixes:

```yaml
design:
  num_vhh_designs: 10
  num_fab_designs: 10
  fab_scaffolds:
    - adalimumab
    - belimumab
```

### 5. Future Enhancements (Non-blocking)

- Enable Boltz-2 IC50 prediction (already available, not antibody-validated)
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
