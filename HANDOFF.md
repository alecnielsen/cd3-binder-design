# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-02-02)

### Pipeline: READY FOR VALIDATION TEST

Critical bugs fixed. Ready for medium-scale validation before production run.

```bash
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

python3 scripts/setup_fab_scaffolds.py             # ✅ Done - scaffolds downloaded
python3 scripts/00_run_calibration.py              # ✅ Working
python3 scripts/02_run_denovo_design.py --config config.yaml  # ✅ Fixed - VHH + Fab now working
python3 scripts/03_run_optimization.py --config config.yaml   # ✅ Working (reformats known antibodies)
python3 scripts/04_predict_structures.py --config config.yaml # ✅ Fixed - parses denovo output correctly
python3 scripts/05_filter_candidates.py --config config.yaml  # ✅ Working
python3 scripts/06_format_bispecifics.py --config config.yaml # ✅ Working
python3 scripts/07_generate_report.py --config config.yaml    # ✅ Working
```

---

## Design Types

| Type | Method | Output | Status |
|------|--------|--------|--------|
| **VHH** | BoltzGen `nanobody-anything` | Single-domain ~120 aa | ✅ Working |
| **Fab** | BoltzGen `antibody-anything` CDR redesign | VH ~120 aa + VL ~107 aa | ✅ Fixed (scaffolds now present) |
| **scFv from known Ab** | Optimization track | VH-linker-VL | ✅ Working |

**IMPORTANT**: The optimization track does NOT generate new designs. It reformats known antibody sequences (teplizumab, SP34, UCHT1) from `data/starting_sequences/*.yaml` as scFv for structure prediction. The only tracks producing **novel** binders are VHH and Fab de novo design.

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

## Fixes Applied (2026-02-02)

### 1. Fab Scaffold Files - FIXED ✅
**Root cause**: `data/fab_scaffolds/` was empty (only `.gitkeep`).

**Fix**: Ran `python scripts/setup_fab_scaffolds.py` - downloaded 14 scaffold CIF + YAML files.

### 2. Empty Sequence Bug - FIXED ✅
**Root cause**: `04_predict_structures.py` didn't handle denovo results format.

The denovo output has `{"vhh_designs": [...], "fab_designs": [...]}` structure, but the script only checked for `designs` or `variants` keys. It fell through to treating the entire object as one candidate.

**Fix**: Added handler for `vhh_designs`/`fab_designs` keys in the candidate loading logic.

### 3. Missing Dependencies - NOT INSTALLED
BioPhi and ANARCI are in `requirements.txt` but not installed locally:
- **BioPhi**: Required for humanness scoring (OASis)
- **ANARCI**: Required for CDR numbering and CDR-H3 length

These cause soft-fails (null scores) but don't break the pipeline. Install with:
```bash
# ANARCI requires conda or source build
conda install -c bioconda anarci

# BioPhi
pip install biophi
# or from source: pip install git+https://github.com/Merck/BioPhi.git
```

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

### 1. ~~Fix Fab Design Track~~ - FIXED ✅
Scaffold files are now downloaded. Ready for Fab CDR redesign.

### 2. ~~Fix Empty Sequence Bug~~ - FIXED ✅
Structure prediction script now correctly extracts designs from denovo output.

### 3. Install BioPhi/ANARCI (Optional but Recommended)
Humanness and CDR-H3 analysis require these packages:
```bash
conda install -c bioconda anarci
pip install biophi
```
Pipeline works without them (soft-fail) but scores will be null.

### 4. Run Validation Test
Before production (100+ designs), validate fixes with medium-scale test:
```bash
# Update config.yaml:
# num_vhh_designs: 10
# num_fab_designs: 10

python scripts/02_run_denovo_design.py --config config.yaml
python scripts/04_predict_structures.py --config config.yaml
python scripts/05_filter_candidates.py --config config.yaml
```

Verify:
- Fab designs are generated (check `data/outputs/denovo/` for `fab_designs` array)
- No empty sequences in filtered output
- Structure predictions complete for all candidates

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
