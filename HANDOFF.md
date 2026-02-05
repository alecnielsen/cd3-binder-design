# Handoff Notes - CD3 Binder Design Pipeline

## Current Status (2026-02-05)

### Pipeline: BOTH TRACKS FULLY WORKING

VHH and Fab de novo design both validated. Latest run (Feb 5): **10 VHH + 9 Fab designs** generated successfully.

```bash
conda activate cd3-binder
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

python3 scripts/setup_fab_scaffolds.py             # ✅ Done - scaffolds downloaded
python3 scripts/00_run_calibration.py              # ✅ Working
python3 scripts/02_run_denovo_design.py --config config.yaml  # ✅ VHH + Fab working
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
| **Fab** | BoltzGen `antibody-anything` CDR redesign | VH ~120 aa + VL ~107 aa | ✅ Working (12 fully-human scaffolds) |
| **scFv from known Ab** | Optimization track | VH-linker-VL | ✅ Working |

**IMPORTANT**: The optimization track does NOT generate new designs. It reformats known antibody sequences (teplizumab, SP34, UCHT1) from `data/starting_sequences/*.yaml` as scFv for structure prediction. The only tracks producing **novel** binders are VHH and Fab de novo design.

---

## Understanding CIF Chain IDs (Important!)

**Why this matters**: BoltzGen uses CIF chain identifiers (`struct_asym.id`) which are different from PDB author chain IDs.

### CIF vs PDB Chain IDs

| Term | Where Used | Example | Meaning |
|------|------------|---------|---------|
| `struct_asym.id` | CIF files, BoltzGen | A, B, C, D | Sequential identifiers assigned to each polymer entity |
| `pdbx_strand_id` | RCSB API, PDB headers | H, L | Author-assigned labels (H=heavy, L=light) |

**Key insight**: The letters A, B, C, D are just sequential identifiers - they don't inherently mean "heavy" or "light" chain. You must check which entity each letter maps to.

### How Chain IDs Are Assigned in CIF Files

When a structure is deposited:
1. Each unique polymer (protein chain) gets an entity ID (1, 2, 3...)
2. Each entity gets a sequential chain label: entity 1 → A, entity 2 → B, etc.
3. If a structure has 2 copies of the Fab, chains might be A,B,C,D (2 copies of VL and VH)

### Example: Adalimumab (6CR1)

```
Entity 1: Light chain (VL-IgE fusion) → Chain A in CIF
Entity 2: Heavy chain (VH-IgE fusion) → Chain B in CIF
Entity 3-5: Ligands/water → Chains C, D, E...
```

So for adalimumab: `vh_chain: "B"`, `vl_chain: "A"`

### All Scaffold Chain Mappings (Verified)

12 fully-human scaffolds (crenezumab and mab1 excluded - see Scaffold Selection Rationale):

| Scaffold | PDB | VH Chain | VL Chain | Notes |
|----------|-----|----------|----------|-------|
| adalimumab | 6cr1 | B | A | VL=entity1, VH=entity2 |
| belimumab | 5y9k | B | A | VL=entity1, VH=entity2 |
| dupilumab | 6wgb | A | B | VH=entity1, VL=entity2 |
| golimumab | 5yoy | G | D | VL=entity2 (D,E,F), VH=entity3 (G,H,I); TNF=entity1 |
| guselkumab | 4m6m | B | A | VL=entity1, VH=entity2 |
| necitumumab | 6b3s | B | C | EGFR=entity1, VH=entity2, VL=entity3 |
| nirsevimab | 5udc | A | B | VH=entity1, VL=entity2 |
| sarilumab | 8iow | D | C | IL6R=entity1, VL=entity2, VH=entity3 |
| secukinumab | 6wio | A | B | VH=entity1, VL=entity2, IL17A=entity3 |
| tezepelumab | 5j13 | C | B | TSLP=entity1, VL=entity2, VH=entity3 |
| tralokinumab | 5l6y | B | C | IL13=entity1, VH=entity2, VL=entity3 |
| ustekinumab | 3hmw | B | A | VL=entity1, VH=entity2 |

### How to Verify Chain IDs for a New Scaffold

```bash
# Download CIF and check struct_asym section
grep "_struct_asym" data/fab_scaffolds/scaffold.cif -A 10 | grep -E "^[A-Z]"
# Output: A N N 1 ?   <- Chain A = Entity 1
#         B N N 2 ?   <- Chain B = Entity 2

# Check what each entity is
grep -A 200 "^_entity_poly.entity_id" scaffold.cif | head -30
# Look for "heavy chain" or "light chain" in descriptions
```

Or use the verification script: `python scripts/verify_fab_chains.py`

---

## Understanding ipTM Scores

**ipTM (interface predicted Template Modeling)** is the primary quality metric from BoltzGen/Boltz-2 for de novo designs.

| ipTM Range | Interpretation |
|------------|----------------|
| < 0.2 | Low confidence - interface may not form as predicted |
| **0.2 - 0.4** | **Typical for de novo designs** - reasonable confidence |
| > 0.4 | High confidence - unusual for de novo, common for known complexes |

**Critical caveat**: ipTM is a structural confidence score, NOT an affinity predictor. A design with ipTM 0.25 could have better actual affinity than one with 0.35. Experimental validation (SPR/BLI) is required.

---

## Latest Run (2026-02-05)

### De Novo Design Results

| Type | Requested | Generated | Best ipTM | Design ID |
|------|-----------|-----------|-----------|-----------|
| **VHH** | 10 | 10 | 0.368 | vhh_1SY6_0000 |
| **Fab** | 10 | 9 | 0.368 | fab_1SY6_0002 |

Top designs by ipTM:
- fab_1SY6_0002: 0.368 (pTM 0.697)
- vhh_1SY6_0000: 0.368
- vhh_1XIW_0001: 0.361
- vhh_1SY6_0003: 0.317
- vhh_1XIW_0003: 0.305

### Key Observations

1. **Fab track now working** - 9/10 Fab designs generated (1 failed RMSD filter)
2. **ipTM scores competitive** - Best designs at 0.368 (high end of typical range)
3. **VHH and Fab comparable** - Both tracks producing similar quality designs

---

## Pipeline Architecture Note

**Why does step 04 (predict_structures) exist if BoltzGen already does structure prediction?**

BoltzGen outputs ipTM/pTM during design, but the filtering step requires additional metrics:
- `interface_area` - buried surface area at interface
- `num_contacts` - residue-level contact count
- `interface_residues_target` - for epitope classification

Step 04 re-runs Boltz-2 to extract these interface metrics. It also processes the optimization track (known antibody scFvs) which weren't designed via BoltzGen.

---

## Earlier Validation Results (2026-02-03)

### VHH Pipeline Run Summary

| Step | Status | Output |
|------|--------|--------|
| De novo design | ✅ | 10 VHH designs generated |
| Structure prediction | ✅ | 22 candidates predicted (10 new VHH + 2 old + 10 optimized) |
| Filtering | ✅ | 10 candidates passed (via fallback relaxation) |
| Bispecific formatting | ✅ | 29 constructs across 5 formats |
| Report generation | ✅ | 12 report files |

### Key Observations (Feb 3)

1. **VHH design reliably produces 10 designs** - Confirmed across 4 consecutive runs
2. **Fab RMSD filter is strict** - Designs may fail 2.5 Å threshold; this is legitimate quality filtering
3. **pDockQ shows 0.000** - Expected; pDockQ is NOT a native Boltz-2 metric. Use pTM/ipTM instead.
4. **Fallback filtering works** - 0 candidates passed strict thresholds, 10 passed after relaxation
5. **De novo VHH designs score competitively** - Within range of known antibody scFvs

---

## Fixes Applied

### Session 2026-02-04 (Latest)

1. **BoltzGen refolding fix** - Upgraded from PyPI v0.2.0 to GitHub main branch which includes Dec 17, 2025 refolding fix
2. **Modal image updated** - Now uses `git+https://github.com/HannesStark/boltzgen.git@main` instead of PyPI package
3. **Full pipeline now works** - All 5 BoltzGen steps (design, inverse_folding, folding, analysis, filtering) complete successfully
4. **Scaffold YAML format fixed** - Updated all 14 scaffolds to official BoltzGen format with `structure_groups`, `exclude`, `design_insertions` (all 6 CDRs), and `reset_res_index` sections
5. **Fab RMSD now passing** - With correct YAML format: 5/20 designs pass 2.5 Å threshold (~25% pass rate). Before fix: 0/20 passed (18-22 Å RMSD)

### Session 2026-02-03

1. **Fab scaffold chain IDs fixed** - Updated all 14 scaffolds with correct CIF chain labels
2. **VH/VL parsing fixed** - Modal app now identifies chains by sequence patterns (EVQL... for VH, DIQM... for VL) instead of assuming entity order
3. **Fab track enabled** - Set `num_fab_designs: 10` in config.yaml

### Session 2026-02-03 (Earlier)

1. **Modal package installation** - Added to conda environment (was missing)
2. **BoltzGen Modal deployment** - Deployed app + downloaded model weights
3. **Path doubling bug** - Fixed in `modal/boltzgen_app.py`
4. **VHH-only validation** - Set `num_fab_designs: 0` to bypass Fab issues

### Session 2026-02-02

1. **Fab Scaffold Files** - Ran `setup_fab_scaffolds.py` to download CIF + YAML files
2. **Empty Sequence Bug** - Added handler for `vhh_designs`/`fab_designs` keys
3. **Conda Environment** - Created `cd3-binder` with ANARCI + Sapiens

---

## Production Config

Both tracks are now ready for production:

```yaml
design:
  num_vhh_designs: 50
  num_fab_designs: 50
  fab_scaffolds:
    - adalimumab
    - belimumab
    - dupilumab
    - golimumab
    - guselkumab
    - necitumumab
    - nirsevimab
    - sarilumab
    - secukinumab
    - tezepelumab
    - tralokinumab
    - ustekinumab
```

### Scaffold Selection Rationale

**12 fully-human scaffolds** are used. Two were excluded:
- **crenezumab**: Humanized (not fully human) - contains murine back-mutation residues in framework regions
- **mab1 (3H42)**: Unclear provenance - PDB structure lacks documentation of human/humanized status

**Why fully-human matters**:
- BoltzGen replaces CDRs but keeps framework regions
- Humanized frameworks retain ~5-10% murine residues (back-mutations for CDR conformation)
- Fully-human frameworks = 100% human germline = lower immunogenicity risk

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
│   ├── denovo_results_20260203_063803.json  # Validation (10 VHH)
│   └── fab_test_designs.json                # Fab test (2 designs)
├── optimized/                # Known antibody scFvs
├── structures/               # Boltz-2 predictions
├── filtered/                 # Filtered candidates
├── formatted/                # Bispecific constructs
└── reports/                  # Scorecards and HTML report
```

---

## Key Documentation

- `docs/reference/calibration-methodology.md` - Explains scFv constructs and metrics vs affinity
- `docs/reference/implementation.md` - Code examples and CIF chain ID explanation
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
