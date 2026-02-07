# Pipeline Scripts

This directory contains the execution scripts for the CD3 binder design pipeline. Each script corresponds to a specific step in the pipeline and can be run individually or orchestrated via `run_full_pipeline.py`.

## Prerequisites

Before running the pipeline, set up Fab scaffolds:

```bash
python scripts/setup_fab_scaffolds.py
```

This downloads 14 human antibody scaffold structures from RCSB PDB and creates BoltzGen YAML specifications.

## Pipeline Order

The scripts are numbered to indicate their execution order:

| Step | Script | Description |
|------|--------|-------------|
| - | `setup_fab_scaffolds.py` | Download Fab scaffolds from PDB (run once) |
| 0 | `00_run_calibration.py` | Calibrate filter thresholds using known binders |
| 1 | `01_setup_targets.py` | Download and prepare CD3 target structures |
| 2 | `02_run_denovo_design.py` | Run BoltzGen VHH + Fab design on Modal |
| 3 | `03_run_optimization.py` | Generate optimized variants of existing binders |
| 4 | `04_predict_structures.py` | Run Boltz-2 complex structure prediction |
| 5 | `05_filter_candidates.py` | Apply filtering cascade (binding, humanness, liabilities) |
| 5b | `05b_validate_candidates.py` | Affinity scoring (ProteinMPNN, AntiFold) + Protenix cross-validation |
| 6 | `06_format_bispecifics.py` | Convert candidates to bispecific antibody formats |
| 7 | `07_generate_report.py` | Generate final HTML and JSON reports |

## Running the Pipeline

### Full Pipeline

```bash
# Setup (once)
python scripts/setup_fab_scaffolds.py

# Run all steps
python scripts/run_full_pipeline.py --config config.yaml

# Run without Modal (local/mock mode)
python scripts/run_full_pipeline.py --config config.yaml --no-modal

# Skip calibration (if already done)
python scripts/run_full_pipeline.py --config config.yaml --skip-calibration

# Start from a specific step
python scripts/run_full_pipeline.py --config config.yaml --start-from 3
```

### Individual Steps

```bash
# Prerequisites: Setup scaffolds (run once)
python scripts/setup_fab_scaffolds.py

# Step 0: Calibration (RUN FIRST)
python scripts/00_run_calibration.py --config config.yaml

# Step 1: Setup targets (does not require config)
python scripts/01_setup_targets.py --output-dir data/targets

# Step 2: De novo design - VHH + Fab (requires Modal)
python scripts/02_run_denovo_design.py --config config.yaml

# etc.
```

## Scientific Assumptions

### Calibration (Step 0)

**Critical**: Calibration must run before filtering to set appropriate thresholds.

- Uses known binders (teplizumab, SP34, UCHT1) as positive controls
- For paired antibodies (VH+VL), constructs scFv (VH-linker-VL) for accurate structure prediction
- Runs Boltz-2 on known binders to establish baseline metrics
- Sets filter thresholds to accept all known binders with margin:
  - `min_pdockq = min(known_binder_pdockq) - 0.05`
  - `min_interface_area = min(known_binder_area) - 100`
  - `min_contacts = min(known_binder_contacts) - 2`

**Assumption**: If known binders fail default thresholds, the thresholds are wrong, not the binders.

### Target Setup (Step 1)

- Downloads PDB structures 1XIW (CD3ed) and 1SY6 (CD3eg + OKT3)
- Extracts CD3 epsilon chain for design targeting
- Validates structures are suitable for BoltzGen input

### De Novo Design (Step 2)

**Assumption**: BoltzGen can generate VHH/Fab sequences that bind CD3 epsilon.

Generates two types of binders:

| Type | Protocol | Output | Description |
|------|----------|--------|-------------|
| **VHH** | `nanobody-anything` | ~120 aa single-domain | De novo nanobody |
| **Fab** | `antibody-anything` | VH + VL (~227 aa total) | CDR redesign on human scaffolds |

**Fab CDR redesign** uses 14 proven human scaffolds (adalimumab, belimumab, etc.) and redesigns only the CDR loops, maintaining human frameworks for lower immunogenicity.

- Runs on Modal with A100 GPU
- Generates configurable number of designs (default: 200 VHH + 200 Fab)
- Uses random seeds for reproducibility
- Target structures from step 1 guide design generation
- Fab design requires scaffold files (run `setup_fab_scaffolds.py` first)

**Limitation**: BoltzGen hit rate is unknown for CD3. May require 500-1000 designs for sufficient diversity.

### Optimization (Step 3)

- Loads starting sequences (teplizumab, SP34, UCHT1)
- Applies humanization via BioPhi/Sapiens
- Generates affinity-attenuated variants using literature mutations
- Creates combinatorial variant panel

**Assumption**: Humanization preserves binding activity (30-70% success rate typical).

### Structure Prediction (Step 4)

**Assumption**: Boltz-2 can predict binder-CD3 complex structures accurately enough for filtering.

- Runs on Modal with A100 GPU
- Aggregates candidates from both denovo and optimized directories
- For paired antibodies (VH+VL), constructs scFv (VH-linker-VL) for structure prediction
- Predicts complex for each candidate + CD3 epsilon
- Extracts pDockQ, interface area, contact residues
- Annotates epitope overlap with OKT3

**Critical**: pDockQ is a structural confidence score, NOT an affinity predictor.

### Filtering (Step 5)

Applies filters in order:
1. **Binding quality**: pDockQ, interface area, contacts (calibrated thresholds)
2. **Humanness**: OASis score > 0.8
3. **Sequence liabilities**: Deamidation, isomerization, glycosylation in CDRs
4. **Developability**: CDR-H3 length, net charge, pI, hydrophobic patches
5. **Aggregation**: Aromatic content, consecutive aromatics

**Fallback logic**: If < 10 candidates survive:
1. Relax soft filters (oxidation, hydrophobic patches)
2. Relax thresholds by up to 10%
3. Include borderline candidates with risk flags

### Candidate Validation (Step 5b)

Runs after filtering on the top N candidates. **Results are informational only** — not used for filtering or ranking.

Three validation analyses:
1. **ProteinMPNN log-likelihood** (MIT): Inverse folding affinity proxy, best-validated for de novo antibodies (Spearman r=0.27-0.41). Runs locally on CPU.
2. **AntiFold log-likelihood** (BSD-3): Antibody-specific inverse folding with nanobody support. Runs locally on CPU.
3. **Protenix re-prediction** (Apache 2.0): Cross-validates Boltz-2 structure predictions. Runs on Modal H100.

Cross-validation flags candidates where Boltz-2 and Protenix ipTM disagree by >0.1.

**Output**: `data/outputs/validated/validated_candidates.json` — same as filtered_candidates.json with added `validation` section per candidate.

**Prerequisite**: Deploy Protenix on Modal (`modal deploy modal/protenix_app.py`). ProteinMPNN and AntiFold are optional (`pip install proteinmpnn antifold`); missing tools produce error messages but don't fail the step.

### Bispecific Formatting (Step 6)

Converts candidates into 5 bispecific formats:
- CrossMab (Fab x Fab with CH1-CL swap)
- Asymmetric Fab + scFv (knob-in-hole)
- Asymmetric Fab + VHH (knob-in-hole)
- IgG-(scFv)2 Morrison (symmetric)
- IgG-(VHH)2 Morrison (symmetric)

**Assumption**: Placeholder target arm (trastuzumab) is configurable.

### Report Generation (Step 7)

- Generates HTML and JSON developability scorecards
- Includes ranked candidate table with metrics
- Documents all filter relaxations and risk flags
- Adds provenance metadata for reproducibility

## Dependencies

- Python 3.10+
- Modal account (for GPU compute in steps 2, 4)
- Dependencies in `requirements.txt`

## Error Handling

Each script:
- Validates config and input files before processing
- Returns non-zero exit code on failure
- Logs errors with context for debugging

The full pipeline aborts on any step failure with a clear error message.

## Outputs

| Step | Output Directory |
|------|-----------------|
| 0 | `data/outputs/calibration.json` |
| 1 | `data/targets/` |
| 2 | `data/outputs/denovo/` |
| 3 | `data/outputs/optimized/` |
| 4 | `data/outputs/structures/` |
| 5 | `data/outputs/filtered/` |
| 5b | `data/outputs/validated/` |
| 6 | `data/outputs/formatted/` |
| 7 | `data/outputs/reports/` |
