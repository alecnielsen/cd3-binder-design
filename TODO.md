# CD3 Binder Design Pipeline - Status

Last updated: 2026-01-19

## Completed

- [x] Directory structure and placeholder files
- [x] pyproject.toml and requirements.txt
- [x] CLAUDE.md for AI context
- [x] Starting sequences (teplizumab, SP34, UCHT1, affinity mutations)
- [x] Framework sequences (IgG1 constant regions, linkers, knob-hole mutations)
- [x] src/analysis/ modules (liabilities, numbering, humanness, developability)
- [x] src/formatting/ modules (all 5 bispecific formats: CrossMab, Fab+scFv, Fab+VHH, IgG-scFv, IgG-VHH)
- [x] src/design/ modules (BoltzGen runner, denovo_design, optimization, affinity_variants)
- [x] src/structure/ modules (ABodyBuilder2, Boltz-2, pdb_utils, interface_analysis)
- [x] src/pipeline/ modules (config, filter_cascade, design_pipeline, report_generator)
- [x] Modal deployments (boltzgen_app, boltz2_app, abodybuilder_app)
- [x] Pipeline scripts (00-07 + run_full_pipeline.py)
- [x] README with Known Limitations, Experimental Validation Requirements, Reproducibility sections
- [x] Calibration phase for filter threshold validation
- [x] Epitope annotation comparing designs to OKT3 epitope
- [x] Fallback handling for insufficient candidates
- [x] Default config.yaml at repo root
- [x] Fixed API mismatches (non-existent classes/methods in design_pipeline.py, optimization.py, filter script)
- [x] Fixed config loader to support nested schema from README (backwards compatible)
- [x] Fixed structure prediction to handle VH/VL pairs (creates scFv with linker)
- [x] Fixed f-string format specifier errors in report_generator.py and filter script
- [x] Fixed calibration to use configurable margins from config
- [x] Fixed de novo design count truncation (even distribution with remainder)

## Optional / Future Work

- [ ] Unit tests for core modules (tests/ directory exists but is empty)
- [ ] Example Jupyter notebooks (notebooks/ mentioned in README)
- [ ] Update src/analysis/__init__.py with full exports

## Architecture Notes

### Key Design Decisions
- **Affinity target**: ~50 nM Kd hypothesis (cannot be computationally enforced - see Known Limitations)
- **Licensing**: All permissive (MIT, BSD, Apache, CC-BY) - excluded IgFold, NetMHCIIpan
- **Compute**: Modal for GPU (BoltzGen/Boltz-2 require CUDA)
- **Calibration**: Uses known binders to set filter thresholds before filtering de novo designs

### File Structure
```
src/
├── analysis/     # Liability detection, numbering, humanness, developability
├── design/       # BoltzGen runner, optimization, affinity variants
├── formatting/   # 5 bispecific formats
├── pipeline/     # Config, filtering, orchestration, reports
├── structure/    # ABodyBuilder2, Boltz-2, PDB utils, interface analysis
└── utils/        # Constants

modal/            # GPU compute deployments
scripts/          # Pipeline execution (00-07)
data/             # Sequences, frameworks, outputs
```

### Running the Pipeline
```bash
# Full pipeline
python scripts/run_full_pipeline.py --config config.yaml

# Or step by step
python scripts/00_run_calibration.py    # Set filter thresholds
python scripts/01_setup_targets.py      # Download PDB structures
python scripts/02_run_denovo_design.py  # BoltzGen
python scripts/03_run_optimization.py   # Humanization + affinity variants
python scripts/04_predict_structures.py # Boltz-2
python scripts/05_filter_candidates.py  # Apply filters
python scripts/06_format_bispecifics.py # Generate formats
python scripts/07_generate_report.py    # HTML/JSON reports
```
