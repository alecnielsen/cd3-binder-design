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

### Bug Fixes (2026-01-19 code review session)

- [x] **Critical: JSON serialization** - Fixed LiabilitySite objects causing TypeError on json.dump()
  - Changed to store positions as integers in CandidateScore
  - Files: `design_pipeline.py:350-396`, `scripts/05_filter_candidates.py:107-152`

- [x] **High: VL ignored in step-by-step structure prediction** - scripts/04_predict_structures.py now constructs scFv for VH/VL pairs
  - File: `scripts/04_predict_structures.py:80-126`

- [x] **High: Calibration used VH-only** - Now constructs scFv for paired antibodies during calibration
  - Files: `scripts/00_run_calibration.py:59-73`, `design_pipeline.py:99-113`

- [x] **Medium: Liability filtering was whole-sequence** - Now checks CDR-specific counts when `allow_*_cdr=False`
  - Added `cdr_deamidation_count`, etc. fields to CandidateScore
  - File: `filter_cascade.py:53-57, 213-236`

- [x] **Medium: De novo scFv treated as VHH** - Added scFv parsing to extract VH/VL from concatenated sequences
  - Added `parse_scfv()` and `is_likely_scfv()` utilities
  - Files: `src/utils/constants.py:82-136`, `src/formatting/__init__.py:86-94`

- [x] **Medium: Humanness filtering disabled in step-by-step** - Added `score_humanness_pair()` call
  - File: `scripts/05_filter_candidates.py:155-162`

- [x] **BioPhi position indexing** - `_apply_mutations()` now tries both 0-based and 1-based indexing
  - File: `src/design/optimization.py:262-295`

### Bug Fixes (2026-01-19 second code review)

- [x] **Critical: scFv treated as VHH** - Added scFv parsing in analysis/filtering to detect concatenated scFv
  - Files: `design_pipeline.py:324-340`, `scripts/05_filter_candidates.py:79-95`

- [x] **High: Contact counts were atomic** - Changed to residue-level contact counting
  - Files: `pdb_utils.py:186-236`, `modal/boltz2_app.py:91-107`

- [x] **Medium: Epitope config not wired** - InterfaceAnalyzer now accepts custom OKT3 residues from config
  - Files: `interface_analysis.py:92-104`, `design_pipeline.py:315-317`, `scripts/05_filter_candidates.py:71-73`

- [x] **Medium: HTML report showed raw enums** - Fixed filter result display with proper CSS classes
  - File: `report_generator.py:259-277, 295-310`

- [x] **Low: Aggregation filter was no-op** - Implemented aromatic content and consecutive aromatics check
  - File: `filter_cascade.py:264-300`

- [x] **Low: Fallback config unused** - Wired `relax_soft_filters_first` and `max_threshold_relaxation`
  - File: `filter_cascade.py:437-503`

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
