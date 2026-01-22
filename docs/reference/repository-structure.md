# Repository Structure

```
cd3-binder-design/
├── README.md                            # Main documentation
├── CLAUDE.md                            # AI assistant context
├── pyproject.toml                       # Python package configuration
├── requirements.txt                     # Dependencies
├── config.yaml                          # Default pipeline configuration
│
├── data/
│   ├── targets/
│   │   ├── 1XIW.pdb                     # CD3εδ structure for design
│   │   ├── 1SY6.pdb                     # CD3εγ + OKT3 complex
│   │   └── README.md                    # Structure preparation notes
│   │
│   ├── starting_sequences/
│   │   ├── teplizumab.yaml              # hOKT3γ1(Ala-Ala) VH/VL
│   │   ├── sp34.yaml                    # SP34 VH/VL
│   │   ├── ucht1.yaml                   # UCHT1 VH/VL
│   │   ├── affinity_mutations.yaml      # Literature affinity-reducing mutations
│   │   └── README.md                    # Sequence sources
│   │
│   ├── frameworks/
│   │   ├── igg1_constant_regions.yaml   # CH1, CH2, CH3, hinge sequences
│   │   ├── knob_hole_mutations.yaml     # Fc engineering mutations
│   │   ├── linkers.yaml                 # scFv and fusion linkers
│   │   └── placeholder_targets.yaml     # Trastuzumab (HER2) sequences
│   │
│   └── outputs/                         # Generated during pipeline runs
│       ├── denovo/                      # BoltzGen designs
│       ├── optimized/                   # Optimized existing binders
│       ├── structures/                  # Predicted structures
│       ├── filtered/                    # Post-filtering candidates
│       └── formatted/                   # Final bispecific sequences
│
├── src/
│   ├── __init__.py
│   │
│   ├── design/
│   │   ├── __init__.py
│   │   ├── boltzgen_runner.py           # BoltzGen interface for Modal
│   │   ├── denovo_design.py             # De novo VHH/scFv generation
│   │   ├── optimization.py              # Mutagenesis of existing binders
│   │   └── affinity_variants.py         # Generate affinity-attenuated panel
│   │
│   ├── structure/
│   │   ├── __init__.py
│   │   ├── abodybuilder.py              # ABodyBuilder2 wrapper
│   │   ├── boltz_complex.py             # Boltz-2 for complex prediction
│   │   ├── pdb_utils.py                 # PDB parsing and manipulation
│   │   └── interface_analysis.py        # Binding interface metrics
│   │
│   ├── analysis/
│   │   ├── __init__.py
│   │   ├── numbering.py                 # ANARCI wrapper
│   │   ├── humanness.py                 # BioPhi/Sapiens scoring
│   │   ├── liabilities.py               # Sequence liability detection
│   │   └── developability.py            # Combined developability scoring
│   │
│   ├── formatting/
│   │   ├── __init__.py
│   │   ├── base.py                      # Base bispecific formatter
│   │   ├── crossmab.py                  # CrossMab Fab×Fab assembly
│   │   ├── fab_scfv.py                  # Asymmetric Fab + scFv
│   │   ├── fab_vhh.py                   # Asymmetric Fab + VHH
│   │   ├── igg_scfv.py                  # Morrison IgG-(scFv)₂
│   │   └── igg_vhh.py                   # Morrison IgG-(VHH)₂
│   │
│   ├── pipeline/
│   │   ├── __init__.py
│   │   ├── config.py                    # Pipeline configuration
│   │   ├── design_pipeline.py           # End-to-end orchestration
│   │   ├── filter_cascade.py            # Multi-stage filtering
│   │   └── report_generator.py          # Developability scorecards
│   │
│   └── utils/
│       ├── __init__.py
│       └── constants.py                 # Amino acid properties, scFv parsing
│
├── modal/
│   ├── __init__.py
│   ├── boltzgen_app.py                  # BoltzGen Modal deployment
│   ├── boltz2_app.py                    # Boltz-2 Modal deployment
│   └── abodybuilder_app.py              # ABodyBuilder2 Modal deployment
│
├── scripts/
│   ├── 00_run_calibration.py            # Calibrate thresholds
│   ├── 01_setup_targets.py              # Download structures
│   ├── 02_run_denovo_design.py          # Run BoltzGen on Modal
│   ├── 03_run_optimization.py           # Generate optimized variants
│   ├── 04_predict_structures.py         # Run structure prediction
│   ├── 05_filter_candidates.py          # Apply filtering cascade
│   ├── 06_format_bispecifics.py         # Convert to bispecific formats
│   ├── 07_generate_report.py            # Generate final report
│   └── run_full_pipeline.py             # Execute all steps
│
├── tests/
│   ├── __init__.py
│   ├── test_scfv_offsets.py             # scFv parsing tests
│   └── test_filter_aggregation.py       # Aggregation filter tests
│
└── docs/
    └── reference/
        ├── bispecific-formats.md        # Format diagrams
        ├── implementation.md            # Code examples
        └── repository-structure.md      # This file
```

## Module Summaries

### src/analysis/
Sequence analysis tools: liability detection (deamidation, glycosylation), ANARCI numbering, BioPhi humanness scoring, developability assessment.

### src/design/
Design generation: BoltzGen runner for de novo design, optimization of existing binders, affinity variant generation.

### src/formatting/
Bispecific assembly: 5 formats (CrossMab, Fab+scFv, Fab+VHH, IgG-scFv, IgG-VHH) with linker insertion and constant region grafting.

### src/pipeline/
Orchestration: configuration loading, multi-stage filtering cascade, end-to-end pipeline, HTML/JSON report generation.

### src/structure/
Structure tools: PDB parsing, interface analysis, ABodyBuilder2 and Boltz-2 wrappers.

### modal/
GPU deployments for Modal cloud compute. Required for BoltzGen and Boltz-2 (CUDA only).

### scripts/
Pipeline execution scripts numbered 00-07 for step-by-step or full pipeline runs.
