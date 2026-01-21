# Claude Project Context: CD3 Binder Design

## Project Overview

This is a computational pipeline for designing CD3-binding domains (VHH nanobodies and scFvs) for use in bispecific antibody therapeutics. The pipeline uses open-source, permissively-licensed tools to generate, filter, and format CD3 binders with drug-like properties.

**Goal**: Design ~10 CD3 binders for experimental testing in cell lines and animal models.

## Key Design Decisions

### Target Affinity: ~50 nM Kd
- High-affinity CD3 binders (< 1 nM) cause more cytokine release syndrome (CRS)
- Literature supports ~50-100 nM range for balancing efficacy with safety
- We include affinity-attenuated variants (10x, 100x weaker than wild-type)

### Licensing: All Permissive
All tools must have MIT, BSD, Apache, or CC-BY licenses for commercial use:
- **Use**: BoltzGen (MIT), Boltz-2 (MIT), ABodyBuilder2 (BSD), ANARCI (BSD), BioPhi (MIT), AbLang (CC-BY)
- **Avoid**: IgFold (JHU Academic - non-commercial only), NetMHCIIpan (DTU academic)

### Bispecific Formats
Five formats are supported:
1. CrossMab (Fab × Fab with CH1-CL swap)
2. Asymmetric Fab + scFv (knob-in-hole)
3. Asymmetric Fab + VHH (knob-in-hole)
4. IgG-(scFv)₂ Morrison (symmetric, 2x CD3)
5. IgG-(VHH)₂ Morrison (symmetric, 2x CD3)

### Compute: Modal for GPU
BoltzGen and Boltz-2 require CUDA GPUs. Modal provides on-demand A100/H100 access.

## Repository Structure

```
cd3-binder-design/
├── data/
│   ├── targets/              # CD3 structures (1XIW, 1SY6)
│   ├── starting_sequences/   # Teplizumab, SP34, UCHT1
│   ├── frameworks/           # IgG constant regions, linkers
│   └── outputs/              # Generated designs
├── src/
│   ├── design/               # BoltzGen runner, optimization
│   ├── structure/            # ABodyBuilder2, Boltz-2 wrappers
│   ├── analysis/             # Numbering, humanness, liabilities
│   ├── formatting/           # Bispecific assembly
│   ├── pipeline/             # Orchestration
│   └── utils/                # Helpers
├── modal/                    # Modal GPU deployments
├── scripts/                  # Pipeline execution scripts (01-07)
└── tests/
```

## Key Files

### Data Files
- `data/starting_sequences/teplizumab.yaml` - Humanized OKT3 VH/VL sequences
- `data/starting_sequences/affinity_mutations.yaml` - Literature-known affinity variants
- `data/frameworks/igg1_constant_regions.yaml` - CH1, CH2, CH3, hinge, CL sequences
- `data/frameworks/knob_hole_mutations.yaml` - T366W (knob), T366S/L368A/Y407V (hole)

### Core Modules
- `src/analysis/liabilities.py` - Deamidation, glycosylation, oxidation detection
- `src/analysis/numbering.py` - ANARCI wrapper for CDR identification
- `src/analysis/humanness.py` - BioPhi/OASis humanness scoring
- `src/formatting/crossmab.py` - CrossMab assembly with CH1-CL swap
- `src/formatting/igg_scfv.py` - Morrison format (symmetric 2x scFv)

### Modal Deployments
- `modal/boltzgen_app.py` - BoltzGen for de novo design
- `modal/boltz2_app.py` - Boltz-2 for complex prediction

## Common Tasks

### Running the Pipeline
```bash
python scripts/run_full_pipeline.py --config config.yaml
```

### Running Individual Steps
```bash
python scripts/01_setup_targets.py      # Download structures
python scripts/02_run_denovo_design.py  # BoltzGen (Modal)
python scripts/05_filter_candidates.py  # Apply filters
python scripts/06_format_bispecifics.py # Generate formats
```

### Adding a New Bispecific Format
1. Create `src/formatting/new_format.py`
2. Implement assembly function following `base.py` interface
3. Register in `src/formatting/__init__.py`
4. Add to `config.yaml` formats list

### Testing Liability Detection
```python
from src.analysis.liabilities import LiabilityScanner
scanner = LiabilityScanner()
report = scanner.scan("EVQLVESGGGLVQ...")  # Returns LiabilityReport object
print(report.deamidation_sites)  # List of LiabilitySite objects
print(report.cdr_liabilities)    # Count of CDR liabilities
```

### Humanness Scoring
```python
from src.analysis.humanness import score_humanness_pair
report = score_humanness_pair(vh_sequence, vl_sequence)  # vl optional for VHH
print(report.mean_score)         # Mean OASis score
print(report.vh_report.oasis_score)
```

### Developability Assessment
```python
from src.analysis.developability import DevelopabilityAssessor
assessor = DevelopabilityAssessor()
report = assessor.assess(vh_sequence, vl_sequence)  # vl optional for VHH
print(report.physicochemical.net_charge)
print(report.aggregation.hydrophobic_patches)
```

## Important Constants

### Filtering Thresholds (default)
```python
FILTER_THRESHOLDS = {
    "min_pdockq": 0.5,
    "min_oasis_score": 0.8,
    "cdr_h3_length_range": (8, 20),
    "net_charge_range": (-2, 4),
    "max_hydrophobic_patches": 2,
}
```

### Liability Motifs
```python
DEAMIDATION = ['NG', 'NS', 'NT', 'ND', 'NH']  # N followed by small/flexible residue
ISOMERIZATION = ['DG', 'DS', 'DT', 'DD', 'DH', 'DN']  # D followed by small/flexible residue
GLYCOSYLATION = r'N[^P][ST]'  # N-X-S/T where X != P
```

### Linker Sequences
```python
SCFV_LINKER = "GGGGSGGGGSGGGGS"  # (G4S)₃
FC_FUSION_LINKER = "GGGGSGGGGS"   # (G4S)₂
```

### Knob-in-Hole Mutations (EU numbering)
```python
KNOB = {"T366W": "W"}
HOLE = {"T366S": "S", "L368A": "A", "Y407V": "V"}
```

## Source Antibody Sequences

### Teplizumab (hOKT3γ1 Ala-Ala)
- FDA approved 2022 for Type 1 diabetes delay
- Humanized IgG1 with Fc mutations (Ala-Ala) to reduce FcR binding
- CDRs derived from murine OKT3
- Does NOT cross-react with cynomolgus CD3

### SP34
- Murine antibody from 1985
- DOES cross-react with cynomolgus CD3 (important for preclinical)
- Needs humanization before use

### UCHT1
- Murine antibody from 1981
- Does NOT cross-react with cynomolgus
- Alternative epitope to OKT3

## Testing

```bash
pytest tests/ -v
pytest tests/test_liabilities.py -v  # Specific module
```

## Notes for Claude

1. **All tools must be permissively licensed** - never suggest IgFold or NetMHCIIpan
2. **Affinity target is ~50 nM** - not maximizing affinity, balancing with CRS risk
3. **Morrison format is symmetric** - both heavy chains have scFv/VHH fused to C-terminus
4. **CrossMab uses CH1-CL swap** - on the CD3 arm only, for correct light chain pairing
5. **Modal is required for GPU** - BoltzGen won't work locally without NVIDIA GPU
6. **Placeholder target is trastuzumab (HER2)** - configurable for actual use
7. **VH/VL pairs modeled as scFv** - structure prediction uses VH + linker + VL for paired chains
8. **Analysis modules use functions, not classes** - e.g., `score_humanness_pair()` not `HumannessScorer`
9. **Config supports nested schema** - `filtering.binding.min_pdockq` or flat `filtering.min_pdockq`
10. **scFv parsing available** - `parse_scfv()` and `is_likely_scfv()` in `src/utils/constants.py` extract VH/VL from concatenated scFv sequences
11. **CDR-specific liability filtering** - `CandidateScore` has `cdr_*_count` fields; `allow_deamidation_cdr=False` only rejects CDR liabilities
12. **Calibration uses full scFv** - For VH/VL pairs, calibration constructs scFv (not VH-only) for accurate threshold setting
13. **Contact counts are residue-level** - `num_contacts` counts unique residue pairs, not atomic contacts
14. **Epitope residues configurable** - `config.epitope.okt3_epitope_residues` overrides hardcoded defaults
15. **Aggregation filter active** - CDR-specific: >20% aromatic or 2+ consecutive aromatics. Fallback (no CDRs): >15% aromatic or 3+ consecutive aromatics.
16. **Fallback config wired** - `relax_soft_filters_first` and `max_threshold_relaxation` are now used
17. **CRITICAL: BoltzGen extracts target chain** - Multi-chain PDBs (1XIW has UCHT1 Fab, 1SY6 has OKT3 Fab) would bias designs. `boltzgen_app.py` extracts only the specified target chain before design.
18. **CRITICAL: 1XIW numbering starts at 12** - CD3ε chain A in 1XIW uses PDB residue numbers starting at 12, not 1. OKT3 epitope comparison handles this via sequence alignment. Never assume 1-indexed numbering for PDB data.
19. **1SY6 chain A is a fusion** - 1SY6 uses a CD3γ/ε fusion construct, NOT pure CD3ε. Use for epitope reference but prefer 1XIW for canonical CD3ε sequence.

## Dependencies

Core:
- Python 3.10+
- numpy, pandas, biopython, pyyaml
- anarci, abnumber (antibody numbering)
- biophi (humanness scoring)

Optional:
- torch, transformers, ablang2 (for local embeddings)
- modal (for GPU compute)
- pytest, black, ruff (development)
- jupyter, matplotlib, py3Dmol (notebooks)

## Ralph Review System

The `review/` directory contains a "Ralph Wiggum Loop" for automated scientific code review.

### How It Works
```
1. Spawn   → fresh Claude session (claude --print)
2. Inject  → prompt + ALL source code + iteration history
3. Attempt → Claude reviews, fixes issues, commits
4. Validate→ NO_ISSUES + git status clean required to exit
5. Log     → record findings to history file
6. Kill    → terminate session (context wiped)
7. Loop    → repeat until clean or max iterations
```

### Key Files
- `review/ralph_review.sh` - Main loop script
- `review/prompts/holistic.md` - Review criteria and instructions
- `review/logs/holistic_history.md` - Iteration history
- `review/tracking.yaml` - State tracking

### Usage
```bash
./review/ralph_review.sh              # Run review loop
./review/ralph_review.sh --status     # Check status
./review/ralph_review.sh --reset      # Start fresh
MAX_ITERATIONS=5 ./review/ralph_review.sh  # Limit iterations
```

### Known Issues / TODO
1. **Context budget**: Full source (~94k tokens) + system prompt (~23k) = ~117k tokens.
   Leaves limited room for response. May need to split repo or compress code.
2. **Token tracking**: JSON parsing for `--output-format json` has edge cases.
   The `extract_summary()` function can fail with multiline results.
3. **Coverage enforcement**: Prompt tells Claude to review all files but doesn't
   programmatically verify. Claude may skip files.
4. **Checklist validation**: Added requirement for review checklist before NO_ISSUES,
   but not yet tested thoroughly.

### Design Decisions
- **Full code in prompt** (not agent-style): Ensures Claude sees all code, can't skip files
- **Git-based validation**: Claude must commit fixes; uncommitted changes block NO_ISSUES
- **Fresh context each iteration**: Prevents confirmation bias from accumulated context
- **Structured history**: Logs what was reviewed/fixed so next iteration knows what's done

## References

Key papers for context:
- Staflin et al. (2021) - CD3 affinity tuning, Sci Rep
- Schaefer et al. (2011) - CrossMab format, PNAS
- BoltzGen (2024) - MIT Jameel Clinic
