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
results = scanner.scan_sequence("EVQLVESGGGLVQ...", region="CDR")
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
DEAMIDATION = ['NG', 'NS', 'NT', 'ND']
ISOMERIZATION = ['DG', 'DS', 'DT', 'DD']
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

## References

Key papers for context:
- Staflin et al. (2021) - CD3 affinity tuning, Sci Rep
- Schaefer et al. (2011) - CrossMab format, PNAS
- BoltzGen (2024) - MIT Jameel Clinic
