# CD3 Binder Design Pipeline

A computational pipeline for designing CD3-binding domains (VHH nanobodies and scFvs) suitable for use in bispecific antibody therapeutics. This project uses open-source, permissively-licensed tools to generate, filter, and format CD3 binders with drug-like properties.

## Table of Contents

1. [Project Overview](#project-overview)
2. [Scientific Background](#scientific-background)
3. [Design Strategy](#design-strategy)
4. [Computational Pipeline](#computational-pipeline)
5. [Bispecific Antibody Formats](#bispecific-antibody-formats)
6. [Tool Stack & Licensing](#tool-stack--licensing)
7. [Repository Structure](#repository-structure)
8. [Implementation Details](#implementation-details)
9. [Expected Outputs](#expected-outputs)
10. [Installation & Usage](#installation--usage)
11. [Technical Considerations](#technical-considerations)
12. [References](#references)

---

## Project Overview

### Goal

Design ~10 CD3-binding protein sequences with favorable developability properties for experimental testing in cell lines and animal models. These binders will serve as the T-cell engaging arm in bispecific antibody constructs.

### Approach

1. **De novo design** using BoltzGen to generate novel VHH and scFv binders against CD3ε
2. **Optimization of existing binders** starting from known anti-CD3 antibodies (teplizumab/OKT3, SP34, UCHT1)
3. **Computational filtering** for developability, humanness, and sequence liabilities
4. **Format conversion** into multiple bispecific antibody architectures

### Key Constraints

- All tools must have **permissive licenses** (MIT, BSD, Apache, CC-BY) for commercial use
- Target **~50 nM Kd** for CD3 binding to balance efficacy with reduced cytokine release syndrome (CRS) risk
- Output sequences must be experimentally testable (expressible, stable, non-immunogenic)

---

## Scientific Background

### CD3 as a Bispecific Target

CD3 is a protein complex that associates with the T-cell receptor (TCR) and is essential for T-cell activation. In bispecific antibody therapeutics, one arm binds a tumor antigen while the other binds CD3, redirecting T-cells to kill tumor cells.

**Key considerations:**
- **Target subunit**: CD3ε (epsilon) is the standard target for T-cell engagers
- **Epitope**: The OKT3 epitope on CD3ε is well-validated clinically
- **Cross-reactivity**: SP34 antibody cross-reacts with cynomolgus CD3 (important for preclinical studies)

### CD3 Affinity and Cytokine Release Syndrome

Higher-affinity CD3 binders cause more T-cell activation and cytokine release:

| Affinity Class | Kd Range | CRS Risk | Clinical Examples |
|----------------|----------|----------|-------------------|
| High | < 1 nM | Higher | Early bispecifics, more CRS |
| Medium | 10-50 nM | Moderate | Modern designs |
| Low | 50-200+ nM | Lower | Intentionally attenuated |

**Literature evidence:**
- The threshold for T-cell cytotoxicity is more sensitive than the threshold for cytokine release
- Many modern bispecifics intentionally use ~50-100 nM CD3 arms to reduce CRS while maintaining efficacy
- Affinity attenuation of 10-100x from high-affinity binders can maintain killing while reducing cytokines

**Design target**: ~50 nM Kd as a starting point, with variants spanning 10-500 nM range.

### Source Structures

Two PDB structures will be used as design targets:

1. **1XIW**: CD3εδ heterodimer - provides the CD3ε surface for de novo design
2. **1SY6**: CD3εγ heterodimer in complex with OKT3 Fab - shows validated binding mode

### Starting Sequences for Optimization

All starting sequences are from expired patents or public domain publications:

| Antibody | Origin | Patent Status | Cross-Reactivity | Notes |
|----------|--------|---------------|------------------|-------|
| **Teplizumab** | Humanized OKT3 (hOKT3γ1 Ala-Ala) | Expired (1980s-90s) | Human only | FDA approved 2022 for T1D; Fc-silenced |
| **SP34** | Murine, 1985 | Expired | Human + Cynomolgus | Important for preclinical |
| **UCHT1** | Murine, 1981 | Expired | Human + Chimpanzee | No macaque cross-reactivity |

---

## Design Strategy

### Track 1: De Novo Design (BoltzGen)

BoltzGen is an open-source (MIT license) generative model for protein binder design. It can design:
- VHH nanobodies (~120 amino acids, single domain)
- scFvs (~250 amino acids, VH-linker-VL)

**Approach:**
1. Provide CD3ε structure as target (from 1XIW and 1SY6)
2. Generate 100-500 candidate binders per target structure
3. Filter computationally to top candidates
4. Predict structures and binding modes with Boltz-2

**Limitations:**
- BoltzGen cannot design full-length IgG directly
- Affinity is not explicitly controlled; filtering required
- Novel epitopes may differ from validated OKT3 site

### Track 2: Optimization of Existing Binders

Starting from known anti-CD3 CDR sequences, we will:
1. Graft CDRs onto human frameworks (for non-humanized sources)
2. Apply humanization via BioPhi/Sapiens
3. Generate affinity-attenuated variants using literature-known mutations
4. Screen for developability improvements

**Affinity attenuation strategy:**
- Include curated mutations from literature that reduce OKT3 affinity by ~10x, ~100x
- Generate combinatorial variants
- Predict effects on binding via structure modeling

### Track 3: Hybrid Approach

For scFv designs:
1. Use BoltzGen to design VH or VL individually
2. Pair with known complementary chain from existing binders
3. Optimize the combination

---

## Computational Pipeline

```
┌─────────────────────────────────────────────────────────────────────────┐
│                           INPUT PREPARATION                              │
├─────────────────────────────────────────────────────────────────────────┤
│  • Download CD3 structures (1XIW, 1SY6)                                 │
│  • Extract CD3ε chain sequences                                         │
│  • Prepare starting sequences (teplizumab, SP34, UCHT1)                 │
│  • Define affinity-attenuating mutations                                │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                         DESIGN GENERATION                                │
├──────────────────────────────┬──────────────────────────────────────────┤
│      De Novo (BoltzGen)      │         Optimization Track               │
│  ┌─────────────────────────┐ │ ┌──────────────────────────────────────┐ │
│  │ • VHH designs (100+)    │ │ │ • Humanize SP34, UCHT1 CDRs          │ │
│  │ • scFv designs (100+)   │ │ │ • Graft onto human frameworks        │ │
│  │ • Target: 1XIW, 1SY6    │ │ │ • Generate affinity variants         │ │
│  │ • Run on Modal (GPU)    │ │ │   (WT, 10x weaker, 100x weaker)      │ │
│  └─────────────────────────┘ │ │ • Combinatorial mutations            │ │
│                              │ └──────────────────────────────────────┘ │
└──────────────────────────────┴──────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     CDR IDENTIFICATION & NUMBERING                       │
├─────────────────────────────────────────────────────────────────────────┤
│  Tool: ANARCI (BSD license)                                             │
│  • Assign IMGT/Chothia/Kabat numbering                                  │
│  • Extract CDR-H1, H2, H3, L1, L2, L3 regions                          │
│  • Identify framework regions                                           │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       STRUCTURE PREDICTION                               │
├─────────────────────────────────────────────────────────────────────────┤
│  Step 1: Antibody Structure (ABodyBuilder2, BSD license)                │
│  • Predict VHH or Fv structure                                          │
│  • High accuracy for CDR loops                                          │
│                                                                         │
│  Step 2: Complex Structure (Boltz-2, MIT license)                       │
│  • Predict binder + CD3ε complex                                        │
│  • Assess binding mode and interface                                    │
│  • Extract pDockQ confidence score                                      │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        FILTERING CASCADE                                 │
├─────────────────────────────────────────────────────────────────────────┤
│  Filter 1: BINDING QUALITY                                              │
│  • Boltz-2 pDockQ > 0.5 (or pTM > 0.6)                                 │
│  • Interface buried surface area > 800 Å²                               │
│  • Contact residues with CD3ε > 10                                      │
│                                                                         │
│  Filter 2: HUMANNESS                                                    │
│  • BioPhi OASis score > 0.8 (for VH/VL)                                │
│  • Sapiens humanization suggestions applied                             │
│                                                                         │
│  Filter 3: SEQUENCE LIABILITIES                                         │
│  • Deamidation: No NG, NS motifs in CDRs                               │
│  • Isomerization: No DG, DS, DT motifs in CDRs                         │
│  • N-glycosylation: No N-X-S/T motifs in CDRs (where X ≠ P)            │
│  • Oxidation: Flag exposed Met/Trp in CDRs (soft filter)               │
│  • Unpaired cysteines: 0 (must be even count)                          │
│                                                                         │
│  Filter 4: DEVELOPABILITY PROPERTIES                                    │
│  • CDR-H3 length: 8-20 residues (typical range)                        │
│  • Net charge: -2 to +4 (avoid extremes)                               │
│  • Isoelectric point: 6.0-9.0                                          │
│  • Hydrophobic patches: ≤ 2 stretches of 5+ hydrophobic residues       │
│  • CDR hydrophobicity: Mean Kyte-Doolittle < 1.0                       │
│                                                                         │
│  Filter 5: AGGREGATION PROPENSITY                                       │
│  • Low aromatic content in CDRs (< 20%)                                │
│  • No consecutive aromatic residues (FF, WW, YY, FW, etc.)             │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       CANDIDATE RANKING                                  │
├─────────────────────────────────────────────────────────────────────────┤
│  Composite score based on:                                              │
│  • Binding confidence (pDockQ, 30% weight)                              │
│  • Humanness (OASis score, 25% weight)                                  │
│  • Liability count (inverse, 25% weight)                               │
│  • Developability metrics (20% weight)                                  │
│                                                                         │
│  Select top 15-20 candidates for format conversion                      │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                      FORMAT CONVERSION                                   │
├─────────────────────────────────────────────────────────────────────────┤
│  Convert each candidate into multiple bispecific formats:               │
│  • CrossMab (Fab × Fab)                                                 │
│  • Asymmetric Fab + scFv                                                │
│  • Asymmetric Fab + VHH                                                 │
│  • IgG-scFv Morrison (2x scFv, symmetric)                               │
│  • IgG-VHH Morrison (2x VHH, symmetric)                                 │
│                                                                         │
│  Include:                                                               │
│  • Human IgG1 Fc with knob-in-hole mutations (for asymmetric)          │
│  • Standard linkers (G4S)x3 for scFv, (G4S)x2 for Fc-VHH fusion        │
│  • Placeholder tumor target arm (anti-HER2 trastuzumab, configurable)  │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                         FINAL OUTPUT                                     │
├─────────────────────────────────────────────────────────────────────────┤
│  ~10 final candidates with:                                             │
│  • Protein sequences for all formats                                    │
│  • Predicted structures (PDB files)                                     │
│  • Developability scorecard (PDF/HTML report)                           │
│  • Ranked by composite score with rationale                             │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Bispecific Antibody Formats

### Overview

Five bispecific formats will be generated for each CD3 binder candidate. The tumor-targeting arm uses a placeholder (trastuzumab anti-HER2) that can be swapped for any target.

### Format 1: CrossMab (Fab × Fab)

```
         Fab (HER2)                 Fab (CD3)
              │                        │
              ▼                        ▼
         ┌────────┐              ┌────────┐
         │   VH   │              │   VH   │
         │   CH1  │              │   CH1* │ ← CrossMab swap
         ├────────┤              ├────────┤
         │   VL   │              │   VL   │
         │   CL   │              │   CL*  │ ← CrossMab swap
         └────┬───┘              └────┬───┘
              │     Knob-in-Hole      │
              │    (T366W / T366S,    │
              │     L368A, Y407V)     │
              └──────────┬────────────┘
                         │
                    ┌────┴────┐
                    │   Fc    │
                    │  (IgG1) │
                    └─────────┘
```

**Components:**
- Heavy chain 1: VH(HER2) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2: VH(CD3) - CL - Hinge - CH2 - CH3(Hole) [CrossMab swap]
- Light chain 1: VL(HER2) - CL
- Light chain 2: VL(CD3) - CH1 [CrossMab swap]

**Knob-in-hole mutations (EU numbering):**
- Knob: T366W
- Hole: T366S, L368A, Y407V

**CrossMab CH1-CL swap:**
On the CD3 arm, CH1 and CL domains are swapped between heavy and light chains to ensure correct light chain pairing.

### Format 2: Asymmetric Knob-in-Hole (Fab + scFv)

```
         Fab (HER2)                 scFv (CD3)
              │                        │
              ▼                        ▼
         ┌────────┐              ┌──────────────┐
         │   VH   │              │ VH-(G4S)3-VL │
         │   CH1  │              └──────┬───────┘
         ├────────┤                     │
         │   VL   │                     │
         │   CL   │                     │
         └────┬───┘                     │
              │     Knob-in-Hole        │
              └──────────┬──────────────┘
                         │
                    ┌────┴────┐
                    │   Fc    │
                    └─────────┘
```

**Components:**
- Heavy chain 1: VH(HER2) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2: scFv(CD3) - Hinge - CH2 - CH3(Hole)
- Light chain: VL(HER2) - CL (only one light chain needed)

**scFv linker:** (GGGGS)₃ = GGGGSGGGGSGGGGS (15 aa)

### Format 3: Asymmetric Knob-in-Hole (Fab + VHH)

```
         Fab (HER2)                 VHH (CD3)
              │                        │
              ▼                        ▼
         ┌────────┐              ┌──────────┐
         │   VH   │              │   VHH    │
         │   CH1  │              └────┬─────┘
         ├────────┤                   │
         │   VL   │                   │
         │   CL   │                   │
         └────┬───┘                   │
              │     Knob-in-Hole      │
              └──────────┬────────────┘
                         │
                    ┌────┴────┐
                    │   Fc    │
                    └─────────┘
```

**Components:**
- Heavy chain 1: VH(HER2) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2: VHH(CD3) - Hinge - CH2 - CH3(Hole)
- Light chain: VL(HER2) - CL

### Format 4: IgG-scFv Morrison (Symmetric, 2x CD3)

```
         Fab (HER2)              Fab (HER2)
              │                        │
              ▼                        ▼
         ┌────────┐              ┌────────┐
         │   VH   │              │   VH   │
         │   CH1  │              │   CH1  │
         ├────────┤              ├────────┤
         │   VL   │              │   VL   │
         │   CL   │              │   CL   │
         └────┬───┘              └────┬───┘
              │                        │
              └──────────┬─────────────┘
                         │
                    ┌────┴────┐
                    │   Fc    │
                    └────┬────┘
                    ╱         ╲
              ┌────┴──┐    ┌──┴────┐
              │ scFv  │    │ scFv  │  ← 2x CD3 (symmetric)
              │ (CD3) │    │ (CD3) │
              └───────┘    └───────┘
```

**Components:**
- Heavy chain: VH(HER2) - CH1 - Hinge - CH2 - CH3 - Linker - scFv(CD3)
- Light chain: VL(HER2) - CL

**C-terminal fusion linker:** (GGGGS)₂ = GGGGSGGGGS (10 aa)

**Valency:**
- 2x HER2 binding (bivalent)
- 2x CD3 binding (bivalent)

### Format 5: IgG-VHH Morrison (Symmetric, 2x CD3)

```
         Fab (HER2)              Fab (HER2)
              │                        │
              ▼                        ▼
         ┌────────┐              ┌────────┐
         │   VH   │              │   VH   │
         │   CH1  │              │   CH1  │
         ├────────┤              ├────────┤
         │   VL   │              │   VL   │
         │   CL   │              │   CL   │
         └────┬───┘              └────┬───┘
              │                        │
              └──────────┬─────────────┘
                         │
                    ┌────┴────┐
                    │   Fc    │
                    └────┬────┘
                    ╱         ╲
              ┌────┴──┐    ┌──┴────┐
              │  VHH  │    │  VHH  │  ← 2x CD3 (symmetric)
              │ (CD3) │    │ (CD3) │
              └───────┘    └───────┘
```

**Components:**
- Heavy chain: VH(HER2) - CH1 - Hinge - CH2 - CH3 - Linker - VHH(CD3)
- Light chain: VL(HER2) - CL

---

## Tool Stack & Licensing

All computational tools have permissive licenses suitable for commercial use:

| Tool | Purpose | License | Commercial OK |
|------|---------|---------|---------------|
| **BoltzGen** | De novo binder design | MIT | Yes |
| **Boltz-2** | Complex structure prediction | MIT | Yes |
| **ABodyBuilder2** | Antibody structure prediction | BSD 3-Clause | Yes |
| **ANARCI** | Antibody numbering | BSD 3-Clause | Yes |
| **BioPhi/Sapiens** | Humanness scoring | MIT | Yes |
| **AbLang** | Sequence embeddings | CC BY 4.0 | Yes (with attribution) |

**Excluded tools (non-permissive licenses):**
- IgFold (JHU Academic License - non-commercial only)
- NetMHCIIpan (DTU academic license)
- CamSol (web server only, no standalone code)

### Compute Infrastructure

**Modal** (https://modal.com) is used for GPU-intensive operations:
- BoltzGen design generation
- Boltz-2 complex prediction
- ABodyBuilder2 structure prediction (optional, can run locally)

Modal provides on-demand A100/H100 GPUs with simple Python deployment. This is preferred over local execution because:
1. BoltzGen requires CUDA GPUs (not compatible with Apple MPS)
2. No local GPU setup required
3. Pay-per-use pricing (~$0.001-0.01 per design)
4. Scales to parallel batch processing

---

## Repository Structure

```
cd3-binder-design/
├── README.md                            # This document
├── CLAUDE.md                            # AI assistant context
├── pyproject.toml                       # Python package configuration
├── requirements.txt                     # Dependencies
│
├── data/
│   ├── targets/
│   │   ├── cd3_epsilon_delta_1XIW.pdb   # CD3εδ structure for design
│   │   ├── cd3_epsilon_gamma_1SY6.pdb   # CD3εγ + OKT3 complex
│   │   ├── cd3_epsilon_chain.fasta      # Extracted CD3ε sequence
│   │   └── README.md                    # Structure preparation notes
│   │
│   ├── starting_sequences/
│   │   ├── teplizumab.yaml              # hOKT3γ1(Ala-Ala) VH/VL
│   │   ├── sp34.yaml                    # SP34 VH/VL (murine + humanized)
│   │   ├── ucht1.yaml                   # UCHT1 VH/VL (murine + humanized)
│   │   ├── affinity_mutations.yaml      # Literature affinity-reducing mutations
│   │   └── README.md                    # Sequence sources and references
│   │
│   ├── frameworks/
│   │   ├── human_vh_frameworks.fasta    # IGHV germline frameworks
│   │   ├── human_vl_kappa.fasta         # IGKV germline frameworks
│   │   ├── human_vl_lambda.fasta        # IGLV germline frameworks
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
│   │   ├── developability.py            # Combined developability scoring
│   │   └── embeddings.py                # AbLang embeddings (optional)
│   │
│   ├── formatting/
│   │   ├── __init__.py
│   │   ├── base.py                      # Base bispecific formatter
│   │   ├── crossmab.py                  # CrossMab Fab×Fab assembly
│   │   ├── fab_scfv.py                  # Asymmetric Fab + scFv
│   │   ├── fab_vhh.py                   # Asymmetric Fab + VHH
│   │   ├── igg_scfv.py                  # Morrison IgG-(scFv)₂
│   │   ├── igg_vhh.py                   # Morrison IgG-(VHH)₂
│   │   └── sequence_utils.py            # Linker insertion, grafting
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
│       ├── io.py                        # File I/O utilities
│       ├── logging.py                   # Logging configuration
│       └── constants.py                 # Amino acid properties, etc.
│
├── modal/
│   ├── __init__.py
│   ├── boltzgen_app.py                  # BoltzGen Modal deployment
│   ├── boltz2_app.py                    # Boltz-2 Modal deployment
│   └── abodybuilder_app.py              # ABodyBuilder2 Modal deployment
│
├── scripts/
│   ├── 01_setup_targets.py              # Download and prepare structures
│   ├── 02_run_denovo_design.py          # Run BoltzGen on Modal
│   ├── 03_run_optimization.py           # Generate optimized variants
│   ├── 04_predict_structures.py         # Run structure prediction
│   ├── 05_filter_candidates.py          # Apply filtering cascade
│   ├── 06_format_bispecifics.py         # Convert to bispecific formats
│   ├── 07_generate_report.py            # Generate final report
│   └── run_full_pipeline.py             # Execute all steps
│
├── notebooks/
│   ├── 01_explore_cd3_structures.ipynb  # Analyze target structures
│   ├── 02_analyze_designs.ipynb         # Review generated designs
│   ├── 03_compare_candidates.ipynb      # Compare final candidates
│   └── 04_visualize_complexes.ipynb     # 3D visualization
│
└── tests/
    ├── __init__.py
    ├── test_design.py
    ├── test_analysis.py
    ├── test_formatting.py
    └── fixtures/                        # Test data
```

---

## Implementation Details

### Configuration

Pipeline behavior is controlled by a YAML configuration file:

```yaml
# config.yaml
design:
  denovo:
    num_vhh_designs: 200          # VHH candidates from BoltzGen
    num_scfv_designs: 200         # scFv candidates from BoltzGen
    target_structures:
      - data/targets/cd3_epsilon_delta_1XIW.pdb
      - data/targets/cd3_epsilon_gamma_1SY6.pdb

  optimization:
    starting_sequences:
      - teplizumab                # Primary optimization target
      - sp34                      # Cynomolgus cross-reactive
      - ucht1                     # Alternative epitope
    affinity_variants:
      - wild_type
      - 10x_weaker
      - 100x_weaker

filtering:
  binding:
    min_pdockq: 0.5
    min_interface_area: 800       # Å²
    min_contacts: 10

  humanness:
    min_oasis_score: 0.8

  liabilities:
    allow_deamidation_cdr: false
    allow_isomerization_cdr: false
    allow_glycosylation_cdr: false
    max_oxidation_sites: 2        # Soft filter

  developability:
    cdr_h3_length_range: [8, 20]
    net_charge_range: [-2, 4]
    pi_range: [6.0, 9.0]
    max_hydrophobic_patches: 2

formatting:
  tumor_target: trastuzumab       # Placeholder, configurable
  formats:
    - crossmab
    - fab_scfv
    - fab_vhh
    - igg_scfv
    - igg_vhh

  linkers:
    scfv: "GGGGSGGGGSGGGGS"       # (G4S)₃
    fc_fusion: "GGGGSGGGGS"       # (G4S)₂

output:
  num_final_candidates: 10
  include_structures: true
  generate_report: true
```

### Sequence Liability Detection

```python
# src/analysis/liabilities.py

DEAMIDATION_MOTIFS = ['NG', 'NS', 'NT', 'ND']  # N followed by small residue
ISOMERIZATION_MOTIFS = ['DG', 'DS', 'DT', 'DD']  # D followed by small residue
OXIDATION_RESIDUES = ['M', 'W']  # Methionine, Tryptophan

def find_glycosylation_sites(sequence: str) -> list[int]:
    """Find N-X-S/T motifs where X != P."""
    sites = []
    for i in range(len(sequence) - 2):
        if (sequence[i] == 'N' and
            sequence[i+1] != 'P' and
            sequence[i+2] in ['S', 'T']):
            sites.append(i)
    return sites

def count_unpaired_cysteines(sequence: str) -> int:
    """Odd cysteine count suggests unpaired."""
    return sequence.count('C') % 2
```

### Humanness Scoring

```python
# src/analysis/humanness.py
from biophi import Sapiens, OASis

def score_humanness(vh_sequence: str, vl_sequence: str = None) -> dict:
    """Score humanness using BioPhi OASis."""
    oasis = OASis()

    vh_score = oasis.score(vh_sequence, chain_type='H')

    result = {'vh_oasis': vh_score}

    if vl_sequence:
        vl_score = oasis.score(vl_sequence, chain_type='L')
        result['vl_oasis'] = vl_score
        result['mean_oasis'] = (vh_score + vl_score) / 2

    return result
```

### Bispecific Assembly

```python
# src/formatting/crossmab.py

def assemble_crossmab(
    target_vh: str, target_vl: str,
    cd3_vh: str, cd3_vl: str,
    fc_knob: str, fc_hole: str,
    ch1: str, cl: str, hinge: str
) -> dict:
    """Assemble CrossMab bispecific sequences."""

    # Target arm: standard Fab
    heavy_chain_1 = target_vh + ch1 + hinge + fc_knob
    light_chain_1 = target_vl + cl

    # CD3 arm: CrossMab swap (CH1-CL exchanged)
    heavy_chain_2 = cd3_vh + cl + hinge + fc_hole  # VH fused to CL
    light_chain_2 = cd3_vl + ch1                    # VL fused to CH1

    return {
        'heavy_chain_1': heavy_chain_1,  # Target arm (knob)
        'heavy_chain_2': heavy_chain_2,  # CD3 arm (hole, CrossMab)
        'light_chain_1': light_chain_1,  # Target light chain
        'light_chain_2': light_chain_2,  # CD3 light chain (CrossMab)
    }
```

### Modal Deployment for BoltzGen

```python
# modal/boltzgen_app.py
import modal

app = modal.App("boltzgen-cd3")

image = modal.Image.debian_slim().pip_install(
    "boltzgen",
    "torch",
    "biopython",
)

@app.function(
    image=image,
    gpu="A100",
    timeout=3600,
)
def run_boltzgen(
    target_pdb_path: str,
    output_type: str = "vhh",  # or "scfv"
    num_designs: int = 100,
    seed: int = 42,
) -> list[dict]:
    """Run BoltzGen to design binders."""
    from boltzgen import BoltzGen

    model = BoltzGen.load()

    designs = model.design(
        target_structure=target_pdb_path,
        binder_type=output_type,
        num_samples=num_designs,
        seed=seed,
    )

    return [
        {"sequence": d.sequence, "confidence": d.confidence}
        for d in designs
    ]
```

---

## Expected Outputs

### Per-Candidate Output

For each of the ~10 final candidates:

1. **Sequences** (FASTA/YAML)
   - VHH or scFv sequence (CD3 binder)
   - Full bispecific sequences for all 5 formats

2. **Structures** (PDB files)
   - Predicted VHH/Fv structure (ABodyBuilder2)
   - Predicted complex with CD3ε (Boltz-2)

3. **Developability Scorecard** (JSON/PDF)
   ```json
   {
     "candidate_id": "denovo_vhh_042",
     "source": "boltzgen_1XIW",
     "binding": {
       "pdockq": 0.72,
       "interface_area": 1245.3,
       "contacts": 18
     },
     "humanness": {
       "oasis_score": 0.87
     },
     "liabilities": {
       "deamidation_sites": [],
       "glycosylation_sites": [],
       "isomerization_sites": [],
       "oxidation_sites": [{"position": 55, "residue": "M"}],
       "unpaired_cys": 0
     },
     "developability": {
       "cdr_h3_length": 14,
       "net_charge": 1,
       "isoelectric_point": 7.8,
       "hydrophobic_patches": 1
     },
     "composite_score": 0.82,
     "rank": 3
   }
   ```

### Summary Report

HTML/PDF report containing:
- Pipeline configuration and parameters
- Summary statistics (designs generated, filtered, selected)
- Ranked candidate table with key metrics
- Structure visualizations
- Sequence alignments
- Recommendations for experimental testing

---

## Installation & Usage

### Prerequisites

- Python 3.10+
- Modal account (for GPU compute): https://modal.com
- ~10GB disk space for structures and models

### Installation

```bash
# Clone repository
git clone https://github.com/alecnielsen/cd3-binder-design.git
cd cd3-binder-design

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -e .

# Install Modal CLI
pip install modal
modal setup  # Follow prompts to authenticate

# Deploy Modal functions
modal deploy modal/boltzgen_app.py
modal deploy modal/boltz2_app.py
```

### Running the Pipeline

```bash
# Step-by-step execution
python scripts/01_setup_targets.py      # Download CD3 structures
python scripts/02_run_denovo_design.py  # Run BoltzGen (Modal)
python scripts/03_run_optimization.py   # Generate optimized variants
python scripts/04_predict_structures.py # Predict structures (Modal)
python scripts/05_filter_candidates.py  # Apply filters
python scripts/06_format_bispecifics.py # Generate bispecific sequences
python scripts/07_generate_report.py    # Create final report

# Or run everything
python scripts/run_full_pipeline.py --config config.yaml
```

### Configuration

Edit `config.yaml` to customize:
- Number of designs to generate
- Filtering thresholds
- Bispecific formats to include
- Tumor target arm (default: trastuzumab/HER2)

---

## Technical Considerations

### Limitations

1. **Affinity prediction**: Computational tools predict binding mode but not precise Kd. Experimental validation required.

2. **BoltzGen novelty**: De novo designs may bind different epitopes than OKT3. This could be advantageous (novel IP) or problematic (different biology).

3. **Cross-reactivity**: BoltzGen designs may not cross-react with cynomolgus CD3. SP34-derived sequences recommended for preclinical studies.

4. **Developability predictions**: Computational filters reduce risk but don't guarantee expressibility or stability. Experimental screening required.

5. **Immunogenicity**: T-cell epitope prediction (NetMHCIIpan) excluded due to license restrictions. Humanness scoring (BioPhi) provides partial mitigation.

### Assumptions

1. CD3ε is the binding target (most common for bispecifics)
2. ~50 nM Kd is acceptable for efficacy with reduced CRS
3. All starting sequences are off-patent and freely usable
4. Placeholder target arm (HER2) will be replaced for actual use

### Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| BoltzGen fails to generate binders | Optimization track provides fallback |
| All designs fail filters | Relax thresholds iteratively |
| Poor binding predictions | Include known binders as positive controls |
| Aggregation issues | Include aggregation-prone sequence filters |

---

## References

### CD3 Biology and Bispecifics

1. Wolf E, et al. (2005). BiTEs: bispecific antibody constructs with unique anti-tumor activity. Drug Discovery Today.
2. Goebeler ME, Bargou RC (2020). T cell-engaging therapies — BiTEs and beyond. Nature Reviews Clinical Oncology.

### CD3 Affinity and Safety

3. Staflin K, et al. (2021). Generation of T-cell-redirecting bispecific antibodies with differentiated profiles of cytokine release and biodistribution by CD3 affinity tuning. Scientific Reports. https://www.nature.com/articles/s41598-021-93842-0
4. Hernandez-Hoyos G, et al. (2016). MOR209/ES414, a novel bispecific antibody targeting PSMA for the treatment of metastatic castration-resistant prostate cancer. Molecular Cancer Therapeutics.

### Source Antibodies

5. Kung P, et al. (1979). Monoclonal antibodies defining distinctive human T cell surface antigens. Science. (OKT3)
6. Pessano S, et al. (1985). The T3/T cell receptor complex: antigenic distinction between the two 20-kd T3 (T3-delta and T3-epsilon) subunits. EMBO J. (SP34)
7. Callard RE, et al. (1981). A new approach to the generation of cytotoxic T cells. Clin Exp Immunol. (UCHT1)

### Computational Tools

8. Abanades B, et al. (2023). ImmuneBuilder: Deep-Learning models for predicting the structures of immune proteins. Communications Biology. https://github.com/oxpig/ImmuneBuilder
9. Dunbar J, Deane CM (2016). ANARCI: antigen receptor numbering and receptor classification. Bioinformatics. https://github.com/oxpig/ANARCI
10. Prihoda D, et al. (2022). BioPhi: A platform for antibody design, humanization, and humanness evaluation. mAbs. https://github.com/Merck/BioPhi
11. BoltzGen (2024). MIT Jameel Clinic. https://jclinic.mit.edu/boltzgen/

### Bispecific Formats

12. Brinkmann U, Kontermann RE (2017). The making of bispecific antibodies. mAbs.
13. Schaefer W, et al. (2011). Immunoglobulin domain crossover as a generic approach for the production of bispecific IgG antibodies. PNAS. (CrossMab)

---

## License

MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
