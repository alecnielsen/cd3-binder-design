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
11. [Known Limitations](#known-limitations)
12. [Experimental Validation Requirements](#experimental-validation-requirements)
13. [Reproducibility](#reproducibility)
14. [Technical Considerations](#technical-considerations)
15. [References](#references)

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
- **Design hypothesis**: ~50 nM Kd may balance efficacy with reduced CRS (see [Known Limitations](#known-limitations) for why this cannot be computationally enforced)
- Output sequences must be experimentally testable (expressible, stable, non-immunogenic)
- **Experimental validation is required** - computational filtering reduces risk but does not guarantee success

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

**Design hypothesis**: ~50 nM Kd may be optimal, with variants spanning 10-500 nM range. However, computational tools cannot predict absolute Kd values (see [Known Limitations](#known-limitations)). This hypothesis must be validated experimentally across an affinity panel.

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
│                    CALIBRATION PHASE (Run Once)                          │
├─────────────────────────────────────────────────────────────────────────┤
│  Purpose: Establish filter thresholds using KNOWN BINDERS               │
│                                                                         │
│  Step 1: Run Boltz-2 on teplizumab/SP34/UCHT1 + CD3ε complexes         │
│  Step 2: Record pDockQ, interface area, contact counts for each         │
│  Step 3: Set thresholds to ACCEPT all known binders:                    │
│          • min_pdockq = min(known_binder_pdockq) - 0.05                │
│          • min_interface_area = min(known_binder_area) - 100           │
│          • min_contacts = min(known_binder_contacts) - 2               │
│  Step 4: Document calibration results in config.yaml                    │
│                                                                         │
│  CRITICAL: If known binders fail default thresholds, thresholds are     │
│  wrong, not the binders. Calibration prevents false negative filtering. │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        FILTERING CASCADE                                 │
├─────────────────────────────────────────────────────────────────────────┤
│  NOTE: All thresholds below are DEFAULTS. Use calibrated thresholds     │
│  from the calibration phase for actual filtering runs.                  │
│                                                                         │
│  Filter 1: BINDING QUALITY (calibrated thresholds)                      │
│  • Boltz-2 pDockQ > [calibrated, default 0.5]                          │
│  • Interface buried surface area > [calibrated, default 800 Å²]        │
│  • Contact residues with CD3ε > [calibrated, default 10]               │
│  NOTE: pDockQ is a structural confidence score, NOT an affinity         │
│  predictor. High pDockQ ≠ high affinity. See Known Limitations.         │
│                                                                         │
│  Filter 2: HUMANNESS                                                    │
│  • BioPhi OASis score > 0.8 (for VH/VL)                                │
│  • Sapiens humanization suggestions applied                             │
│  • If humanization significantly changes CDRs, FLAG for back-mutation   │
│    testing (humanized + key CDR residues restored)                      │
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
│                       EPITOPE ANNOTATION                                 │
├─────────────────────────────────────────────────────────────────────────┤
│  For de novo designs, annotate predicted epitope vs OKT3 epitope:       │
│  • Extract CD3ε contact residues from Boltz-2 complex                   │
│  • Compare to OKT3 epitope residues (from 1SY6 structure)              │
│  • Calculate epitope overlap percentage                                 │
│  • Flag designs with < 50% overlap as "NOVEL EPITOPE"                  │
│                                                                         │
│  NOTE: Novel epitopes may have different biology (efficacy, safety).    │
│  OKT3-like epitopes have clinical validation. Both are included in      │
│  final candidates, but epitope type is tracked for prioritization.      │
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
│  DIVERSITY REQUIREMENT:                                                 │
│  • Cluster candidates by CDR-H3 sequence similarity (70% identity)     │
│  • Select top candidate from each cluster                               │
│  • Ensure mix of: de novo + optimized, VHH + scFv, OKT3 + novel epitope │
│                                                                         │
│  INSUFFICIENT CANDIDATES FALLBACK:                                      │
│  If < 10 candidates survive filtering:                                  │
│  1. Relax soft filters (oxidation, hydrophobic patches)                │
│  2. Relax thresholds by 10% toward calibration baseline                │
│  3. Include "borderline" candidates with explicit risk flags           │
│  4. Document relaxations in output report                              │
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
│   ├── 00_run_calibration.py            # Calibrate thresholds with known binders
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

Pipeline behavior is controlled by `config.yaml` at the repository root. Copy and modify this file to customize pipeline behavior:

```yaml
# config.yaml (included in repo)
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

calibration:
  # Known binders for threshold calibration
  positive_controls:
    - teplizumab
    - sp34
    - ucht1
  # Calibrated thresholds (set after calibration run)
  calibrated_thresholds:
    min_pdockq: null              # Set by calibration
    min_interface_area: null      # Set by calibration
    min_contacts: null            # Set by calibration
  calibration_margin:
    pdockq: 0.05                  # Subtract from min known binder
    interface_area: 100           # Subtract from min known binder
    contacts: 2                   # Subtract from min known binder

filtering:
  binding:
    # DEFAULTS - use calibrated values if available
    min_pdockq: 0.5               # Default, override with calibrated
    min_interface_area: 800       # Å², override with calibrated
    min_contacts: 10              # Override with calibrated
    use_calibrated: true          # Use calibrated thresholds if available

  humanness:
    min_oasis_score: 0.8
    generate_back_mutations: true # Generate CDR back-mutation variants

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

  fallback:
    min_candidates: 10            # Trigger fallback if fewer survive
    relax_soft_filters_first: true
    max_threshold_relaxation: 0.1 # 10% toward calibration baseline

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
  include_provenance: true        # Add metadata to all outputs
  output_dir: data/outputs

reproducibility:
  boltzgen_seed: 42
  sampling_seed: 12345
  clustering_seed: 0

epitope_annotation:
  okt3_epitope_residues: [23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 38, 39, 40, 41, 42, 45, 47]  # CD3ε residues from 1SY6
  overlap_threshold: 0.5          # Flag as "novel epitope" if < 50% overlap
```

### Sequence Liability Detection

```python
# src/analysis/liabilities.py

DEAMIDATION_MOTIFS = ['NG', 'NS', 'NT', 'ND', 'NH']  # N followed by small/flexible residue
ISOMERIZATION_MOTIFS = ['DG', 'DS', 'DT', 'DD', 'DH', 'DN']  # D followed by small/flexible residue
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
       "pdockq_note": "Structural confidence, NOT affinity predictor",
       "interface_area": 1245.3,
       "contacts": 18,
       "passed_calibrated_threshold": true
     },
     "epitope": {
       "cd3_contact_residues": [23, 25, 26, 28, 29, 31, 45, 47],
       "okt3_overlap_percent": 0.75,
       "epitope_class": "OKT3-like"
     },
     "humanness": {
       "oasis_score": 0.87,
       "sapiens_mutations_applied": 3,
       "back_mutation_variant_generated": true
     },
     "liabilities": {
       "deamidation_sites": [],
       "glycosylation_sites": [],
       "isomerization_sites": [],
       "oxidation_sites": [{"position": 55, "residue": "M", "severity": "soft"}],
       "unpaired_cys": 0
     },
     "developability": {
       "cdr_h3_length": 14,
       "net_charge": 1,
       "isoelectric_point": 7.8,
       "hydrophobic_patches": 1
     },
     "risk_flags": [],
     "composite_score": 0.82,
     "rank": 3,
     "_provenance": {
       "pipeline_version": "1.0.0",
       "git_commit": "abc1234",
       "run_timestamp": "2026-01-15T10:30:00Z"
     }
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
# Step 0: Calibration (RUN FIRST - sets filter thresholds)
python scripts/00_run_calibration.py    # Calibrate with known binders

# Step-by-step execution
python scripts/01_setup_targets.py      # Download CD3 structures
python scripts/02_run_denovo_design.py  # Run BoltzGen (Modal)
python scripts/03_run_optimization.py   # Generate optimized variants
python scripts/04_predict_structures.py # Predict structures (Modal)
python scripts/05_filter_candidates.py  # Apply filters (uses calibrated thresholds)
python scripts/06_format_bispecifics.py # Generate bispecific sequences
python scripts/07_generate_report.py    # Create final report

# Or run everything (includes calibration)
python scripts/run_full_pipeline.py --config config.yaml
```

**Important**: Always run calibration before filtering. The calibration step uses known binders (teplizumab, SP34, UCHT1) to set appropriate filter thresholds. Skipping calibration may result in incorrect filtering.

### Configuration

Edit `config.yaml` to customize:
- Number of designs to generate
- Filtering thresholds
- Bispecific formats to include
- Tumor target arm (default: trastuzumab/HER2)

---

## Critical Assumptions for Scientific Review

This section documents the key assumptions underlying this pipeline. **Reviewers should challenge these assumptions** and identify gaps in the reasoning.

### Assumption 1: BoltzGen Can Design Functional CD3 Binders

**Claim**: BoltzGen can generate VHH/scFv sequences that bind CD3ε with reasonable affinity.

**Evidence supporting this**:
- BoltzGen has been validated on other protein targets in published benchmarks
- The model was trained on protein-protein interface data that includes antibody-antigen complexes

**Evidence against / Unknowns**:
- BoltzGen has not been specifically validated on CD3 or T-cell surface proteins
- CD3ε forms a complex (with δ or γ); the binding site may require both subunits
- No published data on BoltzGen success rates for immune receptor targets
- De novo designs may have very low hit rates (<1% functional binders)

**What would falsify this assumption**:
- Zero designs pass calibration thresholds after filtering
- All designs fail to bind in SPR/ELISA
- Designs bind but don't activate T-cells

**Fallback if assumption fails**: Rely entirely on optimization track (teplizumab, SP34, UCHT1 variants).

### Assumption 2: Boltz-2 Complex Predictions Are Meaningful

**Claim**: Boltz-2 can predict binder-CD3ε complex structures accurately enough for filtering.

**Evidence supporting this**:
- Boltz-2 achieves high accuracy on protein complex benchmarks
- pDockQ correlates with DockQ (actual structural accuracy)

**Evidence against / Unknowns**:
- Antibody-antigen complexes may be harder than general protein complexes
- Boltz-2 was not specifically benchmarked on VHH/scFv-antigen complexes
- The model may hallucinate interfaces that look plausible but don't exist

**What would falsify this assumption**:
- Known binders (teplizumab, SP34) score poorly or fail calibration
- Interface residues predicted by Boltz-2 don't match known OKT3 epitope
- No correlation between pDockQ and experimental binding

**How to detect early**: Calibration phase should catch this. If known binders fail, the model is unreliable.

### Assumption 3: ~50 nM Kd Balances Efficacy and Safety

**Claim**: CD3 binders with ~50 nM Kd provide sufficient T-cell killing with reduced cytokine release.

**Evidence supporting this**:
- Staflin et al. (2021) showed affinity-attenuated CD3 arms reduce cytokines while maintaining killing
- Multiple clinical programs have shifted to lower-affinity CD3 arms
- The cytokine threshold appears higher than the cytotoxicity threshold

**Evidence against / Unknowns**:
- Optimal Kd may be context-dependent (tumor type, tumor antigen density)
- 50 nM may be too weak for some indications
- The therapeutic window varies by patient (genetics, tumor burden, prior treatment)

**What would falsify this assumption**:
- All ~50 nM binders fail to kill tumor cells in co-culture assays
- Cytokine release is not reduced compared to high-affinity binders

**Key question**: This pipeline cannot enforce 50 nM Kd computationally. Is the affinity panel approach sufficient?

### Assumption 4: Humanization Preserves Function

**Claim**: BioPhi/Sapiens humanization produces sequences that retain binding activity.

**Evidence supporting this**:
- BioPhi was validated on therapeutic antibodies
- Framework humanization is generally safer than CDR modifications
- Back-mutation variants provide a safety net

**Evidence against / Unknowns**:
- Humanization success rate is typically 30-70%, not 100%
- Some CDRs are sensitive to even conservative mutations
- VHH humanization is less validated than conventional antibody humanization

**What would falsify this assumption**:
- All humanized variants lose binding
- Back-mutation variants have poor humanness scores

### Assumption 5: Developability Filters Are Predictive

**Claim**: Sequences passing developability filters (liabilities, charge, hydrophobicity) are more likely to succeed experimentally.

**Evidence supporting this**:
- Literature supports correlation between sequence features and developability
- Deamidation, glycosylation, and aggregation motifs are well-characterized risks

**Evidence against / Unknowns**:
- Filters are probabilistic, not deterministic
- Some "problematic" sequences express well; some "clean" sequences fail
- CDR context matters (buried vs. exposed positions)

**What would falsify this assumption**:
- No correlation between filter scores and experimental success
- High-scoring candidates fail at similar rates to low-scoring candidates

---

## Open Questions for Reviewers

The following questions represent gaps in the pipeline design that reviewers should consider:

### Biological Questions

1. **Is CD3ε sufficient as a binding target?**
   - CD3 exists as CD3εδ and CD3εγ heterodimers. Should the pipeline target specific heterodimers?
   - Does binding to isolated CD3ε predict binding to CD3 in the context of the full TCR complex?

2. **Is the OKT3 epitope optimal?**
   - OKT3/teplizumab binds a specific epitope on CD3ε. Is this the best epitope for T-cell engagers?
   - Could alternative epitopes have better safety profiles or efficacy?

3. **What about CD3ε conformational states?**
   - CD3ε may adopt different conformations when associated with TCR vs. free
   - The static structures (1XIW, 1SY6) may not capture relevant conformational diversity

4. **Cynomolgus cross-reactivity**
   - Only SP34 cross-reacts with cynomolgus CD3
   - Should cyno cross-reactivity be a hard requirement for preclinical candidates?

### Computational Questions

5. **Are 200 BoltzGen designs sufficient?**
   - If hit rate is 1%, 200 designs yields ~2 hits
   - Should the design count be higher (500-1000)?

6. **Is pDockQ the right filtering metric?**
   - pDockQ measures structural confidence, not binding affinity
   - Should alternative metrics (interface energy, shape complementarity) be included?

7. **How should VHH vs scFv designs be prioritized?**
   - VHHs are smaller and may have better tissue penetration
   - scFvs have two binding loops and may have higher affinity potential
   - The pipeline generates both but doesn't provide guidance on prioritization

8. **What if calibration fails?**
   - If known binders score poorly in Boltz-2, what are the options?
   - Should the pipeline fall back to sequence similarity instead of structural prediction?

### Practical Questions

9. **What is the expected computational cost?**
   - BoltzGen: ~$0.01 per design × 400 designs = ~$4
   - Boltz-2: ~$0.05 per complex × 400 candidates = ~$20
   - Total: ~$50-100 for a full pipeline run

10. **What is the expected timeline?**
    - GPU compute: 2-4 hours for design + structure prediction
    - Analysis and filtering: <1 hour
    - Total: 4-8 hours for computational phase

11. **How many candidates should go to experimental validation?**
    - Pipeline outputs ~10 candidates
    - Is this sufficient diversity? Should it be 20-30?

---

## Success Criteria and Decision Points

### Computational Success Criteria

| Checkpoint | Success Criterion | Action if Fail |
|------------|-------------------|----------------|
| Calibration | All 3 known binders (teplizumab, SP34, UCHT1) pass with pDockQ > threshold - margin | Investigate Boltz-2 predictions; consider alternative structure prediction |
| De novo design | ≥20 designs pass initial pDockQ filter | Increase design count or relax threshold |
| Filtering | ≥10 candidates survive full filter cascade | Apply fallback (relax soft filters, then thresholds) |
| Diversity | Final candidates include ≥2 VHH, ≥2 scFv, ≥1 OKT3-like epitope | Manually adjust selection to ensure diversity |
| Epitope coverage | ≥3 candidates have OKT3-like epitope (>50% overlap) | Prioritize optimization track; de novo designs may have novel epitopes |

### Experimental Go/No-Go Criteria

| Milestone | Go Criterion | No-Go Criterion |
|-----------|--------------|-----------------|
| Binding confirmation | ≥3/10 candidates bind CD3ε by SPR with Kd < 500 nM | 0/10 candidates bind |
| Functional activity | ≥2 candidates activate T-cells (CD69+ >20% vs. control) | No candidates show T-cell activation |
| Affinity range | Measured Kd values span at least 10-fold range | All candidates have similar affinity (no panel diversity) |
| Developability | ≥5 candidates express at >10 mg/L transient | <3 candidates express adequately |

### Failure Mode Detection

| Failure Mode | How to Detect | Mitigation |
|--------------|---------------|------------|
| BoltzGen doesn't generate CD3 binders | All designs fail binding assays | Use optimization track only |
| Boltz-2 gives false positives | High pDockQ designs don't bind | Add orthogonal filtering (e.g., sequence-based scores) |
| All designs bind same (wrong) epitope | Epitope binning shows single cluster | Redesign against different CD3ε face or use 1SY6 constraints |
| Humanization destroys all binders | All humanized variants lose binding | Use back-mutation variants; accept lower humanness |
| Affinity too high (CRS risk) | All designs have Kd < 10 nM | Generate affinity-attenuated variants with known mutations |
| Affinity too low (no killing) | All designs have Kd > 500 nM | Optimize by affinity maturation (in silico or experimental) |

---

## Known Limitations

This section explicitly states what this computational pipeline **cannot do**. These limitations are fundamental to the current state of computational antibody design and should inform experimental planning.

### 1. Affinity Cannot Be Predicted Computationally

**Critical limitation**: No computational tool can reliably predict absolute Kd values.

- **Boltz-2 pDockQ** is a structural confidence score (how confident the model is in the predicted structure), NOT an affinity predictor
- High pDockQ does NOT mean high affinity; low pDockQ does NOT mean low affinity
- The ~50 nM Kd target is a **design hypothesis based on literature**, not a computationally enforced constraint
- Interface metrics (buried surface area, contact counts) correlate weakly with affinity but have wide error margins

**Implication**: The entire affinity range must be validated experimentally. Computational filtering removes structurally implausible binders, not affinity-inappropriate ones.

### 2. Filter Thresholds Are Not Validated

**Critical limitation**: Default filter thresholds are heuristics from literature, not validated on this specific system.

- pDockQ > 0.5 is a commonly used cutoff, but may exclude valid binders or include non-binders
- Humanness thresholds (OASis > 0.8) are reasonable but not validated for CD3 binders specifically
- The calibration phase (using known binders) addresses this, but calibration quality depends on Boltz-2 accuracy

**Implication**: Always run calibration with known binders before filtering. If known binders fail filters, adjust thresholds.

### 3. Humanization May Destroy Binding

**Critical limitation**: BioPhi/Sapiens humanization suggests mutations to increase humanness, but these may disrupt binding.

- CDR mutations during humanization can eliminate antigen binding
- Framework mutations can alter CDR loop conformations
- No computational method reliably predicts which humanization mutations are safe

**Mitigation strategy**:
1. For each humanized candidate, generate a "back-mutation" variant with original CDR residues restored
2. Include both humanized and back-mutated versions in experimental testing
3. Prioritize humanization mutations in framework regions over CDRs

### 4. De Novo Designs May Have Different Biology

**Critical limitation**: BoltzGen designs may bind epitopes different from the clinically validated OKT3 epitope.

- Novel epitopes may have different T-cell activation profiles
- Novel epitopes may have different safety profiles (CRS, neurotoxicity)
- Cross-reactivity patterns (human vs cynomolgus) are unpredictable

**Mitigation strategy**:
1. Annotate epitope overlap with OKT3 for all de novo designs
2. Prioritize OKT3-like epitopes for initial testing
3. Include novel epitope designs with explicit risk acknowledgment

### 5. Immunogenicity Cannot Be Fully Predicted

**Critical limitation**: T-cell epitope prediction tools (NetMHCIIpan) are excluded due to license restrictions.

- BioPhi OASis scores humanness (similarity to human antibody repertoire)
- Humanness is a proxy for low immunogenicity, but not a guarantee
- Some human-like sequences are still immunogenic; some non-human sequences are tolerated

**Implication**: Immunogenicity must be assessed experimentally or with licensed prediction tools outside this pipeline.

### 6. Developability Is Probabilistic, Not Deterministic

**Critical limitation**: Computational filters reduce risk but don't guarantee success.

- Sequences passing all filters may still aggregate, express poorly, or be unstable
- Sequences failing soft filters may still be developable
- Context matters: CDR liabilities in buried positions may be acceptable

**Implication**: Plan for 50-70% attrition in experimental validation, even for computationally "clean" candidates.

---

## Experimental Validation Requirements

Computational design produces **candidates for experimental testing**, not finished therapeutics. The following experimental validation is required.

### Minimum Viable Validation

| Assay | Purpose | Required For |
|-------|---------|--------------|
| SPR/BLI binding | Confirm binding, measure Kd | All candidates |
| CD3ε-Fc ELISA | Confirm specificity | All candidates |
| T-cell activation (CD69/CD25) | Confirm functional engagement | Top 5 candidates |
| Cytokine release (IL-2, IFN-γ, TNF-α) | Assess CRS potential | Top 5 candidates |
| Expression titer (transient) | Assess manufacturability | All candidates |
| SEC-HPLC | Assess aggregation | All candidates |
| DSF/DSC | Assess thermal stability | Top 5 candidates |

### Validation of Affinity Hypothesis

The ~50 nM Kd hypothesis requires an **affinity panel** experiment:

1. Select 2-3 candidates with varying predicted binding metrics
2. Generate affinity variants (if using optimized track): WT, 10x weaker, 100x weaker
3. Measure Kd for all variants by SPR
4. Measure T-cell killing and cytokine release for all variants
5. Plot Kd vs. killing efficacy and Kd vs. cytokine release
6. Identify optimal Kd range empirically

### Validation of Epitope Effects

If de novo designs with novel epitopes are advanced:

1. Perform epitope binning with OKT3
2. Compare T-cell activation kinetics to OKT3-like binders
3. Assess cross-blocking with endogenous TCR signaling

---

## Reproducibility

### Version Pinning

All tool versions are pinned in `requirements.txt` and `pyproject.toml`:

```
boltzgen==0.2.1
boltz2==0.3.0
anarci==1.3.1
biophi==1.0.5
abnumber==0.3.5
abodybuilder2==1.0.2
```

### Random Seeds

All stochastic operations use explicit seeds documented in `config.yaml`:

```yaml
reproducibility:
  boltzgen_seed: 42
  sampling_seed: 12345
  train_test_split_seed: 0
```

### Modal Container Hashes

GPU compute containers are versioned by hash:

```yaml
modal:
  boltzgen_image_hash: "sha256:abc123..."
  boltz2_image_hash: "sha256:def456..."
```

### Output Provenance

Each output file includes provenance metadata:

```yaml
# In each output YAML/JSON
_provenance:
  pipeline_version: "1.0.0"
  git_commit: "abc1234"
  run_timestamp: "2026-01-15T10:30:00Z"
  config_hash: "sha256:..."
  tool_versions:
    boltzgen: "0.2.1"
    boltz2: "0.3.0"
```

### Reproducing Results

```bash
# Clone at specific commit
git clone https://github.com/alecnielsen/cd3-binder-design.git
cd cd3-binder-design
git checkout <commit_hash>

# Install exact versions
pip install -r requirements.txt --no-deps

# Run with same config
python scripts/run_full_pipeline.py --config config.yaml
```

---

## Technical Considerations

### Core Assumptions

1. **CD3ε is the binding target** - Most bispecifics target CD3ε. Alternative approaches targeting CD3γ or CD3δ are less validated.

2. **~50 nM Kd balances efficacy and safety** - Based on Staflin et al. (2021) and clinical trends. May not hold for all tumor types.

3. **Static structures are sufficient** - Using crystal structures (1XIW, 1SY6) assumes CD3ε conformation in crystals represents the relevant state in vivo.

4. **Sequence → structure → function** - The pipeline assumes that sequence properties predict structure, which predicts function. This chain has uncertainty at each step.

5. **Starting sequences are freely usable** - Teplizumab, SP34, and UCHT1 are from expired patents/public domain. Legal review recommended for commercial use.

### Questions This Pipeline Cannot Answer

| Question | Why Not | Alternative Approach |
|----------|---------|---------------------|
| What is the exact Kd of this binder? | No computational method predicts absolute affinity | SPR/BLI measurement required |
| Will this binder activate T-cells? | Binding ≠ function; signaling depends on geometry, avidity, context | T-cell activation assays required |
| Will this binder cause CRS? | CRS depends on patient, tumor burden, dosing, not just binder properties | Cynomolgus studies, clinical observation |
| Is this sequence immunogenic? | T-cell epitope tools excluded; humanness is a proxy | Licensed tools (NetMHCIIpan) or clinical data |
| Will this express well in CHO? | Sequence features weakly predict expression | Transient expression testing |
| What is the optimal format? | All 5 formats have tradeoffs; optimal depends on indication | Test multiple formats in parallel |
| Will this work in patients? | Preclinical models have limited predictive value | Clinical trials |

### Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| BoltzGen fails to generate binders | Low | High | Optimization track provides fallback; use both target structures |
| All designs fail filters | Medium | High | Calibration phase sets appropriate thresholds; iterative relaxation |
| Humanization destroys binding | Medium | Medium | Include back-mutation variants; test both humanized and original |
| De novo epitopes have poor biology | Medium | High | Annotate epitope overlap; prioritize OKT3-like; include SP34/UCHT1 track |
| Insufficient candidates survive | Medium | Medium | Fallback protocol: relax soft filters, then hard filters with documentation |
| Poor Kd range for CRS balance | High | High | Generate affinity panel; validate hypothesis experimentally |
| Cross-reactivity issues for preclinical | Medium | Medium | Include SP34-derived sequences; test cynomolgus binding early |

### What Could Go Wrong (Detailed Failure Scenarios)

**Scenario 1: Boltz-2 predictions are unreliable for this system**
- *Symptoms*: Known binders score poorly; predicted interfaces don't match known epitopes
- *Detection*: Calibration phase (step 0)
- *Response*: Fall back to sequence-based filtering (humanness, liabilities) without structural metrics; increase reliance on optimization track

**Scenario 2: All de novo designs bind wrong epitope**
- *Symptoms*: Epitope overlap with OKT3 is <20% for all designs
- *Detection*: Epitope annotation in filtering step
- *Response*: Constrain BoltzGen to OKT3 epitope region (if supported); rely on optimization track; test novel epitopes with explicit risk acknowledgment

**Scenario 3: Affinity panel shows all binders are too tight (CRS risk)**
- *Symptoms*: All measured Kd values are <10 nM
- *Detection*: SPR measurements in experimental validation
- *Response*: Apply known affinity-attenuating mutations from literature; generate additional variants

**Scenario 4: Humanized sequences all fail binding**
- *Symptoms*: 0% retention of binding after BioPhi humanization
- *Detection*: Binding assays on humanized vs. original sequences
- *Response*: Use back-mutation variants; apply conservative humanization (fewer mutations); accept lower humanness score

**Scenario 5: Expression failures across all formats**
- *Symptoms*: <1 mg/L transient expression for all candidates
- *Detection*: Transient expression testing
- *Response*: Screen additional candidates; investigate sequence features correlated with expression; consider different expression system

### Computational Resources

| Step | GPU Type | Time (est.) | Cost (est.) |
|------|----------|-------------|-------------|
| BoltzGen (400 designs) | A100 | 1-2 hours | $5-15 |
| Boltz-2 (400 complexes) | A100 | 2-4 hours | $10-30 |
| ABodyBuilder2 | A100 (optional) | 30 min | $2-5 |
| **Total** | | **4-8 hours** | **$20-50** |

Local execution requires NVIDIA GPU with CUDA support. Apple MPS is not compatible with BoltzGen.

See [Known Limitations](#known-limitations) for detailed discussion of fundamental constraints.

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
