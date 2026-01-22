# CD3 Binder Design Pipeline

A computational pipeline for designing CD3-binding domains (VHH nanobodies and scFvs) suitable for use in bispecific antibody therapeutics. This project uses open-source, permissively-licensed tools to generate, filter, and format CD3 binders with drug-like properties.

## Table of Contents

1. [Project Overview](#project-overview)
2. [Scientific Background](#scientific-background)
3. [Design Strategy](#design-strategy)
4. [Computational Pipeline](#computational-pipeline)
5. [Tool Stack & Licensing](#tool-stack--licensing)
6. [Installation & Usage](#installation--usage)
7. [Known Limitations](#known-limitations)
8. [Experimental Validation Requirements](#experimental-validation-requirements)
9. [Critical Assumptions for Scientific Review](#critical-assumptions-for-scientific-review)
10. [Success Criteria and Decision Points](#success-criteria-and-decision-points)
11. [References](#references)

> **Additional documentation**: See `docs/reference/` for detailed implementation examples, bispecific format diagrams, and repository structure.

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

**Design hypothesis**: ~50 nM Kd may be optimal, with variants spanning 10-500 nM range. However, computational tools cannot predict absolute Kd values (see [Known Limitations](#known-limitations)). This hypothesis must be validated experimentally across an affinity panel.

### Source Structures

Two PDB structures are used as design targets:

1. **1XIW**: CD3εδ heterodimer - provides the CD3ε surface for de novo design
2. **1SY6**: CD3εγ heterodimer in complex with OKT3 Fab - shows validated binding mode

### Starting Sequences

All starting sequences are from expired patents or public domain publications:

| Antibody | Origin | Cross-Reactivity | Notes |
|----------|--------|------------------|-------|
| **Teplizumab** | Humanized OKT3 | Human only | FDA approved 2022 for T1D |
| **SP34** | Murine, 1985 | Human + Cynomolgus | Important for preclinical |
| **UCHT1** | Murine, 1981 | Human + Chimpanzee | Alternative epitope |

---

## Design Strategy

### Track 1: De Novo Design (BoltzGen)

BoltzGen generates VHH nanobodies (~120 aa) and scFvs (~250 aa) targeting CD3ε from 1XIW and 1SY6 structures.

**Limitations:**
- Cannot design full-length IgG directly
- Affinity is not explicitly controlled
- Novel epitopes may differ from validated OKT3 site

### Track 2: Optimization of Existing Binders

Starting from known anti-CD3 CDR sequences:
1. Graft CDRs onto human frameworks
2. Apply humanization via BioPhi/Sapiens
3. Generate affinity-attenuated variants (WT, 10x weaker, 100x weaker)

### Track 3: Hybrid Approach

Use BoltzGen to design VH or VL individually, then pair with known complementary chains.

---

## Computational Pipeline

```
INPUT PREPARATION
├── Download CD3 structures (1XIW, 1SY6)
├── Extract target chains
└── Prepare starting sequences

DESIGN GENERATION
├── De Novo (BoltzGen on Modal GPU)
│   └── VHH + scFv designs (100+ each)
└── Optimization Track
    └── Humanize + generate affinity variants

STRUCTURE PREDICTION
├── ABodyBuilder2 (antibody structure)
└── Boltz-2 (complex with CD3ε)

CALIBRATION (run first!)
└── Set filter thresholds using known binders

FILTERING CASCADE
├── Binding quality (pDockQ, interface area, contacts)
├── Humanness (OASis > 0.8)
├── Liabilities (deamidation, glycosylation, etc.)
├── Developability (charge, CDR-H3 length, etc.)
└── Aggregation propensity

FORMAT CONVERSION
└── 5 bispecific formats (CrossMab, Fab+scFv, etc.)

OUTPUT: ~10 candidates with sequences, structures, scorecards
```

---

## Tool Stack & Licensing

All tools have permissive licenses suitable for commercial use:

| Tool | Purpose | License |
|------|---------|---------|
| **BoltzGen** | De novo binder design | MIT |
| **Boltz-2** | Complex structure prediction | MIT |
| **ABodyBuilder2** | Antibody structure prediction | BSD 3-Clause |
| **ANARCI** | Antibody numbering | BSD 3-Clause |
| **BioPhi/Sapiens** | Humanness scoring | MIT |

**Excluded (non-permissive):** IgFold (JHU Academic), NetMHCIIpan (DTU academic)

**Compute:** Modal provides on-demand A100/H100 GPUs for BoltzGen and Boltz-2.

---

## Installation & Usage

### Prerequisites

- Python 3.10+
- Modal account: https://modal.com

### Installation

```bash
git clone https://github.com/alecnielsen/cd3-binder-design.git
cd cd3-binder-design
pip install -e .

# Modal setup
pip install modal
modal setup
modal deploy modal/boltzgen_app.py
modal deploy modal/boltz2_app.py
```

### Running the Pipeline

```bash
# Step 0: Calibration (RUN FIRST)
python scripts/00_run_calibration.py

# Full pipeline
python scripts/run_full_pipeline.py --config config.yaml

# Or step-by-step
python scripts/01_setup_targets.py      # Download structures
python scripts/02_run_denovo_design.py  # BoltzGen (Modal)
python scripts/03_run_optimization.py   # Humanization
python scripts/04_predict_structures.py # Boltz-2 (Modal)
python scripts/05_filter_candidates.py  # Apply filters
python scripts/06_format_bispecifics.py # Generate formats
python scripts/07_generate_report.py    # Final report
```

---

## Known Limitations

This section explicitly states what this computational pipeline **cannot do**.

### 1. Affinity Cannot Be Predicted Computationally

**Critical**: No computational tool can reliably predict absolute Kd values.

- **Boltz-2 pDockQ** is a structural confidence score, NOT an affinity predictor
- High pDockQ does NOT mean high affinity
- The ~50 nM Kd target is a design hypothesis, not a computationally enforced constraint

**Implication**: The entire affinity range must be validated experimentally.

### 2. Filter Thresholds Are Not Validated

Default thresholds are heuristics from literature, not validated on this specific system. The calibration phase uses known binders to set appropriate thresholds.

**Implication**: Always run calibration before filtering.

### 3. Humanization May Destroy Binding

BioPhi/Sapiens mutations may disrupt binding. The pipeline generates back-mutation variants as a safety net.

### 4. De Novo Designs May Have Different Biology

BoltzGen designs may bind epitopes different from the clinically validated OKT3 epitope, with unknown T-cell activation and safety profiles.

### 5. Immunogenicity Cannot Be Fully Predicted

T-cell epitope prediction tools (NetMHCIIpan) are excluded due to license restrictions. BioPhi OASis scores humanness as a proxy.

### 6. Developability Is Probabilistic

Computational filters reduce risk but don't guarantee success. Plan for 50-70% attrition in experimental validation.

---

## Experimental Validation Requirements

Computational design produces **candidates for experimental testing**, not finished therapeutics.

### Minimum Viable Validation

| Assay | Purpose | Required For |
|-------|---------|--------------|
| SPR/BLI binding | Confirm binding, measure Kd | All candidates |
| CD3ε-Fc ELISA | Confirm specificity | All candidates |
| T-cell activation (CD69/CD25) | Confirm functional engagement | Top 5 candidates |
| Cytokine release (IL-2, IFN-γ, TNF-α) | Assess CRS potential | Top 5 candidates |
| Expression titer (transient) | Assess manufacturability | All candidates |
| SEC-HPLC | Assess aggregation | All candidates |

### Validation of Affinity Hypothesis

The ~50 nM Kd hypothesis requires an **affinity panel** experiment:
1. Select 2-3 candidates with varying predicted binding metrics
2. Generate affinity variants: WT, 10x weaker, 100x weaker
3. Measure Kd, T-cell killing, and cytokine release for all
4. Identify optimal Kd range empirically

---

## Critical Assumptions for Scientific Review

Reviewers should challenge these assumptions:

### Assumption 1: BoltzGen Can Design Functional CD3 Binders

- **For**: Validated on other targets
- **Against**: Not specifically validated on CD3 or immune receptors
- **Fallback**: Rely on optimization track (teplizumab, SP34, UCHT1)

### Assumption 2: Boltz-2 Complex Predictions Are Meaningful

- **For**: High accuracy on protein complex benchmarks
- **Against**: May not generalize to antibody-antigen complexes
- **Detection**: Calibration phase should catch failures

### Assumption 3: ~50 nM Kd Balances Efficacy and Safety

- **For**: Literature supports this range (Staflin et al. 2021)
- **Against**: Optimal Kd may be context-dependent
- **Key question**: Pipeline cannot enforce 50 nM computationally

### Assumption 4: Humanization Preserves Function

- **For**: BioPhi validated on therapeutic antibodies
- **Against**: 30-70% success rate, not 100%
- **Mitigation**: Back-mutation variants

### Assumption 5: Developability Filters Are Predictive

- **For**: Literature supports correlation
- **Against**: Probabilistic, not deterministic
- **Detection**: Experimental validation

---

## Success Criteria and Decision Points

### Computational Success Criteria

| Checkpoint | Success Criterion | Action if Fail |
|------------|-------------------|----------------|
| Calibration | All 3 known binders pass | Investigate Boltz-2 predictions |
| De novo design | ≥20 designs pass pDockQ filter | Increase design count |
| Filtering | ≥10 candidates survive | Apply fallback (relax filters) |
| Diversity | Mix of VHH/scFv, OKT3-like/novel | Adjust selection |

### Experimental Go/No-Go Criteria

| Milestone | Go Criterion | No-Go Criterion |
|-----------|--------------|-----------------|
| Binding confirmation | ≥3/10 bind with Kd < 500 nM | 0/10 bind |
| Functional activity | ≥2 activate T-cells | None activate |
| Affinity range | Kd spans ≥10-fold range | All similar affinity |
| Developability | ≥5 express at >10 mg/L | <3 express |

### Failure Mode Detection

| Failure Mode | Detection | Mitigation |
|--------------|-----------|------------|
| BoltzGen doesn't generate binders | All fail binding assays | Use optimization track only |
| Boltz-2 gives false positives | No correlation pDockQ vs binding | Add sequence-based scores |
| All designs bind wrong epitope | Epitope binning shows single cluster | Redesign with constraints |
| Humanization destroys all binders | All humanized variants lose binding | Use back-mutation variants |

---

## References

### CD3 Biology and Bispecifics
1. Wolf E, et al. (2005). BiTEs: bispecific antibody constructs. Drug Discovery Today.
2. Goebeler ME, Bargou RC (2020). T cell-engaging therapies. Nature Reviews Clinical Oncology.

### CD3 Affinity and Safety
3. Staflin K, et al. (2021). CD3 affinity tuning. Scientific Reports. https://www.nature.com/articles/s41598-021-93842-0

### Source Antibodies
4. Kung P, et al. (1979). OKT3. Science.
5. Pessano S, et al. (1985). SP34. EMBO J.
6. Callard RE, et al. (1981). UCHT1. Clin Exp Immunol.

### Computational Tools
7. Abanades B, et al. (2023). ImmuneBuilder. Communications Biology.
8. Dunbar J, Deane CM (2016). ANARCI. Bioinformatics.
9. Prihoda D, et al. (2022). BioPhi. mAbs.
10. BoltzGen (2024). MIT Jameel Clinic.

### Bispecific Formats
11. Brinkmann U, Kontermann RE (2017). Bispecific antibodies. mAbs.
12. Schaefer W, et al. (2011). CrossMab. PNAS.

---

## License

MIT License - See LICENSE file for details.
