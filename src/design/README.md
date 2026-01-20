# Design Module - Scientific Context

This module handles:
- De novo binder design using BoltzGen
- Optimization of existing binders (teplizumab, SP34, UCHT1)
- Generation of affinity-attenuated variants

## Design Strategy

### Track 1: De Novo Design (BoltzGen)

BoltzGen is an open-source (MIT license) generative model for protein binder design. It can design:
- VHH nanobodies (~120 amino acids, single domain)
- scFvs (~250 amino acids, VH-linker-VL)

**Approach:**
1. Provide CD3 epsilon structure as target (from 1XIW and 1SY6)
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

## Design Generation Pipeline

```
De Novo (BoltzGen)              Optimization Track
- VHH designs (100+)            - Humanize SP34, UCHT1 CDRs
- scFv designs (100+)           - Graft onto human frameworks
- Target: 1XIW, 1SY6            - Generate affinity variants
- Run on Modal (GPU)              (WT, 10x weaker, 100x weaker)
                                - Combinatorial mutations
```

## Configuration

```yaml
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
```

## Critical Assumptions

### Assumption 1: BoltzGen Can Design Functional CD3 Binders

**Claim**: BoltzGen can generate VHH/scFv sequences that bind CD3 epsilon with reasonable affinity.

**Evidence supporting this**:
- BoltzGen has been validated on other protein targets in published benchmarks
- The model was trained on protein-protein interface data that includes antibody-antigen complexes

**Evidence against / Unknowns**:
- BoltzGen has not been specifically validated on CD3 or T-cell surface proteins
- CD3 epsilon forms a complex (with delta or gamma); the binding site may require both subunits
- No published data on BoltzGen success rates for immune receptor targets
- De novo designs may have very low hit rates (<1% functional binders)

**What would falsify this assumption**:
- Zero designs pass calibration thresholds after filtering
- All designs fail to bind in SPR/ELISA
- Designs bind but don't activate T-cells

**Fallback if assumption fails**: Rely entirely on optimization track (teplizumab, SP34, UCHT1 variants).

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

### Known Limitation: Affinity Cannot Be Predicted Computationally

**Critical limitation**: No computational tool can reliably predict absolute Kd values.

- The ~50 nM Kd target is a **design hypothesis based on literature**, not a computationally enforced constraint
- Interface metrics (buried surface area, contact counts) correlate weakly with affinity but have wide error margins

**Implication**: The entire affinity range must be validated experimentally. Computational filtering removes structurally implausible binders, not affinity-inappropriate ones.

## Starting Sequences

All starting sequences are from expired patents or public domain publications:

| Antibody | Origin | Cross-Reactivity | Notes |
|----------|--------|------------------|-------|
| **Teplizumab** | Humanized OKT3 | Human only | Already humanized, FDA approved |
| **SP34** | Murine, 1985 | Human + Cynomolgus | Important for preclinical |
| **UCHT1** | Murine, 1981 | Human + Chimpanzee | Alternative epitope |

## Implementation Notes

- De novo and optimization tracks are clearly separated
- Outputs are annotated with their source track
- Humanization is applied only to murine sequences (not teplizumab)
- Seeds must be passed for reproducibility
- Both target structures (1XIW, 1SY6) should be used for design diversity
