# Calibration Methodology

## Purpose

Calibration establishes filter thresholds by running Boltz-2 structure prediction on **known CD3 binders** with experimentally validated binding. The logic: de novo designs that score worse than validated binders on structural metrics are likely poor candidates.

## Important Clarifications

### What We're Measuring

**These are STRUCTURAL CONFIDENCE metrics, not binding affinity predictors:**

| Metric | Meaning | NOT |
|--------|---------|-----|
| **pTM** | Confidence in overall fold | Binding strength |
| **pLDDT** | Per-residue position confidence | Stability |
| **ipTM** | Interface prediction confidence | Affinity (Kd) |
| **num_contacts** | Residue pairs within 5Å | Binding energy |
| **interface_area** | Estimated buried surface | Binding free energy |

**pDockQ** is NOT a native Boltz-2 metric—it's an AlphaFold-Multimer post-processing score. Boltz-2 always returns pDockQ = 0.0, so the binding filter only applies pDockQ thresholds when both the threshold and value are non-zero. Interface area and contacts are the primary binding quality filters.

Boltz-2 does have a separate **affinity prediction module** that outputs log10(IC50), but this is not validated for antibody-antigen systems and is not currently used.

### What We're Testing

The calibration uses **scFv constructs derived from parental antibodies**, not full IgG:

| Name in Code | Actual Construct | Parental Antibody |
|--------------|------------------|-------------------|
| `teplizumab` | scFv (VH-linker-VL) | Teplizumab (humanized OKT3) |
| `sp34` | scFv (VH-linker-VL) | SP34 monoclonal |
| `ucht1` | scFv (VH-linker-VL) | UCHT1 monoclonal |

The scFv linker used is `GGGGSGGGGSGGGGS` (G4S)3.

### Why scFv Instead of Full IgG?

1. **Boltz-2 input constraints**: Complex prediction works best with smaller systems
2. **Consistency with designs**: BoltzGen outputs VHH (~120 aa) or scFv (~250 aa), not full IgG
3. **Binding domain equivalence**: The VH/VL domains contain all CDRs responsible for binding

## Known Binders and Expected Affinities

| Antibody | Reported Kd | Source | Notes |
|----------|-------------|--------|-------|
| Teplizumab | ~5-10 nM | Herold et al. | Humanized, FDA-approved |
| SP34 | ~1-10 nM | Pessano et al. 1985 | Cross-reacts with cyno CD3 |
| UCHT1 | ~1 nM | Beverley & Callard 1981 | High affinity |

**Important**: These affinities are for the parental IgG antibodies. The scFv constructs may have different affinities due to:
- Loss of avidity (monovalent vs bivalent)
- Linker effects on domain orientation
- Absence of Fc-mediated stabilization

## Calibration Results Interpretation

### Example Results (2026-01-26)

| Binder | pTM | pLDDT | Contacts | Interface Area |
|--------|-----|-------|----------|----------------|
| Teplizumab-scFv | 0.944 | 96.4 | 42 | 2560 Å² |
| SP34-scFv | 0.753 | 91.7 | 30 | 2160 Å² |
| UCHT1-scFv | 0.784 | 94.1 | 32 | 2240 Å² |

### Why Teplizumab Scores Higher

Teplizumab is a **humanized** antibody with a framework optimized for stability. Boltz-2's training data includes many human antibody structures, so it's more confident predicting human-like sequences. This does NOT mean teplizumab binds better than UCHT1 (which actually has higher reported affinity).

### Why SP34 Scores Lower on pTM

SP34 is a **murine** antibody with a less human-like framework. Lower pTM reflects Boltz-2's reduced confidence in predicting murine antibody folds, not weaker binding.

## Threshold Calculation

Thresholds are set as: `minimum observed value - margin`

```
min_contacts = min(42, 30, 32) - 2 = 28
min_interface_area = min(2560, 2160, 2240) - 100 = 2060 Å²
```

This ensures that de novo designs must score at least as well as the weakest-scoring known binder (with a small tolerance).

## Limitations

### 1. No Affinity Validation

We cannot validate that Boltz-2 metrics correlate with actual binding affinity. The assumption is that structural confidence is a necessary (but not sufficient) condition for binding.

### 2. scFv ≠ IgG

Calibration uses scFv constructs. Results may not translate directly to full IgG or other formats.

### 3. Single Target Structure

Calibration uses CD3ε from 1XIW (chain A). Different CD3 conformations or epitopes are not tested.

### 4. Interface Metrics Are Approximate

- **interface_area** is estimated as `(binder_residues + target_residues) × 80 Å²`
- True buried surface area requires solvent-accessible surface calculation
- **num_contacts** counts residue pairs, not atomic contacts

## Boltz-2 Affinity Prediction (Not Currently Used)

Boltz-2 has an affinity prediction module that outputs **log10(IC50)** values. This could theoretically be enabled with:

```bash
boltz predict input.yaml --sampling_steps_affinity 200 --diffusion_samples_affinity 5
```

And adding to the YAML:
```yaml
properties:
  affinity:
    ligand: B  # binder chain
```

**Why we don't use it:**
1. Not validated for antibody-antigen systems
2. IC50 is assay-dependent (not intrinsic Kd)
3. Would require separate validation against SPR/BLI data

## Validation Baselines

When `calibration.run_validation_baselines: true` (default), step 00 also runs ProteinMPNN, AntiFold, and Protenix on the known binder controls. This establishes baselines for the validation scores used in ranking.

**Process:**
1. Save Boltz-2 CIF files for controls to `data/outputs/calibration_cif/`
2. Run `batch_score_affinity()` on control CIFs (ProteinMPNN + AntiFold, local CPU)
3. Run Protenix `predict_complex()` on controls via Modal H100
4. Store baselines in `calibration.json` under `validation_baselines`:

```json
{
  "validation_baselines": {
    "proteinmpnn_ll": {"min": ..., "max": ..., "mean": ..., "by_control": {...}},
    "antifold_ll": {"min": ..., "max": ..., "mean": ..., "by_control": {...}},
    "protenix_iptm": {"min": ..., "max": ..., "mean": ..., "by_control": {...}}
  }
}
```

These baselines provide context for interpreting de novo design scores — a design scoring better than known binders on ProteinMPNN/AntiFold log-likelihood is a positive signal (though not a guarantee of binding).

## Recommendations

1. **Always run calibration** before filtering de novo designs
2. **Don't over-interpret** high pTM/pLDDT as strong binding
3. **Use calibrated thresholds** to filter candidates
4. **Review validation baselines** to contextualize de novo design scores
5. **Validate experimentally** with SPR/BLI for actual affinity
6. **Consider enabling affinity prediction** if validation data becomes available
