# Structure Module - Scientific Context

This module handles:
- Antibody structure prediction (ABodyBuilder2)
- Complex structure prediction (Boltz-2)
- Interface analysis (contact residues, buried surface area, pDockQ)
- Modal GPU deployment for compute-intensive predictions

## Structure Prediction Pipeline

**Step 1: Antibody Structure (ABodyBuilder2, BSD license)**
- Predict VHH or Fv structure
- High accuracy for CDR loops

**Step 2: Complex Structure (Boltz-2, MIT license)**
- Predict binder + CD3 epsilon complex
- Assess binding mode and interface
- Extract pDockQ confidence score

## Calibration Phase

**Purpose:** Establish filter thresholds using KNOWN BINDERS

1. Run Boltz-2 on teplizumab/SP34/UCHT1 + CD3 epsilon complexes
2. Record pDockQ, interface area, contact counts for each
3. Set thresholds to ACCEPT all known binders:
   - min_pdockq = min(known_binder_pdockq) - 0.05
   - min_interface_area = min(known_binder_area) - 100
   - min_contacts = min(known_binder_contacts) - 2
4. Document calibration results in config.yaml

**CRITICAL:** If known binders fail default thresholds, thresholds are wrong, not the binders. Calibration prevents false negative filtering.

## Configuration

```yaml
calibration:
  positive_controls:
    - teplizumab
    - sp34
    - ucht1
  calibrated_thresholds:
    min_pdockq: null              # Set by calibration
    min_interface_area: null      # Set by calibration
    min_contacts: null            # Set by calibration
  calibration_margin:
    pdockq: 0.05                  # Subtract from min known binder
    interface_area: 100           # Subtract from min known binder
    contacts: 2                   # Subtract from min known binder
```

## pDockQ Interpretation

**Critical caveat:** pDockQ is a structural confidence score, NOT an affinity predictor.

- pDockQ > 0.5: Higher confidence in predicted interface
- pDockQ < 0.3: Lower confidence, prediction may be unreliable
- High pDockQ does NOT mean high affinity
- Low pDockQ does NOT mean low affinity

Filtering uses pDockQ as a structural quality filter, not an affinity predictor.

## Critical Assumptions

### Assumption 2: Boltz-2 Complex Predictions Are Meaningful

**Claim**: Boltz-2 can predict binder-CD3 epsilon complex structures accurately enough for filtering.

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

## Modal Deployment for BoltzGen

```python
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

## Interface Analysis

**Metrics:**
- **Contact residues**: Pairs of residues within a distance cutoff (5 Å default, appropriate for heavy atom contacts)
- **Interface area (estimated)**: Approximated as ~80 Å² per interface residue. This is NOT true BSA (buried surface area) calculated from SASA. For accurate BSA, use FreeSASA or similar tools.
- **Contact count**: Number of unique residue-residue contacts across the interface (residue-level, not atomic)

**Important:** Contact counts are residue-level - `num_contacts` counts unique residue pairs, not atomic contacts.

**Note on interface area:** The simplified estimate of 80 Å² per residue is within literature range (60-100 Å²) and sufficient for relative comparisons between candidates. Absolute values should not be compared to crystallographic BSA measurements.

## Implementation Notes

- For VH/VL pairs, calibration constructs scFv (not VH-only) for accurate threshold setting
- Seeds are passed to Modal functions for reproducibility
- GPU type: A100 recommended for Boltz-2
- Timeout should be sufficient for complex prediction (can take 5-10 minutes per complex)
- Error handling needed for Modal failures (network, GPU OOM, timeout)
