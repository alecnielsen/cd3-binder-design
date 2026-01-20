# Scientific Review: Structure Module

You are reviewing the `src/structure/` and `modal/` modules of a CD3 binder design pipeline. These modules handle:
- Antibody structure prediction (ABodyBuilder2)
- Complex structure prediction (Boltz-2)
- Interface analysis (contact residues, buried surface area, pDockQ)
- Modal GPU deployment for compute-intensive predictions

## Important: Module README Updates

The scientific context for this module lives in `src/structure/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. Boltz-2 Complex Prediction

Boltz-2 predicts the structure of binder-antigen complexes.

**Key outputs:**
- Predicted complex structure (PDB format)
- pDockQ score (structural confidence metric)
- Interface residues and contacts

**Critical caveat from README:**
> "pDockQ is a structural confidence score (how confident the model is in the predicted structure), NOT an affinity predictor"
> "High pDockQ does NOT mean high affinity; low pDockQ does NOT mean low affinity"

**Verify:**
- Is pDockQ correctly extracted from Boltz-2 output?
- Are there any comments or variable names that incorrectly suggest pDockQ predicts affinity?
- Is the interpretation documented correctly in output metadata?

### 2. Interface Analysis Correctness

Interface analysis extracts binding site information:

**Metrics to verify:**
- **Contact residues**: Pairs of residues within a distance cutoff (typically 4-5 A for heavy atoms)
- **Buried surface area (BSA)**: Solvent-accessible surface area buried upon complex formation
- **Contact count**: Number of residue-residue contacts across the interface

**Verify:**
- Distance cutoff is scientifically appropriate (4-5 A for contacts)
- BSA calculation uses correct atom radii and probe radius (typically 1.4 A)
- Contact counting is residue-level, not atom-level (as noted in CLAUDE.md)

**From CLAUDE.md:**
> "Contact counts are residue-level - `num_contacts` counts unique residue pairs, not atomic contacts"

### 3. Calibration Phase Logic

The calibration phase uses known binders (teplizumab, SP34, UCHT1) to set filter thresholds.

**Calibration logic:**
1. Run Boltz-2 on known binder + CD3 epsilon complexes
2. Record pDockQ, interface area, contact count for each
3. Set thresholds = minimum(known binders) - margin

**Verify:**
- All three known binders are used for calibration
- Margins are applied correctly (subtract from minimum, not maximum)
- If any known binder fails, it indicates a problem with the prediction model, not the binder

**Critical principle from README:**
> "If known binders fail default thresholds, thresholds are wrong, not the binders"

### 4. pDockQ Interpretation

pDockQ is Boltz-2's confidence score for the predicted complex.

**Score interpretation:**
- pDockQ > 0.5: Higher confidence in predicted interface
- pDockQ < 0.3: Lower confidence, prediction may be unreliable

**What pDockQ does NOT indicate:**
- Binding affinity (Kd)
- Whether the binder will work experimentally
- Relative ranking of binders by affinity

**Verify:**
- Code does not conflate pDockQ with affinity
- Filtering uses pDockQ as a structural quality filter, not affinity predictor
- Documentation is clear about this distinction

### 5. Modal Deployment

GPU-intensive operations run on Modal (cloud GPU platform).

**Verify:**
- GPU type is appropriate (A100 recommended for Boltz-2)
- Timeout is sufficient for complex prediction (can take 5-10 minutes per complex)
- Error handling for Modal failures (network, GPU OOM, timeout)
- Reproducibility: seeds are passed to Modal functions

### 6. ABodyBuilder2 Usage

ABodyBuilder2 predicts antibody Fv structures (VH/VL or VHH alone).

**Verify:**
- Correct chain types are passed (VH, VL, or VHH)
- For scFvs, the linker is handled appropriately
- Output structures are in standard PDB format

### 7. Scientific Assumptions to Challenge

**From Assumption 2 (Boltz-2):**
> "Antibody-antigen complexes may be harder than general protein complexes"
> "Boltz-2 was not specifically benchmarked on VHH/scFv-antigen complexes"
> "The model may hallucinate interfaces that look plausible but don't exist"

**Questions:**
- Does the code have any validation that predicted interfaces are sensible?
- Is there comparison to the known OKT3 epitope (from 1SY6) as a sanity check?
- What happens if calibration reveals known binders score poorly?

## Output Format

If you find **code issues**, list them with:
1. File path and line number
2. Description of the issue
3. Scientific rationale for why it matters
4. Suggested fix (if applicable)

If no **code issues** are found, respond with exactly:
NO_ISSUES

**Important notes:**
- Do not invent issues if the code is correct. Only report genuine scientific or implementation problems.
- README updates do NOT count as issues. If you update the module README.md to fix documentation, still output "NO_ISSUES" if no code problems were found.
- "NO_ISSUES" means the code is scientifically correct, even if you made documentation improvements.
