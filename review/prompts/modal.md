# Scientific Review: Modal Module

You are reviewing the `modal/` directory of a CD3 binder design pipeline. This module contains Modal GPU deployments for:
- BoltzGen de novo binder design
- Boltz-2 protein complex structure prediction
- ABodyBuilder2 antibody structure prediction (optional)

## Important: Module README Updates

The scientific context for this module lives in `modal/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. BoltzGen Configuration (boltzgen_app.py)

**Verify:**
- Container image includes all required dependencies
- GPU type is appropriate (A100 recommended)
- Timeout is sufficient for batch operations
- Random seed is passed through for reproducibility

**Parameters to check:**
- `binder_type`: Should support "vhh" and "scfv"
- `num_designs`: Default should be reasonable (100-200)
- `hotspot_residues`: Optional epitope targeting
- `temperature`: Sampling diversity control

### 2. Boltz-2 Configuration (boltz2_app.py)

**Verify:**
- Container image includes Boltz dependencies
- Complex prediction handles binder + target sequences
- Interface metrics are calculated correctly:
  - `num_contacts`: Unique residue pairs within cutoff
  - `interface_area`: Estimated buried surface area
  - `interface_residues_target`: CD3 contact residues

**Critical distinction:**
> "pDockQ is a structural confidence score (how confident the model is in the predicted structure), NOT an affinity predictor. High pDockQ does NOT mean high affinity."

### 3. Interface Metric Calculations

**Verify contact counting:**
- Distance cutoff is reasonable (5.0 A is standard)
- Contacts count unique residue pairs, not atomic contacts
- Interface area estimation uses reasonable per-residue contribution

**From CLAUDE.md:**
> "Contact counts are residue-level - `num_contacts` counts unique residue pairs, not atomic contacts"

### 4. Calibration Function (boltz2_app.py)

**Verify:**
- `run_calibration()` processes known binders correctly
- Calibrated thresholds are calculated with margin:
  - `min_pdockq = min(known_binder_pdockq) - 0.05`
  - `min_interface_area = min(known_binder_area) - 100`
  - `min_contacts = max(0, min(known_binder_contacts) - 2)`
- Statistics (min, max, mean) are returned for inspection

### 5. GPU Configuration

**Verify:**
- GPU type matches computational requirements:
  - BoltzGen: A100 (40GB) recommended
  - Boltz-2: A100 (40GB) recommended
- Timeout is sufficient:
  - Single design: 30 min - 1 hour
  - Batch operations: 2+ hours
- Retries are configured for transient failures

### 6. Error Handling

**Verify:**
- GPU unavailability produces clear error
- Timeout errors are distinguishable from other failures
- Batch operations continue on single-item failure
- Failed items are returned with error information

### 7. Reproducibility

**Verify:**
- Random seeds are used consistently
- Batch operations increment seeds appropriately (e.g., `seed + i`)
- Container versions can be pinned (image hashes)

### 8. PDB Parsing and Sequence Extraction

**Verify `parse_pdb_sequence()`:**
- Handles standard amino acids correctly
- Filters by chain ID
- Handles insertion codes and gaps
- Returns sequences in residue number order

### 9. Scientific Assumptions to Challenge

**BoltzGen assumptions:**
- Can BoltzGen design VHH/scFv that bind CD3 epsilon?
- Are 100-200 designs sufficient for hit discovery?
- Does hotspot targeting improve binding to OKT3 epitope?

**Boltz-2 assumptions:**
- Are complex predictions reliable for VHH/scFv-antigen systems?
- Does pDockQ correlate with actual binding?
- Are interface metrics predictive of binding quality?

**Questions:**
- Is the interface area estimation reasonable (~80 A^2 per residue)?
- Should the distance cutoff be configurable?
- Are there edge cases in PDB parsing (modified residues, multiple models)?

### 10. Licensing Verification

All tools must have permissive licenses:
- BoltzGen: MIT (verify)
- Boltz-2: MIT (verify)
- ABodyBuilder2: BSD 3-Clause (verify)

**From README:**
> "All tools must have permissive licenses (MIT, BSD, Apache, CC-BY) for commercial use"

## Output Format

If you find **code issues**:
1. **FIX THEM** by editing the code files directly
2. Then report what you fixed:
   - File path and line number
   - Description of what was wrong
   - What you changed and why
3. Do NOT say "NO_ISSUES" - another iteration will verify your fixes

If you reviewed the code and found **no issues to fix**, respond with exactly:
NO_ISSUES

**Important notes:**
- Only say "NO_ISSUES" if you made zero edits this iteration
- If you fixed anything, report it and let the next iteration verify
- Do not invent issues if the code is correct. Only report genuine scientific or implementation problems.
- README updates do NOT count as issues. If you update the module README.md to fix documentation, still output "NO_ISSUES" if no code problems were found.
