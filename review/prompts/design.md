# Scientific Review: Design Module

You are reviewing the `src/design/` module of a CD3 binder design pipeline. This module handles:
- De novo binder design using BoltzGen
- Optimization of existing binders (teplizumab, SP34, UCHT1)
- Generation of affinity-attenuated variants

## Important: Module README Updates

The scientific context for this module lives in `src/design/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. BoltzGen Integration Correctness

BoltzGen is used to generate de novo VHH and scFv designs against CD3 epsilon.

**Verify:**
- Target structure is provided correctly (PDB format, CD3 epsilon chain)
- Design type is specified correctly ("vhh" vs "scfv")
- Seed is passed for reproducibility
- Output sequences are validated (correct length, no invalid characters)

**Expected behavior:**
- VHH designs: ~120 amino acids, single domain
- scFv designs: ~250 amino acids, VH-linker-VL format

**Questions:**
- Does the code validate that generated sequences are plausible antibodies?
- Is there handling for BoltzGen failures or timeouts?

### 2. Affinity Attenuation Strategy

The pipeline aims for ~50 nM Kd based on literature suggesting this balances efficacy and reduced CRS.

**Known affinity-attenuating mutations for OKT3-derived binders:**
- Mutations should be sourced from peer-reviewed literature
- Typical strategy: mutate CDR contact residues to reduce binding

**Verify:**
- Are literature-derived mutations correctly implemented?
- Are mutation positions specified in a consistent numbering scheme (IMGT preferred)?
- Does the code generate a panel spanning the intended affinity range (WT, 10x weaker, 100x weaker)?

**Critical caveat from README:**
> "No computational tool can reliably predict absolute Kd values."

**Questions:**
- Does the code make any claims about predicted affinity?
- Is it clear that affinity must be validated experimentally?

### 3. Starting Sequence Handling

The pipeline uses three starting sequences:
1. **Teplizumab** (humanized OKT3) - already humanized, FDA approved
2. **SP34** - murine, cross-reacts with cynomolgus CD3
3. **UCHT1** - murine, alternative epitope

**Verify:**
- Sequences are loaded from correct YAML files
- Chain types (VH, VL) are correctly identified
- Humanization is applied only to murine sequences (not teplizumab)

### 4. De Novo vs Optimization Track

**De novo (BoltzGen):**
- No starting sequence; designs generated from scratch
- May bind different epitopes than OKT3

**Optimization:**
- Starts from known binders
- Preserves epitope specificity
- Generates variants with mutations

**Verify:**
- These tracks are clearly separated in the code
- Outputs are annotated with their source track
- Both tracks contribute to final candidate pool

### 5. Scientific Assumptions to Challenge

**From Assumption 1 (BoltzGen):**
> "BoltzGen has not been specifically validated on CD3 or T-cell surface proteins"
> "De novo designs may have very low hit rates (<1% functional binders)"

**From Assumption 3 (50 nM Kd):**
> "Optimal Kd may be context-dependent (tumor type, tumor antigen density)"
> "This pipeline cannot enforce 50 nM Kd computationally"

**Questions:**
- Does the code acknowledge these limitations?
- Is there a fallback if de novo designs fail calibration?
- Are affinity assumptions documented in output metadata?

### 6. Target Structure Usage

Two PDB structures are used:
- **1XIW**: CD3 epsilon-delta heterodimer
- **1SY6**: CD3 epsilon-gamma with OKT3 Fab bound

**Verify:**
- Correct chain is extracted from each structure (CD3 epsilon)
- Both structures are used for design diversity
- The OKT3 epitope from 1SY6 can be used for epitope annotation

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
