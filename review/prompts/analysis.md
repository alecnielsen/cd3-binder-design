# Scientific Review: Analysis Module

You are reviewing the `src/analysis/` module of a CD3 binder design pipeline. This module handles:
- Sequence liability detection (deamidation, isomerization, glycosylation, oxidation)
- CDR numbering and identification (via ANARCI)
- Humanness scoring (via BioPhi/OASis)
- Developability assessment (charge, hydrophobicity, aggregation propensity)

## Important: Module README Updates

The scientific context for this module lives in `src/analysis/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. Liability Motif Correctness
Verify the liability detection patterns match established biochemistry:

**Deamidation motifs** - Should detect N followed by small residues:
- NG, NS, NT, ND (asparagine followed by glycine, serine, threonine, aspartate)
- NG is the most labile; rate depends on sequence context

**Isomerization motifs** - Should detect D followed by small residues:
- DG, DS, DT, DD (aspartate followed by small residues)
- DG is the most labile

**N-glycosylation motifs** - Should detect N-X-S/T where X is NOT proline:
- Pattern: N[^P][ST] in regex
- Proline disrupts the beta-turn required for glycosylation

**Oxidation** - Should flag:
- Methionine (M) - oxidizes to methionine sulfoxide
- Tryptophan (W) - can oxidize, especially when solvent-exposed

**Questions to verify:**
- Are the motif patterns complete? Any missing high-risk motifs?
- Does the code correctly handle edge cases (sequence boundaries)?
- Is CDR-specific filtering implemented correctly (liabilities in CDRs vs frameworks)?

### 2. CDR Numbering Accuracy
Verify the ANARCI wrapper correctly:
- Supports IMGT, Chothia, and Kabat numbering schemes
- Handles VH, VL, and VHH chains appropriately
- Correctly identifies CDR boundaries for each scheme:
  - IMGT: CDR1 (27-38), CDR2 (56-65), CDR3 (105-117)
  - Kabat: CDR-H1 (31-35B), CDR-H2 (50-65), CDR-H3 (95-102)
- Returns consistent data structures

### 3. Humanness Scoring
Verify the BioPhi/OASis integration:
- Correctly passes chain_type ('H' for VH/VHH, 'L' for VL)
- Handles VHH scoring (single domain, no light chain)
- Aggregates scores appropriately for VH/VL pairs
- Score interpretation: higher OASis = more human-like

### 4. Developability Calculations
Verify physical property calculations:
- Net charge calculation at pH 7.4 (pKa values for D, E, H, K, R)
- Isoelectric point estimation
- Hydrophobic patch detection (consecutive hydrophobic residues)
- CDR-H3 length extraction

**Hydrophobic residues for patches**: A, V, I, L, M, F, Y, W

### 5. Scientific Assumptions to Challenge

From Assumption 5 in the README:
> "Sequences passing developability filters (liabilities, charge, hydrophobicity) are more likely to succeed experimentally."

**Questions:**
- Are the thresholds (e.g., max 2 hydrophobic patches, charge -2 to +4) scientifically justified?
- Does the code correctly distinguish CDR vs framework liabilities?
- Is there appropriate handling of "soft" vs "hard" filters?

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
