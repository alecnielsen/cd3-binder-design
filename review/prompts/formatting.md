# Scientific Review: Formatting Module

You are reviewing the `src/formatting/` module of a CD3 binder design pipeline. This module handles assembly of bispecific antibody constructs in five formats:
1. CrossMab (Fab x Fab with CH1-CL domain swap)
2. Asymmetric Fab + scFv (knob-in-hole)
3. Asymmetric Fab + VHH (knob-in-hole)
4. IgG-scFv Morrison (symmetric, 2x scFv at C-terminus)
5. IgG-VHH Morrison (symmetric, 2x VHH at C-terminus)

## Important: Module README Updates

The scientific context for this module lives in `src/formatting/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. CrossMab Domain Swap Correctness

The CrossMab format uses a CH1-CL domain swap on ONE arm to ensure correct light chain pairing:

**Standard Fab arm (tumor target):**
- Heavy: VH - CH1 - Hinge - CH2 - CH3
- Light: VL - CL

**CrossMab arm (CD3):**
- Heavy: VH - CL - Hinge - CH2 - CH3 (note: CL instead of CH1)
- Light: VL - CH1 (note: CH1 instead of CL)

**Verify:**
- Is the domain swap applied to the correct arm (CD3, not target)?
- Are the sequences concatenated in the correct order?
- Is the hinge sequence placed correctly after CL in the swapped heavy chain?

### 2. Knob-in-Hole Mutations (EU Numbering)

For asymmetric formats, verify correct Fc engineering mutations:

**Knob (on one heavy chain):**
- T366W (threonine to tryptophan at position 366)

**Hole (on the other heavy chain):**
- T366S (threonine to serine)
- L368A (leucine to alanine)
- Y407V (tyrosine to valine)

**Verify:**
- Are mutations applied to the correct positions (EU numbering)?
- Is the knob consistently on one arm and hole on the other?
- For Morrison format (symmetric), knob-in-hole should NOT be used

### 3. Linker Sequences

Standard linker sequences:

**scFv linker (VH to VL):** (GGGGS)x3 = "GGGGSGGGGSGGGGS" (15 aa)
- Flexible, hydrophilic
- Sufficient length for VH-VL orientation

**Fc fusion linker (CH3 to scFv/VHH):** (GGGGS)x2 = "GGGGSGGGGS" (10 aa)
- C-terminal fusion in Morrison format

**Verify:**
- Are the correct linker sequences used for each format?
- Is the scFv constructed as VH-linker-VL (not VL-linker-VH)?
- For VHH formats, no scFv linker should be present

### 4. Chain Assembly Order

Verify the correct domain order for each chain type:

**Full IgG heavy chain:** VH - CH1 - Hinge - CH2 - CH3
**Light chain:** VL - CL (kappa or lambda)
**scFv:** VH - (G4S)3 - VL
**Morrison heavy:** VH - CH1 - Hinge - CH2 - CH3 - (G4S)2 - scFv/VHH

**Verify:**
- No missing domains in the assembly
- Correct order of constant region domains (CH1, hinge, CH2, CH3)
- Hinge is included between CH1 and CH2

### 5. Valency and Format Correctness

**Asymmetric formats (CrossMab, Fab+scFv, Fab+VHH):**
- 1x tumor target binding
- 1x CD3 binding
- Heterodimeric Fc (knob-in-hole)

**Morrison formats (IgG-scFv, IgG-VHH):**
- 2x tumor target binding (both Fabs)
- 2x CD3 binding (both C-terminal fusions)
- Homodimeric Fc (no knob-in-hole)

**Verify:**
- Asymmetric formats have different heavy chains
- Symmetric formats have identical heavy chains
- Light chain count is correct (2 for CrossMab, 1 for Fab+scFv/VHH, 2 for Morrison)

### 6. Scientific Assumptions to Challenge

**Questions:**
- Does the code handle VHH inputs correctly (single domain, no light chain)?
- Are constant region sequences from the correct IgG subclass (IgG1)?
- Is there validation that input sequences are the correct chain type?

## Output Format

If you find **code issues**:
1. **FIX THEM** by editing the code files directly
2. Then report what you fixed:
   - File path and line number
   - Description of what was wrong
   - What you changed and why

If no **code issues** are found (or after fixing all issues), respond with exactly:
NO_ISSUES

**Important notes:**
- Do not invent issues if the code is correct. Only report genuine scientific or implementation problems.
- README updates do NOT count as issues. If you update the module README.md to fix documentation, still output "NO_ISSUES" if no code problems were found.
- "NO_ISSUES" means the code is scientifically correct, even if you made documentation improvements.
