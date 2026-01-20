# Formatting Module - Scientific Context

This module handles assembly of bispecific antibody constructs in five formats:
1. CrossMab (Fab x Fab with CH1-CL domain swap)
2. Asymmetric Fab + scFv (knob-in-hole)
3. Asymmetric Fab + VHH (knob-in-hole)
4. IgG-scFv Morrison (symmetric, 2x scFv at C-terminus)
5. IgG-VHH Morrison (symmetric, 2x VHH at C-terminus)

## Format Conversion Pipeline Step

Convert each candidate into multiple bispecific formats:
- CrossMab (Fab x Fab)
- Asymmetric Fab + scFv
- Asymmetric Fab + VHH
- IgG-scFv Morrison (2x scFv, symmetric)
- IgG-VHH Morrison (2x VHH, symmetric)

Include:
- Human IgG1 Fc with knob-in-hole mutations (for asymmetric)
- Standard linkers (G4S)x3 for scFv, (G4S)x2 for Fc-VHH fusion
- Placeholder tumor target arm (anti-HER2 trastuzumab, configurable)

## Bispecific Antibody Formats

### Format 1: CrossMab (Fab x Fab)

**Components:**
- Heavy chain 1: VH(HER2) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2: VH(CD3) - CL - Hinge - CH2 - CH3(Hole) [CrossMab swap]
- Light chain 1: VL(HER2) - CL
- Light chain 2: VL(CD3) - CH1 [CrossMab swap]

**CrossMab CH1-CL swap:**
On the CD3 arm, CH1 and CL domains are swapped between heavy and light chains to ensure correct light chain pairing. This prevents light chain mispairing in the bispecific.

### Format 2: Asymmetric Knob-in-Hole (Fab + scFv)

**Components:**
- Heavy chain 1: VH(HER2) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2: scFv(CD3) - Hinge - CH2 - CH3(Hole)
- Light chain: VL(HER2) - CL (only one light chain needed)

**scFv linker:** (GGGGS)x3 = GGGGSGGGGSGGGGS (15 aa)

### Format 3: Asymmetric Knob-in-Hole (Fab + VHH)

**Components:**
- Heavy chain 1: VH(HER2) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2: VHH(CD3) - Hinge - CH2 - CH3(Hole)
- Light chain: VL(HER2) - CL

### Format 4: IgG-scFv Morrison (Symmetric, 2x CD3)

**Components:**
- Heavy chain: VH(HER2) - CH1 - Hinge - CH2 - CH3 - Linker - scFv(CD3)
- Light chain: VL(HER2) - CL

**C-terminal fusion linker:** (GGGGS)x2 = GGGGSGGGGS (10 aa)

**Valency:**
- 2x HER2 binding (bivalent)
- 2x CD3 binding (bivalent)

### Format 5: IgG-VHH Morrison (Symmetric, 2x CD3)

**Components:**
- Heavy chain: VH(HER2) - CH1 - Hinge - CH2 - CH3 - Linker - VHH(CD3)
- Light chain: VL(HER2) - CL

## Knob-in-Hole Mutations (EU Numbering)

For asymmetric formats, Fc engineering ensures correct heavy chain pairing:

**Knob (on one heavy chain):**
- T366W (threonine to tryptophan at position 366)

**Hole (on the other heavy chain):**
- T366S (threonine to serine)
- L368A (leucine to alanine)
- Y407V (tyrosine to valine)

**Important:** For Morrison format (symmetric), knob-in-hole should NOT be used since both heavy chains are identical.

## Linker Sequences

```yaml
linkers:
  scfv: "GGGGSGGGGSGGGGS"       # (G4S)x3, 15 aa
  fc_fusion: "GGGGSGGGGS"       # (G4S)x2, 10 aa
```

## CrossMab Assembly Logic

The actual implementation uses a `CrossMabFormatter` class that:
1. Retrieves constant region sequences (CH1, CL, hinge, CH2, CH3) from a `SequenceLibrary`
2. Applies knob (T366W) and hole (T366S/L368A/Y407V) mutations to CH3 domains
3. Assembles chains with the CrossMab domain swap on the CD3 arm

```python
# Simplified assembly logic (see crossmab.py for full implementation):

# Target arm: standard Fab with Knob
heavy_chain_1 = target_vh + ch1 + hinge + ch2 + ch3_knob
light_chain_1 = target_vl + cl

# CD3 arm: CrossMab swap (CH1-CL exchanged) with Hole
heavy_chain_2 = cd3_vh + cl + hinge + ch2 + ch3_hole  # VH fused to CL (swap)
light_chain_2 = cd3_vl + ch1                           # VL fused to CH1 (swap)
```

## Implementation Notes

- scFv is constructed as VH-linker-VL (not VL-linker-VH)
- VHH formats have no scFv linker (single domain)
- Light chain count varies: 2 for CrossMab, 1 for Fab+scFv/VHH, 2 for Morrison
- Asymmetric formats have different heavy chains; symmetric formats have identical heavy chains
- Constant regions should be IgG1 subclass
