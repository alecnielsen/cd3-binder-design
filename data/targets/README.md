# CD3 Target Structures

This directory contains the CD3 structures used as design targets.

## Structures

### 1XIW - CD3εδ Heterodimer
- **Resolution**: 1.9 Å
- **Description**: Human CD3 epsilon and delta extracellular domains
- **Use**: Primary target for de novo design
- **Download**: https://www.rcsb.org/structure/1XIW

### 1SY6 - CD3εγ + OKT3 Fab Complex
- **Resolution**: 2.8 Å
- **Description**: Human CD3 epsilon/gamma with OKT3 Fab bound
- **Use**: Reference for validated binding mode; secondary design target
- **Download**: https://www.rcsb.org/structure/1SY6

## Setup

Run the setup script to download structures:

```bash
python scripts/01_setup_targets.py
```

Or download manually:

```bash
# Download from RCSB
curl -o cd3_epsilon_delta_1XIW.pdb "https://files.rcsb.org/download/1XIW.pdb"
curl -o cd3_epsilon_gamma_1SY6.pdb "https://files.rcsb.org/download/1SY6.pdb"
```

## Chain Information

### 1XIW
1XIW contains two copies of the CD3εδ heterodimer plus UCHT1 scFv:
- **Chains A, E**: CD3ε (epsilon) - 104 residues each
- **Chains B, F**: CD3δ (delta) - 79 residues each
- **Chains C, G**: UCHT1 VL (light chain variable region)
- **Chains D, H**: UCHT1 VH (heavy chain variable region)

Use chain A or E for canonical CD3ε sequence (they are identical).

### 1SY6
**Note:** Chain IDs in the actual PDB file differ from some documentation.
- Chain A: CD3ε/γ fusion protein - includes CD3ε extracellular domain with linker
  (168 residues, numbered 0-203). This is NOT pure CD3ε - it's a γ/ε fusion construct.
- Chain H: OKT3 VH (heavy chain variable region)
- Chain L: OKT3 VL (light chain variable region)

**Important:** Residue numbering in 1SY6 chain A differs from 1XIW chains A/E due to
the fusion construct. When comparing epitopes across structures, use sequence
alignment (see `compare_epitopes_aligned()` in InterfaceAnalyzer).

The OKT3 epitope on CD3ε (chain A) spans residues ~139-192 in PDB numbering,
primarily contacting the FG loop region.

## Notes

- For BoltzGen design, extract the CD3ε chain as the target
- The OKT3 complex (1SY6) shows the validated binding epitope
- Both structures are human CD3; cynomolgus cross-reactivity requires SP34-based designs
