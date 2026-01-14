# Starting Sequences for CD3 Binder Optimization

This directory contains sequences for known anti-CD3 antibodies that serve as starting points for the optimization track.

## Files

### teplizumab.yaml
Humanized OKT3 (hOKT3γ1 Ala-Ala) - FDA approved, well-characterized.
- High affinity (~5 nM)
- Does NOT cross-react with cynomolgus
- Use for human-only studies or as CDR source for humanization

### sp34.yaml
Murine SP34 with humanization suggestions.
- Moderate affinity (~10 nM)
- DOES cross-react with cynomolgus (critical for preclinical)
- Requires humanization before therapeutic use

### ucht1.yaml
Murine UCHT1 with humanization suggestions.
- Lower affinity (~20 nM)
- Alternative epitope to OKT3
- Does NOT cross-react with cynomolgus

### affinity_mutations.yaml
Literature-derived mutations for reducing CD3 affinity.
- Mutations to reduce Kd by ~10x, ~50x, ~100x
- Targeting ~50-200 nM range for bispecific use
- Based on OKT3/teplizumab hot spots

## Selection Guide

| Use Case | Recommended Starting Sequence |
|----------|------------------------------|
| Human clinical + cyno preclinical | SP34 (humanized) |
| Human clinical only | Teplizumab |
| Alternative epitope exploration | UCHT1 (humanized) |
| Affinity attenuation study | Teplizumab + affinity_mutations |

## Sequence Format

All YAML files contain:
- `vh_sequence` / `vl_sequence`: Full variable region sequences
- `cdrs`: CDR sequences extracted (IMGT numbering)
- `frameworks`: Framework region sequences
- `binding`: Known binding properties
- `humanization`: Suggested human frameworks (for murine sequences)

## Notes

1. **Patent Status**: All sequences are from expired patents or publications >20 years old
2. **Humanization**: Murine sequences (SP34, UCHT1) require humanization; suggestions provided but verify humanness scores
3. **Cynomolgus Cross-Reactivity**: Only SP34 cross-reacts with cynomolgus CD3
4. **Affinity Targets**: For bispecific use, ~50 nM Kd is often preferred over high affinity
