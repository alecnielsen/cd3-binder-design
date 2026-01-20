# Analysis Module - Scientific Context

This module handles sequence analysis for CD3 binder candidates:
- Sequence liability detection (deamidation, isomerization, glycosylation, oxidation)
- CDR numbering and identification (via ANARCI)
- Humanness scoring (via BioPhi/OASis)
- Developability assessment (charge, hydrophobicity, aggregation propensity)

## Filtering Cascade (Analysis Filters)

The analysis module implements Filters 2-5 of the cascade:

**Filter 2: HUMANNESS**
- BioPhi OASis score > 0.8 (for VH/VL)
- Sapiens humanization suggestions applied
- If humanization significantly changes CDRs, FLAG for back-mutation testing

**Filter 3: SEQUENCE LIABILITIES**
- Deamidation: No NG, NS, NT, ND, NH motifs in CDRs
- Isomerization: No DG, DS, DT, DD, DH, DN motifs in CDRs
- N-glycosylation: No N-X-S/T motifs in CDRs (where X != P)
- Oxidation: Flag exposed Met/Trp in CDRs (soft filter)
- Unpaired cysteines: 0 (must be even count)

**Filter 4: DEVELOPABILITY PROPERTIES**
- CDR-H3 length: 8-20 residues (typical range)
- Net charge: -2 to +4 (avoid extremes)
- Isoelectric point: 6.0-9.0
- Hydrophobic patches: <= 2 stretches of 5+ hydrophobic residues
- CDR hydrophobicity: Mean Kyte-Doolittle < 1.0

**Filter 5: AGGREGATION PROPENSITY**
- Low aromatic content in CDRs (< 20%)
- No consecutive aromatic residues (FF, WW, YY, FW, etc.)

## Sequence Liability Detection

```python
# Liability motif definitions
DEAMIDATION_MOTIFS = ['NG', 'NS', 'NT', 'ND', 'NH']  # N followed by small/flexible residue
ISOMERIZATION_MOTIFS = ['DG', 'DS', 'DT', 'DD', 'DH', 'DN']  # D followed by small/flexible residue
OXIDATION_RESIDUES = ['M', 'W']  # Methionine, Tryptophan

# N-glycosylation: N-X-S/T where X != P
# Pattern: N[^P][ST] in regex
```

The glycosylation motif check must verify the middle residue is NOT proline, as proline disrupts the beta-turn required for N-linked glycosylation.

## Humanness Scoring

BioPhi OASis scores humanness by comparing sequences to the human antibody repertoire (OAS database):
- Chain type must be specified correctly: 'H' for VH/VHH, 'L' for VL
- VHH scoring uses chain_type='H' (single domain, no light chain)
- For VH/VL pairs, the mean OASis score is reported
- Higher OASis score = more human-like

## Critical Assumptions

### Assumption 5: Developability Filters Are Predictive

**Claim**: Sequences passing developability filters (liabilities, charge, hydrophobicity) are more likely to succeed experimentally.

**Evidence supporting this**:
- Literature supports correlation between sequence features and developability
- Deamidation, glycosylation, and aggregation motifs are well-characterized risks

**Evidence against / Unknowns**:
- Filters are probabilistic, not deterministic
- Some "problematic" sequences express well; some "clean" sequences fail
- CDR context matters (buried vs. exposed positions)

**What would falsify this assumption**:
- No correlation between filter scores and experimental success
- High-scoring candidates fail at similar rates to low-scoring candidates

### Known Limitation: Developability Is Probabilistic

Computational filters reduce risk but don't guarantee success:
- Sequences passing all filters may still aggregate, express poorly, or be unstable
- Sequences failing soft filters may still be developable
- Context matters: CDR liabilities in buried positions may be acceptable

**Implication**: Plan for 50-70% attrition in experimental validation, even for computationally "clean" candidates.

## Implementation Notes

- CDR-specific liability filtering: `CandidateScore` has `cdr_*_count` fields; `allow_deamidation_cdr=False` only rejects CDR liabilities
- Hard vs soft filters: Deamidation/isomerization/glycosylation in CDRs are hard filters; oxidation is a soft filter
- Hydrophobic residues for patch detection: A, I, L, M, F, V, W (excludes Y which has Kyte-Doolittle -1.3)
