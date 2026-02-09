# Analysis Module - Scientific Context

This module handles sequence analysis for CD3 binder candidates:
- Sequence liability detection (deamidation, isomerization, glycosylation, oxidation)
- CDR numbering and identification (via ANARCI)
- Humanness scoring (via BioPhi/OASis)
- Developability assessment (charge, hydrophobicity, aggregation propensity)
- Affinity proxy scoring (ProteinMPNN + AntiFold log-likelihoods)

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

**Filter 5: AGGREGATION PROPENSITY (soft metrics)**
The following are tracked as soft filters (flag but don't reject):
- CDR hydrophobicity: Mean Kyte-Doolittle (tracked)
- Aromatic content in CDRs (tracked, threshold 20% in config)
- Consecutive aromatic residues (FF, WW, YY, FW, etc. - counted as aromatic clusters)

These feed into the worst-metric-rank ranking system rather than directly penalizing a composite score.

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
- A score of `0.0` is treated as a real value (i.e., it fails the threshold, not a soft-fail)

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

## CDR Numbering Schemes

The module supports IMGT, Chothia, and Kabat numbering via ANARCI. Key boundaries:
If ANARCI is not available or numbering fails, CDR positions are omitted (no positions are assigned).

**IMGT** (default):
- CDR1: 27-38, CDR2: 56-65, CDR3: 105-117 (same for H and L chains)

**Chothia** (extended/AbM definitions):
- H1: 26-32, H2: 52-58, H3: 95-102
- L1: 24-34, L2: 50-56, L3: 89-97

Note: The original structural Chothia definition uses H2: 52-56, but we use the "extended Chothia" (AbM) definition H2: 52-58 which is more common in antibody engineering tools.

**Kabat**:
- H1: 31-35, H2: 50-65, H3: 95-102
- L1: 24-34, L2: 50-56, L3: 89-97

## Implementation Notes

- CDR-specific liability filtering: `LiabilityReport` has `cdr_liabilities` (hard: deamidation, isomerization, glycosylation) and `cdr_oxidation_count` (soft) fields
- Hard vs soft filters: Deamidation/isomerization/glycosylation in CDRs are hard filters (cause rejection); oxidation in CDRs is a soft filter (flags but doesn't reject)
- Hydrophobic residues for patch detection: A, I, L, M, F, V, W (excludes Y which has Kyte-Doolittle -1.3)
- BioPhi soft-fail: If BioPhi is not installed, humanness scoring returns `None` scores and the pipeline continues (soft-fail behavior). A real score of `0.0` is treated as a valid value (fails the threshold).

## Affinity Proxy Scoring (Step 04a)

`affinity_scoring.py` provides inverse folding log-likelihood scoring as an affinity proxy for de novo designs. Scores are computed in step 04a (pre-filter + scoring) and **used as ranking metrics** in step 05 (proteinmpnn_ll weight=1, antifold_ll weight=2, protenix_iptm weight=1).

### Tools

| Tool | License | Metric | Validation |
|------|---------|--------|------------|
| **ProteinMPNN** | MIT | Log-likelihood of binder sequence given complex structure | Spearman r=0.27-0.41 on AbBiBench (best for de novo) |
| **AntiFold** | BSD-3 | Antibody-specific log-likelihood | Supports nanobodies via `--nanobody_chain` |

### Usage

```python
from src.analysis.affinity_scoring import batch_score_affinity

results = batch_score_affinity(
    cif_dir="data/outputs/structures/cif/",
    design_ids=["fab_1XIW_0040", "vhh_1XIW_0008"],
    binder_types=["scfv", "vhh"],
)
for r in results:
    print(f"{r.design_id}: MPNN={r.proteinmpnn_ll}, AF={r.antifold_ll}")
```

### Implementation Details

- **ProteinMPNN** uses `parse_PDB` → `tied_featurize` → model forward pass → `_scores`. Since `parse_PDB` only handles PDB format, CIF files are converted via `_cif_to_pdb()` (BioPython `MMCIFParser` + `PDBIO`). Model weights: `v_48_020.pt` (vanilla, 48 edges, 0.2Å noise). Cached globally to avoid reloading per candidate.
- **AntiFold** uses `load_model()` → `get_pdbs_logits()` with a DataFrame specifying chain mappings. For VHH/scFv (single binder chain B), requires `custom_chain_mode=True` and `Lchain=None`. For 3-chain Fab files, uses standard mode with `Hchain="B"`, `Lchain="C"`. NLL computed from `pdb_res` column and amino acid logit columns.

### Important Caveats

- Log-likelihoods are NOT direct affinity predictors
- Higher values indicate better structural complementarity, which correlates with binding
- ProteinMPNN operates on the Boltz-2 predicted structure, so errors in prediction propagate
- Both tools run locally on CPU; no Modal/GPU required
- If a tool is not installed, `batch_score_affinity()` returns an error message (not an exception)
- `pip install proteinmpnn antifold` downgrades torch to 2.3.1 — no impact on CPU scoring

### ANTIPASTI (Deprecated)

`antipasti_scoring.py` is deprecated and redirects to `affinity_scoring.py`. ANTIPASTI was not integrated because:
- Requires R dependency (rpy2), adding environment complexity
- No VHH/nanobody support
- Performance degrades on predicted (vs experimental) structures
