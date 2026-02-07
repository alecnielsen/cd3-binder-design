# Pipeline Module - Scientific Context

This module handles:
- Pipeline configuration and validation
- Filter cascade orchestration
- Candidate ranking and selection
- Fallback logic when insufficient candidates survive filtering
- Report generation

## Computational Pipeline Overview

The pipeline follows this flow:

1. **Input Preparation**: Download CD3 structures, prepare starting sequences
2. **Design Generation**: De novo (BoltzGen) + Optimization tracks
3. **CDR Identification**: ANARCI numbering (if ANARCI is unavailable, CDR positions are omitted and CDR-specific filters fall back to sequence-level heuristics)
4. **Structure Prediction**: ABodyBuilder2 + Boltz-2
5. **Calibration Phase**: Set thresholds using known binders
6. **Filtering Cascade**: Binding quality -> Humanness -> Liabilities -> Developability -> Aggregation
7. **Epitope Annotation**: Compare to OKT3 epitope
8. **Candidate Ranking**: Worst-metric-rank + greedy maximin diversity selection
9. **Format Conversion**: Generate bispecific sequences
10. **Final Output**: ~10 candidates with sequences, CIF structures, and scorecards

## Filter Cascade Order

Filters are applied in the documented order:

1. **Binding quality** (calibrated thresholds): interface area, contacts (pDockQ only if non-zero)
2. **Humanness**: OASis score > 0.8
3. **Sequence liabilities**: Deamidation, isomerization, glycosylation, oxidation
4. **Developability**: CDR-H3 length, net charge, pI, hydrophobic patches
5. **Aggregation propensity**: Aromatic content, consecutive aromatics

**Important:** Hard filters (liabilities in CDRs) vs soft filters (oxidation) are distinguished.
**Humanness scoring:** `None` scores (e.g., BioPhi missing) are treated as soft-fail; a real score of `0.0` is treated as a hard fail below threshold.

## Candidate Ranking

**Default: Worst-metric-rank** (replaces legacy composite score)

Each candidate is ranked independently across 6 metrics (ipTM, pTM, pLDDT, interface area, contacts, humanness). Each rank is divided by the metric's importance weight. The quality key is the maximum weighted rank — i.e., the worst (highest) rank across all metrics determines overall quality. Sort ascending (lower = better). Tiebreak by ipTM.

Default metric weights:
| Metric | Weight | Rationale |
|--------|--------|-----------|
| ipTM | 1 (high) | Best interface quality signal |
| pTM | 1 (high) | Fold confidence |
| interface_area | 1 (high) | Binding interface size |
| humanness | 1 (high) | Critical for therapeutics |
| num_contacts | 2 (lower) | Correlated with area |
| pLDDT | 2 (lower) | Typically high, less discriminating |

**Diversity selection: Greedy maximin** (alpha=0.001)

After ranking, iteratively picks candidates maximizing:
`(1 - alpha) * quality + alpha * (1 - max_identity_to_selected)`

This prevents near-duplicate sequences from occupying multiple slots. At alpha=0.001, quality dominates but identical sequences are penalized.

**Legacy: Composite score** (`method: composite` in config)

Available for backwards compatibility. Uses 30% pDockQ + 25% humanness + 25% liabilities + 20% developability. Not recommended because pDockQ is always 0.0 from Boltz-2.

## Fallback Logic

When insufficient candidates survive filtering (< 10), fallback is triggered:

1. Relax soft filters first (oxidation, hydrophobic patches)
2. Relax thresholds by up to 10% toward calibration baseline
3. Include borderline candidates with explicit risk flags
4. Document all relaxations in output

## Configuration

```yaml
ranking:
  method: worst_metric_rank       # or "composite" for legacy behavior
  diversity_alpha: 0.001          # Greedy maximin diversity weight
  use_diversity_selection: true

filtering:
  binding:
    min_pdockq: 0.5               # Only applied if non-zero; Boltz-2 always returns 0.0
    min_interface_area: 800       # A^2, override with calibrated
    min_contacts: 10              # Override with calibrated
    use_calibrated: true          # Use calibrated thresholds if available

  humanness:
    min_oasis_score: 0.8
    generate_back_mutations: true

  liabilities:
    allow_deamidation_cdr: false
    allow_isomerization_cdr: false
    allow_glycosylation_cdr: false
    max_oxidation_sites: 2        # Soft filter

  developability:
    cdr_h3_length_range: [8, 20]
    net_charge_range: [-2, 4]
    pi_range: [6.0, 9.0]
    max_hydrophobic_patches: 2

  fallback:
    min_candidates: 10            # Trigger fallback if fewer survive
    relax_soft_filters_first: true
    max_threshold_relaxation: 0.1 # 10% toward calibration baseline

output:
  export_cif: true                # Save CIF structure files to structures/cif/

epitope_annotation:
  okt3_epitope_residues: null     # null = dynamic extraction from 1SY6; or provide explicit list
  overlap_threshold: 0.5          # Flag as "novel epitope" if < 50% overlap
```

## Success Criteria

| Checkpoint | Success Criterion | Action if Fail |
|------------|-------------------|----------------|
| Calibration | All 3 known binders pass interface area/contacts thresholds | Investigate Boltz-2 predictions |
| De novo design | >= 20 designs pass binding filter (interface area, contacts) | Increase design count or relax threshold |
| Filtering | >= 10 candidates survive full filter cascade | Apply fallback |
| Diversity | Mix of VHH/scFv, de novo/optimized, OKT3-like/novel | Manually adjust selection |
| Epitope coverage | >= 3 candidates have OKT3-like epitope | Prioritize optimization track |

## Failure Mode Detection

| Failure Mode | How to Detect | Mitigation |
|--------------|---------------|------------|
| BoltzGen doesn't generate CD3 binders | All designs fail binding assays | Use optimization track only |
| Boltz-2 gives false positives | High pDockQ designs don't bind | Add orthogonal filtering |
| All designs bind same (wrong) epitope | Epitope binning shows single cluster | Redesign against different CD3 epsilon face |
| Humanization destroys all binders | All humanized variants lose binding | Use back-mutation variants |
| Affinity too high (CRS risk) | All designs have Kd < 10 nM | Generate affinity-attenuated variants |
| Affinity too low (no killing) | All designs have Kd > 500 nM | Affinity maturation |

## Provenance and Reproducibility

Output files include provenance metadata:
- Pipeline version and git commit
- Run timestamp

```yaml
_provenance:
  pipeline_version: "1.0.0"
  git_commit: "abc1234"
  run_timestamp: "2026-01-15T10:30:00Z"
```

**Note:** Config hash is available via `PipelineConfig.config_hash()` but not automatically included in provenance. Tool versions are not currently tracked.

## Implementation Notes

- Config supports nested schema: `filtering.binding.min_pdockq` or flat `filtering.min_pdockq`
- Fallback config is wired: `relax_soft_filters_first` and `max_threshold_relaxation` are used
- Epitope residues default to dynamic extraction from 1SY6; `config.epitope.okt3_epitope_residues` overrides with explicit list
- Aggregation filter is active: CDR-specific: >20% aromatic or 2+ consecutive aromatics. Fallback (no CDRs): >15% aromatic or 3+ consecutive aromatics.
- scFv inputs are parsed into VH/VL when only a full scFv sequence is provided, so all downstream analyses operate on the correct chain boundaries.
