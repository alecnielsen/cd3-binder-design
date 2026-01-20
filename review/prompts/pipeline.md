# Scientific Review: Pipeline Module

You are reviewing the `src/pipeline/` module of a CD3 binder design pipeline. This module handles:
- Pipeline configuration and validation
- Filter cascade orchestration
- Candidate ranking and selection
- Fallback logic when insufficient candidates survive filtering
- Report generation

## Important: Module README Updates

The scientific context for this module lives in `src/pipeline/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. Filter Cascade Order and Logic

The filtering cascade should be applied in the documented order:

1. **Binding quality** (calibrated thresholds): pDockQ, interface area, contacts
2. **Humanness**: OASis score > 0.8
3. **Sequence liabilities**: Deamidation, isomerization, glycosylation, oxidation
4. **Developability**: CDR-H3 length, net charge, pI, hydrophobic patches
5. **Aggregation propensity**: Aromatic content, consecutive aromatics

**Verify:**
- Filters are applied in the correct order
- Hard filters (liabilities in CDRs) vs soft filters (oxidation) are distinguished
- Calibrated thresholds are used when available, defaults otherwise

### 2. Configuration Handling

Configuration is loaded from `config.yaml` with nested structure:

**Verify:**
- Nested keys are accessed correctly (e.g., `filtering.binding.min_pdockq`)
- Default values are provided for optional keys
- Calibrated thresholds override defaults when `use_calibrated: true`

**From CLAUDE.md:**
> "Config supports nested schema - `filtering.binding.min_pdockq` or flat `filtering.min_pdockq`"

### 3. Fallback Logic

When insufficient candidates survive filtering (< 10), fallback is triggered:

**Fallback sequence:**
1. Relax soft filters first (oxidation, hydrophobic patches)
2. Relax thresholds by up to 10% toward calibration baseline
3. Include borderline candidates with explicit risk flags
4. Document all relaxations in output

**Verify:**
- Fallback is triggered at the correct threshold (`min_candidates: 10`)
- Soft filters are relaxed before hard filters
- Relaxation is limited (`max_threshold_relaxation: 0.1`)
- Risk flags are added to borderline candidates

**From CLAUDE.md:**
> "Fallback config wired - `relax_soft_filters_first` and `max_threshold_relaxation` are now used"

### 4. Candidate Ranking

Composite score calculation with documented weights:
- Binding confidence (pDockQ): 30%
- Humanness (OASis): 25%
- Liability count (inverse): 25%
- Developability metrics: 20%

**Verify:**
- Weights sum to 100%
- Score normalization is applied consistently
- Higher scores = better candidates

**Diversity requirement:**
- Cluster by CDR-H3 similarity (70% identity threshold)
- Select top candidate from each cluster
- Ensure mix of: de novo + optimized, VHH + scFv, OKT3 + novel epitope

### 5. Calibration Integration

Calibration results should be loaded and used:

**Verify:**
- Calibrated thresholds are loaded from config
- If calibration not run, defaults are used with warning
- Calibration results include all three known binders

### 6. Epitope Annotation

De novo designs should be annotated with epitope information:

**Verify:**
- CD3 epsilon contact residues are extracted from Boltz-2 predictions
- Overlap with OKT3 epitope is calculated
- Designs with < 50% overlap are flagged as "NOVEL EPITOPE"

**From config:**
> `okt3_epitope_residues: [23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 38, 39, 40, 41, 42, 45, 47]`

### 7. Provenance and Reproducibility

Output files should include provenance metadata:

**Verify:**
- Pipeline version and git commit are recorded
- Run timestamp is included
- Config hash is computed
- Tool versions are documented

### 8. Success Criteria Checks

The pipeline should check computational success criteria:

| Checkpoint | Success Criterion |
|------------|-------------------|
| Calibration | All 3 known binders pass |
| De novo | >= 20 designs pass initial filter |
| Filtering | >= 10 candidates survive cascade |
| Diversity | Mix of VHH/scFv, de novo/optimized, OKT3-like/novel |

**Verify:**
- These checks are implemented
- Appropriate warnings or errors are raised when criteria fail
- Fallback is triggered when appropriate

### 9. Scientific Assumptions to Challenge

**Questions:**
- Does the pipeline correctly handle edge cases (zero candidates, all from one track)?
- Is the ranking formula scientifically justified, or arbitrary?
- Are there race conditions or ordering issues in filter application?
- Does the diversity requirement actually ensure diverse final candidates?

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
