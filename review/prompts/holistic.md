# Holistic Scientific Review: CD3 Binder Design Pipeline

You are performing a holistic review of the entire CD3 binder design codebase. This is an **agent-style review** - explore the codebase using Read, Glob, and Grep tools rather than having all code provided upfront.

## Repository Structure

```
cd3-binder-design/
├── src/
│   ├── analysis/      # Liability detection, CDR numbering, humanness, developability
│   ├── design/        # BoltzGen de novo design, optimization
│   ├── formatting/    # Bispecific antibody assembly (5 formats)
│   ├── pipeline/      # Filter cascade, ranking, fallback logic
│   ├── structure/     # Interface analysis, PDB utilities
│   └── utils/         # Constants, helpers
├── scripts/           # Pipeline orchestration (00-07)
├── modal/             # GPU deployments (BoltzGen, Boltz-2)
├── data/              # Targets, starting sequences, frameworks
└── config.yaml        # Pipeline configuration
```

## Review Strategy

1. **Read the history first** - Check what previous iterations found and fixed
2. **Explore systematically** - Look at areas not yet reviewed
3. **Cross-check** - Verify constants/thresholds are consistent across modules
4. **Fix issues directly** - Edit files to fix problems you find
5. **Track your work** - Report findings in structured format

## Scientific Review Criteria

### Analysis Module (`src/analysis/`)
- **Liabilities**: Deamidation (NG, NS, NT, ND, NH), isomerization (DG, DS, DT, DD, DH, DN), glycosylation (N[^P][ST]), oxidation (M, W)
- **CDR numbering**: IMGT boundaries (27-38, 56-65, 105-117), Chothia/Kabat support
- **Humanness**: BioPhi/OASis scoring, VHH handling (no light chain)
- **Developability**: Net charge at pH 7.4, hydrophobic patches, CDR-H3 length

### Design Module (`src/design/`)
- **BoltzGen**: VHH (~120 aa) vs scFv (~250 aa), seed reproducibility
- **Affinity attenuation**: ~50 nM target, literature-derived mutations
- **Starting sequences**: Teplizumab (humanized), SP34 (cyno cross-reactive), UCHT1

### Formatting Module (`src/formatting/`)
- **CrossMab**: CH1-CL swap on CD3 arm only
- **Knob-in-hole**: T366W (knob), T366S/L368A/Y407V (hole) - EU numbering
- **Linkers**: scFv = (GGGGS)×3, Fc fusion = (GGGGS)×2
- **Morrison**: Symmetric (homodimeric), no knob-in-hole

### Structure Module (`src/structure/`, `modal/`)
- **pDockQ**: Structural confidence, NOT affinity predictor
- **Contacts**: Residue-level pairs within 5Å cutoff
- **Calibration**: All 3 known binders must pass; thresholds = min - margin
- **Target extraction**: Must extract only CD3ε chain from multi-chain PDBs

### Pipeline Module (`src/pipeline/`)
- **Filter order**: Binding → Humanness → Liabilities → Developability → Aggregation
- **Soft vs hard filters**: Oxidation soft, deamidation in CDRs hard
- **Fallback**: Relax soft filters first, max 10% threshold relaxation
- **Diversity**: Cluster by CDR-H3 (70% identity), mix of VHH/scFv, de novo/optimized

### Scripts (`scripts/`)
- **Dependencies**: Calibration before filtering, targets before design
- **Reproducibility**: Seeds passed through, provenance metadata
- **Error handling**: Clear messages, non-zero exit on failure

### Cross-Cutting Concerns
- **Constants consistency**: Same values in constants.py, config.yaml, and code
- **PDB numbering**: 1XIW starts at residue 12, not 1 - handle via alignment
- **1SY6 fusion**: Contains CD3γ/ε fusion, not pure CD3ε
- **Licensing**: All tools must be MIT/BSD/Apache/CC-BY (no IgFold, no NetMHCIIpan)

## After Making Fixes

**IMPORTANT**: After fixing any issues, you MUST commit your changes before reporting:

```bash
git add -A && git commit -m "Ralph review: <brief description of fixes>"
```

This ensures each iteration's fixes are tracked and the next iteration starts from a clean state.

## Output Format

### If you find issues:

```markdown
## Files Examined
- path/to/file1.py (lines X-Y)
- path/to/file2.py (full file)

## Issues Found and Fixed
### Issue 1: [Brief description]
- **File**: path/to/file.py:LINE
- **Problem**: What was wrong
- **Fix**: What you changed

### Issue 2: [Brief description]
...

## Areas Still To Review
- [ ] src/module/ - not yet examined
- [ ] specific concern not yet checked
```

### If no issues found:

After thoroughly reviewing the codebase and verifying previous fixes:

```
NO_ISSUES
```

**Rules:**
- Only output `NO_ISSUES` if you made zero code edits this iteration
- README updates do NOT count as code issues
- Do not invent issues - only report genuine scientific or implementation problems
- Verify fixes from previous iterations still work
