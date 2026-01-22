# CD3 Binder Design - Comprehensive Review Checklist

**Status**: In Progress
**Goal**: Zero code and scientific issues detected

## Phase 1: Discovery & Exploration

### Module Discovery
- [ ] Review `src/analysis/` modules (developability, humanness, liabilities, numbering)
- [ ] Review `src/design/` modules (affinity_variants, optimization)
- [ ] Review `src/formatting/` modules (bispecifics, crossmab, knob_hole, linkers)
- [ ] Review `src/pipeline/` modules (config, design_pipeline, filter_cascade)
- [ ] Review `src/structure/` modules (boltz_complex, interface_analysis, pdb_utils)
- [ ] Review `src/utils/` modules (constants)

### Pipeline Scripts Discovery
- [ ] Review `scripts/00_run_calibration.py`
- [ ] Review `scripts/01_setup_targets.py`
- [ ] Review `scripts/02_generate_variants.py`
- [ ] Review `scripts/03_combine_candidates.py`
- [ ] Review `scripts/04_predict_structures.py`
- [ ] Review `scripts/05_analyze_binding.py`
- [ ] Review `scripts/06_format_bispecifics.py`
- [ ] Review `scripts/07_generate_report.py`

### Test Coverage Discovery
- [ ] Review existing tests in `tests/`
- [ ] Identify gaps in test coverage
- [ ] Document critical paths without tests

### Configuration Discovery
- [ ] Review `config.yaml` for scientific parameters
- [ ] Verify thresholds have documented rationale
- [ ] Check for hardcoded magic numbers

## Phase 2: Deep Investigation

### Code Quality Deep Dive
- [ ] Trace error handling patterns across modules
- [ ] Check input validation at API boundaries
- [ ] Verify edge case handling (empty inputs, zero values, missing data)
- [ ] Check for code duplication opportunities
- [ ] Review logging and debugging support
- [ ] Verify file I/O error handling

### Scientific Accuracy Deep Dive
- [ ] **Binding affinity**: Validate pDockQ thresholds against literature
- [ ] **Interface analysis**: Verify contact counting algorithms
- [ ] **Distance calculations**: Check coordinate math and cutoffs
- [ ] **Epitope mapping**: Verify sequence alignment and residue matching
- [ ] **Humanness scoring**: Validate scoring methodology
- [ ] **Liability detection**: Verify motif patterns and CDR specificity
- [ ] **Aggregation prediction**: Validate aromatic percentage thresholds
- [ ] **Format conversions**: Verify bispecific structure correctness
- [ ] **Sequence parsing**: Check CDR extraction and numbering schemes
- [ ] **Filtering cascade**: Verify threshold relaxation logic

### Algorithm Correctness
- [ ] Verify Boltz-2 structure prediction usage
- [ ] Check sequence alignment algorithms
- [ ] Validate coordinate transformations
- [ ] Verify scoring aggregation logic
- [ ] Check statistical calculations

## Phase 3: Issue Fixing & Testing

### Critical Issues (Fix First)
- [ ] (Issues will be added as discovered)

### Major Issues
- [ ] (Issues will be added as discovered)

### Minor Issues
- [ ] (Issues will be added as discovered)

### Test Coverage Additions
- [ ] (Test gaps will be added as discovered)

## Phase 4: Final Validation

### Pre-Exit Checklist
- [ ] Re-review all modified modules
- [ ] Run full test suite (`pytest tests/ -v`)
- [ ] Verify all TODOs and FIXMEs resolved
- [ ] Check all scientific assumptions documented
- [ ] Confirm zero open issues
- [ ] Verify all commits have clear messages

## Notes & Findings

(Add notes here as you discover issues)

---

**Instructions**:
- Mark items `[x]` as you complete them
- Add new issues to Phase 3 as you discover them
- Document severity: CRITICAL, MAJOR, MINOR
- Link issues to modules: `[CRITICAL] src/structure/interface_analysis.py:123 - incorrect distance calculation`
