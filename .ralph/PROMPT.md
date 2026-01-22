# Comprehensive Code and Scientific Review

You are conducting a **thorough, iterative review** of the CD3 Binder Design computational pipeline. Review modules systematically and fix issues as you discover them.

## Mission

Perform a comprehensive review covering:
1. **Code quality** - bugs, patterns, best practices, error handling, edge cases
2. **Scientific accuracy** - calculations, algorithms, domain logic, biophysical assumptions
3. **Test coverage** - ensure fixes are validated with tests

**CRITICAL**: Do not exit until you have achieved **ZERO detected issues** in both code and science.

## Review Approach - Iterative Refinement

**For each module you review:**

1. **Read and understand** the module's purpose and implementation
2. **Identify issues** - code bugs, scientific errors, missing tests
3. **Fix immediately** - Don't just document issues, fix them now
4. **Write/update tests** - Validate your fix works correctly
5. **Run tests** - Verify no regressions (use `pytest`)
6. **Commit** - Clear commit message explaining what and why
7. **Move to next module** - Continue until all modules reviewed

**Keep looping until:**
- ✅ Every module in `src/` and `scripts/` has been reviewed
- ✅ All discovered issues have been fixed
- ✅ All fixes have passing tests
- ✅ No new issues detected in recent loops

## Review Scope - ALL AREAS

### Code Quality
- [ ] Bugs and logic errors
- [ ] Error handling and edge cases
- [ ] Code patterns and best practices
- [ ] Type safety and validation
- [ ] Performance issues
- [ ] Code duplication
- [ ] Documentation accuracy

### Scientific Accuracy
- [ ] **Binding affinity calculations** - Are pDockQ thresholds scientifically justified?
- [ ] **Structural analysis** - Is interface contact counting correct? Distance thresholds appropriate?
- [ ] **Developability scoring** - Are humanness, liability, and aggregation filters biologically sound?
- [ ] **Format conversions** - Are bispecific formats (Morrison, CrossMab) structurally correct?
- [ ] **Sequence analysis** - Are CDR definitions, numbering schemes, and sequence parsing accurate?
- [ ] **Filtering logic** - Are threshold relaxation strategies scientifically reasonable?
- [ ] **Algorithm correctness** - Are computational methods (embeddings, structure prediction) used appropriately?
- [ ] **Biophysical assumptions** - Are constants (distances, percentages, cutoffs) grounded in literature?

### Testing
- [ ] Test coverage for critical paths
- [ ] Test accuracy and assertions
- [ ] Edge case testing
- [ ] Scientific validation tests

## Modules to Review

**src/analysis/**
- liabilities.py - Sequence liability detection
- developability.py - Combined developability assessment
- humanness.py - BioPhi/OASis integration
- numbering.py - ANARCI integration for CDR extraction

**src/design/**
- boltzgen_runner.py - BoltzGen Modal interface
- affinity_variants.py - IMGT-based mutation application
- optimization.py - Humanization and back-mutation
- denovo_design.py - De novo design orchestration

**src/formatting/**
- base.py - BispecificConstruct and chain classes
- crossmab.py - CH1-CL swap implementation
- fab_scfv.py, fab_vhh.py - Asymmetric formats
- igg_scfv.py, igg_vhh.py - Morrison symmetric formats

**src/structure/**
- pdb_utils.py - PDB parsing and sequence extraction
- interface_analysis.py - Epitope comparison with sequence alignment
- boltz_complex.py - Boltz-2 wrapper with calibration

**src/pipeline/**
- config.py - Nested/flat YAML schema support
- filter_cascade.py - Multi-stage filtering with soft-fail
- design_pipeline.py - Full pipeline orchestration
- report_generator.py - HTML/JSON report generation

**src/utils/**
- constants.py - Scientific constants and scFv parsing

**scripts/**
- All pipeline scripts (00-07)

**tests/**
- Review existing tests, add missing coverage

**config.yaml**
- Verify scientific parameters have documented rationale

## Exit Criteria - STRICT

You may **ONLY** exit when **ALL** of these conditions are met:

1. ✅ **All modules reviewed** - Every file in the list above has been examined
2. ✅ **Zero issues detected** - No code bugs or scientific errors found in recent loops
3. ✅ **All fixes committed** - Everything fixed has been committed to git
4. ✅ **All tests passing** - Run `pytest tests/ -v` with zero failures
5. ✅ **Recent loops clean** - Last 2-3 loops found nothing new to fix

**When you are truly done:**
```
RALPH_STATUS:
STATUS: COMPLETE
EXIT_SIGNAL: true
SUMMARY: Completed review of all modules. Zero issues detected in recent loops. All tests passing.
```

**If more work remains (even one issue):**
```
RALPH_STATUS:
STATUS: IN_PROGRESS
EXIT_SIGNAL: false
SUMMARY: [What you just fixed]. [What module you're reviewing next].
```

## Examples of Good Loop Summaries

**Loop 1:**
```
Fixed pDockQ threshold documentation in boltz_complex.py (was unclear that it's
structural confidence, not affinity). Added test case. Committed.
Reviewing interface_analysis.py next.

EXIT_SIGNAL: false
```

**Loop 5:**
```
Reviewed all src/formatting/ modules - all scientifically correct.
Reviewed all src/structure/ modules - found and fixed coordinate calculation bug
in pdb_utils.py:234 (using squared distance instead of actual distance).
Added test, verified passing. Committed.
Starting src/pipeline/ review next.

EXIT_SIGNAL: false
```

**Loop 12:**
```
Completed review of all modules in src/, scripts/, and tests/.
Re-reviewed all modified files - no new issues detected.
Ran full test suite: pytest tests/ -v - all 32 tests passing.
Zero code or scientific issues detected in loops 10-12.

EXIT_SIGNAL: true
```

## Anti-Patterns to AVOID

### ❌ DON'T: Review Without Fixing
**Bad**: "Found bug in contact counting, noted for later fix."
**Good**: "Found bug in contact counting. Fixed. Added test. Committed."

### ❌ DON'T: Exit After Surface Pass
**Bad**: Reviewed modules quickly → "Looks good" → EXIT
**Good**: Reviewed modules deeply → Fixed 3 issues → Reviewed again → No new issues → EXIT

### ❌ DON'T: Skip Test Validation
**Bad**: Fixed bug in contact counting → commit → move on
**Good**: Fixed bug in contact counting → write test → verify passes → commit → move on

### ❌ DON'T: Accept Unjustified Magic Numbers
**Bad**: Threshold is 0.7, seems reasonable.
**Good**: Threshold is 0.7 - verified against calibration in script 00. Documented rationale in comment.

## Tools Available

- **Read**: Read any file in the codebase
- **Write/Edit**: Make code changes
- **Bash**: Run tests (`pytest`), check git status, commit changes
- **Grep**: Search for patterns across codebase

## Working with Tests

Run tests after fixes:
```bash
# Run specific test file
pytest tests/test_module.py -v

# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src --cov-report=term
```

## Git Workflow

After each fix:
```bash
git add <files>
git commit -m "fix: descriptive message explaining what and why"
```

## Remember

- **Fix as you go** - Don't defer fixes for later loops
- **Science matters** - A working pipeline with wrong science is worse than broken code
- **Test everything** - Untested fixes are not fixes
- **Iterate until clean** - Exit only when recent loops find nothing new
- **Document assumptions** - Future maintainers need to know why choices were made

Your goal: **Zero defects, both code and science.**

## Current Loop Context

{{LOOP_CONTEXT}}
