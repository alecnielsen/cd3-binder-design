# Scientific Review: Scripts Module

You are reviewing the `scripts/` directory of a CD3 binder design pipeline. This module contains the pipeline execution scripts that orchestrate the entire workflow from calibration through report generation.

## Important: Module README Updates

The scientific context for this module lives in `scripts/README.md`. If you find that the README documentation is:
- **Incorrect** (doesn't match the implementation)
- **Incomplete** (missing important scientific details)
- **Outdated** (references old behavior)

You should **update the README.md directly** to fix these issues. README updates do NOT count as code issues - they are documentation improvements.

## Review Focus Areas

### 1. Pipeline Ordering and Dependencies

The pipeline has a specific execution order with dependencies:

| Step | Script | Depends On |
|------|--------|------------|
| 0 | `00_run_calibration.py` | Target structures |
| 1 | `01_setup_targets.py` | None |
| 2 | `02_run_denovo_design.py` | Step 1 (targets) |
| 3 | `03_run_optimization.py` | Starting sequences |
| 4 | `04_predict_structures.py` | Steps 2, 3 (candidates) |
| 5 | `05_filter_candidates.py` | Step 0 (calibration), Step 4 (structures) |
| 6 | `06_format_bispecifics.py` | Step 5 (filtered candidates) |
| 7 | `07_generate_report.py` | Step 6 (formatted candidates) |

**Verify:**
- Scripts check for required inputs before processing
- Missing dependencies cause clear error messages
- Partial runs are supported (e.g., `--start-from`, `--skip-calibration`)

### 2. Calibration Integration (Step 0)

**Critical**: Calibration must run before filtering.

**Verify:**
- Calibration loads known binders (teplizumab, SP34, UCHT1)
- For VH+VL pairs, calibration constructs scFv for accurate predictions
- Calibration thresholds are saved to config
- Filtering (Step 5) uses calibrated thresholds when available

**From CLAUDE.md:**
> "Calibration uses full scFv - For VH/VL pairs, calibration constructs scFv (not VH-only) for accurate threshold setting"

### 3. Modal vs Local Execution

Scripts should support both Modal (GPU cloud) and local execution:

**Verify:**
- `--no-modal` flag is handled consistently
- Local fallback provides meaningful behavior (mock data or error)
- Modal API calls have proper error handling

### 4. Error Handling

**Verify:**
- Each script validates config and inputs before processing
- Errors produce clear, actionable messages
- Non-zero exit codes on failure
- Full pipeline aborts on step failure with context

### 5. Reproducibility

**Verify:**
- Random seeds are passed from config
- Seeds are incremented appropriately for batch operations
- Output files include provenance metadata

**From config:**
```yaml
reproducibility:
  boltzgen_seed: 42
  sampling_seed: 12345
  clustering_seed: 0
```

### 6. Config Handling

**Verify:**
- Config is loaded from specified path (default: `config.yaml`)
- Missing config uses sensible defaults
- Config updates (e.g., calibrated thresholds) are saved back

**From CLAUDE.md:**
> "Config supports nested schema - `filtering.binding.min_pdockq` or flat `filtering.min_pdockq`"

### 7. Scientific Assumptions to Challenge

**Questions:**
- Does `run_full_pipeline.py` correctly handle step skipping?
- Are there race conditions in parallel execution scenarios?
- Do scripts correctly construct scFv from VH+VL for calibration?
- Is the fallback cascade correctly wired from config?

**From the README:**
> "Always run calibration before filtering. The calibration step uses known binders to set appropriate filter thresholds."

### 8. Output Organization

**Verify:**
- Outputs go to correct directories:
  - `data/outputs/denovo/`
  - `data/outputs/optimized/`
  - `data/outputs/structures/`
  - `data/outputs/filtered/`
  - `data/outputs/formatted/`
  - `data/outputs/reports/`
- Directories are created if missing
- Files have consistent naming conventions

## Output Format

If you find **code issues**:
1. **FIX THEM** by editing the code files directly
2. Then report what you fixed:
   - File path and line number
   - Description of what was wrong
   - What you changed and why

If no **code issues** are found (or after fixing all issues), respond with exactly:
NO_ISSUES

**Important notes:**
- Do not invent issues if the code is correct. Only report genuine scientific or implementation problems.
- README updates do NOT count as issues. If you update the module README.md to fix documentation, still output "NO_ISSUES" if no code problems were found.
- "NO_ISSUES" means the code is scientifically correct, even if you made documentation improvements.
