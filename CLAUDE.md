# Claude Project Context: CD3 Binder Design

## Overview

Computational pipeline for designing CD3-binding domains (VHH/Fab) for bispecific antibody therapeutics. See README.md for full documentation.

**Goal**: Design ~10 CD3 binders for experimental testing.

## Key Constraints

- **Licensing**: All tools must be permissive (MIT, BSD, Apache, CC-BY). Never suggest IgFold or NetMHCIIpan.
- **Affinity**: ~50 nM Kd hypothesis (balance efficacy with CRS risk). Cannot be computationally enforced.
- **Compute**: Modal required for GPU (BoltzGen/Boltz-2 need CUDA).

## Design Types

| Type | Method | Output | Status |
|------|--------|--------|--------|
| **VHH** | BoltzGen `nanobody-anything` | Single-domain ~120 aa | ✅ Validated |
| **Fab** | BoltzGen `antibody-anything` CDR redesign | VH ~120 aa + VL ~107 aa | ✅ Working |
| **scFv from known Ab** | Optimization track | VH-linker-VL | ✅ Working |

The Fab CDR redesign uses human antibody scaffolds (adalimumab, belimumab, etc.) and redesigns only the CDR loops, maintaining proper VH/VL structure.

**IMPORTANT**: The optimization track does NOT generate new designs - it reformats known antibody sequences (teplizumab, SP34, UCHT1) for structure prediction and comparison. Only VHH and Fab tracks produce **novel** binders.

## Notes for Claude

Critical implementation details:

1. **All tools must be permissively licensed** - never suggest IgFold or NetMHCIIpan
2. **Affinity target is ~50 nM** - not maximizing affinity, balancing with CRS risk
3. **Morrison format is symmetric** - both heavy chains have scFv/VHH fused to C-terminus
4. **CrossMab uses CH1-CL swap** - on the CD3 arm only, for correct light chain pairing
5. **Modal is required for GPU** - BoltzGen won't work locally without NVIDIA GPU
6. **Placeholder target is trastuzumab (HER2)** - configurable for actual use
7. **VH/VL pairs modeled as scFv** - structure prediction uses VH + linker + VL for paired chains
8. **Analysis modules use functions, not classes** - e.g., `score_humanness_pair()` not `HumannessScorer`
9. **Config supports nested schema** - `filtering.binding.min_pdockq` or flat `filtering.min_pdockq`
10. **scFv parsing available** - `parse_scfv()` and `is_likely_scfv()` in `src/utils/constants.py` extract VH/VL from concatenated scFv sequences
11. **CDR-specific liability filtering** - `CandidateScore` has `cdr_*_count` fields; `allow_deamidation_cdr=False` only rejects CDR liabilities
12. **Calibration uses full scFv** - For VH/VL pairs, calibration constructs scFv (not VH-only) for accurate threshold setting
13. **Contact counts are residue-level** - `num_contacts` counts unique residue pairs, not atomic contacts
14. **Epitope residues configurable** - `config.epitope.okt3_epitope_residues` overrides hardcoded defaults
15. **Aggregation filter active** - CDR-specific: >20% aromatic or 2+ consecutive aromatics. Fallback (no CDRs): >15% aromatic or 3+ consecutive aromatics.
16. **Fallback config wired** - `relax_soft_filters_first` and `max_threshold_relaxation` are now used
17. **CRITICAL: BoltzGen extracts target chain** - Multi-chain PDBs (1XIW has UCHT1 Fab, 1SY6 has OKT3 Fab) would bias designs. `boltzgen_app.py` extracts only the specified target chain before design.
18. **CRITICAL: 1XIW numbering starts at 12** - CD3ε chain A in 1XIW uses PDB residue numbers starting at 12, not 1. OKT3 epitope comparison handles this via sequence alignment. Never assume 1-indexed numbering for PDB data.
19. **1SY6 chain A is a fusion** - 1SY6 uses a CD3γ/ε fusion construct, NOT pure CD3ε. Use for epitope reference but prefer 1XIW for canonical CD3ε sequence.
20. **CRITICAL: BoltzGen YAML format** - Use `sequence: 110..130` (range notation), NOT `sequence: XXX...`. Reference target via `file: { path: target.pdb }`, NOT inline sequence. Invalid keys like `design: true` or `binding: [A]` will cause silent failures.
21. **BoltzGen metrics** - `ipTM` (interface pTM) indicates predicted binding quality, `pTM` indicates binder fold quality. Values 0.2-0.4 are typical for de novo designs.
22. **pDockQ is NOT a native Boltz-2 metric** - It's an AlphaFold-Multimer post-processing score. Boltz-2 outputs pTM, ipTM, pLDDT, PAE instead. Use ipTM for interface quality.
23. **Boltz-2 has affinity prediction (IC50)** - Can be enabled with `--sampling_steps_affinity 200` and YAML `properties: { affinity: { ligand: B } }`. NOT validated for antibodies, not currently used.
24. **Calibration uses scFv-derived constructs** - "teplizumab", "sp34", "ucht1" in code are VH-linker-VL scFvs, NOT full IgG antibodies. See `docs/reference/calibration-methodology.md`.
25. **Structural metrics ≠ affinity** - High pTM/pLDDT/contacts does NOT mean high affinity. These are structural confidence scores. Affinity must be validated experimentally with SPR/BLI.
26. **Fab scaffolds must be downloaded** - Run `python scripts/setup_fab_scaffolds.py` before Fab CDR redesign. Downloads CIF+YAML files to `data/fab_scaffolds/`. Only 12 fully-human scaffolds are used (crenezumab and mab1 excluded).
27. **Denovo output format** - `02_run_denovo_design.py` outputs `{"vhh_designs": [...], "fab_designs": [...]}`. The structure prediction script extracts both arrays.
28. **Use conda environment** - A `cd3-binder` conda env with Python 3.9, ANARCI, and Sapiens is set up. Activate with `conda activate cd3-binder`.
29. **Humanness scoring uses Sapiens** - Neural network model. OASis (database-based) is optional fallback.
30. **ANARCI installed via conda** - CDR numbering and CDR-H3 length analysis work. Located at `/opt/homebrew/Caskroom/miniconda/base/envs/cd3-binder/`.
31. **Optimization track is NOT de novo** - It only reformats known antibody sequences from `data/starting_sequences/*.yaml` as scFv for comparison. It does not generate new binders.
32. **Fab scaffold chain IDs fixed** - `setup_fab_scaffolds.py` now uses correct CIF `struct_asym.id` chain labels (not PDB author chain IDs). 12 fully-human scaffolds verified: adalimumab (B/A), belimumab (B/A), dupilumab (A/B), golimumab (G/D), guselkumab (B/A), necitumumab (B/C), nirsevimab (A/B), sarilumab (D/C), secukinumab (A/B), tezepelumab (C/B), tralokinumab (B/C), ustekinumab (B/A).
33. **Modal must be deployed before running** - Run `modal deploy modal/boltzgen_app.py` and `modal run modal/boltzgen_app.py --download` to set up BoltzGen on Modal.
34. **Modal package required in conda env** - Run `pip install modal` in the cd3-binder conda environment.
35. **Validation test passed (VHH)** - 10 VHH designs generated, predicted, filtered, formatted into 29 bispecific constructs. De novo VHH scores competitively with known antibody scFvs.
36. **Fab CDR redesign working** - VH/VL parsing fixed to identify chains by sequence patterns (EVQL... for VH, DIQM... for VL) rather than entity order.
37. **BoltzGen installed from GitHub main** - PyPI v0.2.0 (Dec 10) had a refolding bug; GitHub main has Dec 17 fix. Modal image uses `git+https://github.com/HannesStark/boltzgen.git@main`.
38. **BoltzGen RMSD filter is strict** - Fab CDR redesigns may fail the 2.5 Å RMSD threshold. With correct scaffold YAML format, ~25% of designs pass. Run 20+ designs to get passing candidates.
39. **Scaffold YAML format is critical** - Must include `structure_groups`, `exclude`, `design_insertions` (all 6 CDRs), and `reset_res_index` sections. Official YAMLs from BoltzGen repo are now used.
40. **ipTM interpretation** - For de novo designs, ipTM 0.2-0.4 is typical. Values at high end (~0.35-0.40) indicate good structural confidence. ipTM is NOT an affinity predictor.
41. **Step 04 extracts interface metrics** - BoltzGen outputs ipTM/pTM, but filtering needs interface_area, num_contacts, and epitope residues. Step 04 re-runs Boltz-2 to extract these.
42. **Fab track validated** - Feb 5 run: 9/10 Fab designs generated with ipTM up to 0.368. One design typically fails BoltzGen's internal RMSD filter.
43. **100x scale run completed** - Feb 6: 192 designs (100 VHH + 92 Fab) → 10 final candidates (6 Fab, 4 VHH). Top score 0.412. Output: `denovo_results_20260205_173718.json`.
44. **Modal timeouts expected at scale** - ~14% failure rate (26/192) due to H100 GPU contention. Failures are intermittent, not sequential. Use batch endpoint or retry.
45. **Filtering uses 5 cascaded filters** - (1) Binding: interface_area ≥2060Å², contacts ≥28. (2) Humanness: Sapiens ≥0.8. (3) Liabilities: no CDR deamidation/isomerization/glycosylation. (4) Developability: CDR-H3 8-20aa, charge -2 to +4. (5) Aggregation: CDR aromatics ≤20%.
46. **Hard vs soft filters** - Hard=FAIL rejects outright, Soft=SOFT_FAIL flags but passes. Oxidation (>2 sites), developability ranges, aggregation are soft filters.
47. **Calibrated thresholds from known binders** - min_interface_area (2060Å²) and min_contacts (28) derived from teplizumab/SP34/UCHT1 scFvs minus margin. pDockQ threshold is 0.0 (disabled).
48. **Liability detection is regex-based** - No ML models. Pattern matching: NG/NS/NT/ND/NH (deamidation), DG/DS/DT/DD/DH/DN (isomerization), N-X-S/T (glycosylation), M/W (oxidation).
49. **Fab designs dominate at scale** - At 100x, 6/10 top candidates are Fab (vs 2/10 at 10x). Fab has better humanness scores due to fully-human scaffolds.
50. **Ranking supports BoltzGen native or worst_metric_rank** - BoltzGen's internal decision tree (ipTM, pTM, PAE, H-bonds, salt bridges, SASA) was experimentally validated at 66% nanobody hit rate. Config: `ranking.method: boltzgen` (preferred) or `worst_metric_rank` (current fallback — used because `boltzgen_rank` data is not available in 100x run). Pipeline applies therapeutic pass/fail filters then sorts by chosen ranking method.
51. **Diversity selection is greedy maximin** - After ranking, iteratively pick candidates maximizing `(1-alpha)*quality + alpha*(1-max_identity_to_selected)` with alpha=0.001. Prevents near-duplicate pairs from occupying multiple slots.
52. **CIF files saved to `data/outputs/structures/cif/`** - Boltz-2 outputs CIF format, now preserved on disk for downstream tools (ANTIPASTI, visualization).
53. **iptm now captured in ComplexPredictionResult** - `iptm` field added to both `ComplexPredictionResult` and `CandidateScore`. Extracted from Boltz-2 results.
54. **Validation step 05b integrates Protenix + affinity scoring** - `scripts/05b_validate_candidates.py` runs after filtering. Uses ProteinMPNN (MIT, inverse folding log-likelihood, Spearman r=0.27-0.41 on AbBiBench) and AntiFold (BSD-3, antibody-specific with nanobody support) for affinity proxy scoring. Protenix (Apache 2.0) provides cross-validation of Boltz-2 structure predictions on Modal H100. Results are informational only — not used for filtering/ranking. All 10 candidates validated; top Protenix ipTM = 0.570 (vhh_1XIW_0006). ANTIPASTI was skipped (R dependency, no VHH support, degrades on predicted structures).
55. **PRODIGY not integrated** - r=0.16 on antibody-antigen complexes (poor). AttABseq requires WT reference (not applicable for de novo). Boltz-2 IC50 is small-molecule only.
56. **Protenix deployed on Modal** - `modal/protenix_app.py` runs Protenix (Apache 2.0) on H100 with CUDA 12.4.1 devel image (required — Protenix JIT-compiles CUDA kernels). Functions: `predict_complex()`, `batch_predict()`, `warmup()`. Deploy: `modal deploy modal/protenix_app.py && modal run modal/protenix_app.py --warmup-flag`. No separate `download_weights` exists; `warmup()` runs a minimal prediction to auto-download and cache weights.
57. **ProteinMPNN + AntiFold for affinity scoring** - `src/analysis/affinity_scoring.py` provides `batch_score_affinity()`. Run locally (CPU). ProteinMPNN (MIT) is best-validated affinity proxy for de novo antibodies. AntiFold (BSD-3) supports nanobodies via `--nanobody_chain`.
58. **Validation is informational only** - Step 05b scores are NOT used for filtering or ranking. They provide additional context for experimental candidate selection. CandidateScore has `proteinmpnn_ll`, `antifold_ll`, `protenix_iptm`, `protenix_ptm`, `protenix_ranking_score` fields.
59. **ANTIPASTI deprecated** - `src/analysis/antipasti_scoring.py` now redirects to `affinity_scoring.py`. Old `AntipastiResult` and `score_affinity()` return deprecation errors.
60. **ValidationConfig in config.yaml** - `validation.enabled`, `validation.run_protenix`, `validation.run_proteinmpnn`, `validation.run_antifold`. `iptm_disagreement_threshold: 0.1` flags candidates where Boltz-2 and Protenix ipTM differ by >0.1.
61. **Protenix requires CUDA devel image** - `nvidia/cuda:12.4.1-devel-ubuntu22.04` base image, not `debian_slim`. Protenix JIT-compiles custom CUDA kernels (e.g., `fast_layer_norm_cuda_v2`) which require `nvcc` compiler. Set `CUDA_HOME=/usr/local/cuda`.
62. **Python 3.9 type syntax** - The conda env uses Python 3.9. Use `Optional[int]` (from `typing`) not `int | None` (requires 3.10+). This applies to all scripts run locally.
63. **Ranking uses worst_metric_rank** - Config currently set to `worst_metric_rank` because `boltzgen_rank` data is not available in the 100x run (predates rank capture). Switch to `boltzgen` for future runs.
64. **Protenix CIF files saved** - `data/outputs/structures/protenix_cif/` contains CIF structures from Protenix cross-validation for all 10 candidates.

## Quick Reference

### Running the Pipeline
```bash
conda activate cd3-binder
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

# One-time setup
modal deploy modal/boltzgen_app.py      # Deploy BoltzGen to Modal
modal run modal/boltzgen_app.py --download  # Download model weights
modal deploy modal/protenix_app.py      # Deploy Protenix to Modal
modal run modal/protenix_app.py --warmup-flag  # Pre-cache Protenix weights
python scripts/setup_fab_scaffolds.py   # Download Fab scaffolds

# Run pipeline
python scripts/00_run_calibration.py    # Run first!
python scripts/run_full_pipeline.py --config config.yaml
```

### Key Directories
- `src/` - Core library modules
- `modal/` - GPU compute deployments (BoltzGen + Protenix)
- `scripts/` - Pipeline execution (00-07, including 05b validation)
- `data/fab_scaffolds/` - Human Fab scaffold files (CIF + YAML)
- `docs/reference/` - Detailed implementation docs

### Fab Scaffolds Available
12 fully-human antibody scaffolds: adalimumab, belimumab, dupilumab, golimumab, guselkumab, necitumumab, nirsevimab, sarilumab, secukinumab, tezepelumab, tralokinumab, ustekinumab

Excluded: crenezumab (humanized, has murine framework residues), mab1 (unclear provenance)

Configure in `config.yaml`:
```yaml
design:
  fab_scaffolds:
    - adalimumab
    - belimumab
    - dupilumab
```

### Future: Affinity Prediction Tools (Not Yet Integrated)
- **Boltz-2 IC50** - Enable with `--sampling_steps_affinity 200` (MIT, not antibody-validated)
- **AttABseq** - Sequence-based ΔΔG prediction (MIT, https://github.com/ruofanjin/AttABseq)
- **Graphinity** - Structure-based ΔΔG prediction (BSD-3, https://github.com/oxpig/Graphinity)
- See HANDOFF.md "Future Enhancements" section for details
