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
54. **Validation step 05b is cross-validation only** - `scripts/05b_validate_candidates.py` runs after filtering. Compares Boltz-2 vs Protenix ipTM, flags disagreements. ProteinMPNN + AntiFold scoring moved to step 04a.
55. **PRODIGY not integrated** - r=0.16 on antibody-antigen complexes (poor). AttABseq requires WT reference (not applicable for de novo). Boltz-2 IC50 is small-molecule only.
56. **Protenix deployed on Modal** - `modal/protenix_app.py` runs Protenix (Apache 2.0) on H100 with CUDA 12.4.1 devel image (required — Protenix JIT-compiles CUDA kernels). Functions: `predict_complex()`, `batch_predict()`, `warmup()`. Deploy: `modal deploy modal/protenix_app.py && modal run modal/protenix_app.py --warmup-flag`. No separate `download_weights` exists; `warmup()` runs a minimal prediction to auto-download and cache weights.
57. **ProteinMPNN + AntiFold for affinity scoring** - `src/analysis/affinity_scoring.py` provides `batch_score_affinity()`. Run locally (CPU). ProteinMPNN uses `parse_PDB` + `tied_featurize` + model forward pass (requires CIF→PDB conversion via BioPython). AntiFold uses `load_model` + `get_pdbs_logits` with DataFrame chain mapping (`custom_chain_mode=True` for VHH/scFv). Both installed: `pip install proteinmpnn antifold`.
58. **Validation scores used for ranking** - ProteinMPNN, AntiFold, and Protenix scores from step 04a are used as ranking metrics in worst_metric_rank (proteinmpnn_ll weight=1, antifold_ll weight=2, protenix_iptm weight=1). CandidateScore has `proteinmpnn_ll`, `antifold_ll`, `protenix_iptm`, `protenix_ptm`, `protenix_ranking_score` fields.
59. **ANTIPASTI deprecated** - `src/analysis/antipasti_scoring.py` now redirects to `affinity_scoring.py`. Old `AntipastiResult` and `score_affinity()` return deprecation errors.
60. **ValidationConfig in config.yaml** - `validation.enabled`, `validation.run_protenix`, `validation.run_proteinmpnn`, `validation.run_antifold`. `iptm_disagreement_threshold: 0.1` flags candidates where Boltz-2 and Protenix ipTM differ by >0.1.
61. **Protenix requires CUDA devel image** - `nvidia/cuda:12.4.1-devel-ubuntu22.04` base image, not `debian_slim`. Protenix JIT-compiles custom CUDA kernels (e.g., `fast_layer_norm_cuda_v2`) which require `nvcc` compiler. Set `CUDA_HOME=/usr/local/cuda`.
62. **Python 3.9 type syntax** - The conda env uses Python 3.9. Use `Optional[int]` (from `typing`) not `int | None` (requires 3.10+). This applies to all scripts run locally.
63. **Ranking uses boltzgen** - Config set to `method: boltzgen` with `secondary_method: worst_metric_rank` fallback. As of Feb 11 run, all designs have `boltzgen_rank` data so primary method is active.
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
modal deploy modal/hudiff_app.py        # Deploy HuDiff to Modal
modal run modal/hudiff_app.py --download  # Download HuDiff checkpoints (~2 GB)
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

### Pipeline Steps (Updated)
```
00 → 01 → 02 → 03 → 04 → 04a → 04b → 05 → 05b → 06 → 07
```
- **Step 04a** (`scripts/04a_score_candidates.py`): Hard pre-filter → ProteinMPNN + AntiFold + Protenix scoring on survivors. Output: `candidates_with_scores.json`.
- **Step 04b** (`scripts/04b_humanize_candidates.py`): Post-hoc humanization of near-miss candidates (humanness 0.70-0.80) using HuDiff on Modal. Re-scores humanness, re-predicts structures, re-scores with validation tools. Output: `candidates_humanized.json`.
- **Step 05** now reads humanized output from 04b (or scored from 04a) and uses validation metrics (proteinmpnn_ll, antifold_ll, protenix_iptm) in worst_metric_rank ranking.
- **Step 05b** is now lightweight cross-validation only (no longer runs affinity scoring — moved to 04a).
- **Dual Fab prediction**: Step 04 predicts Fabs as BOTH scFv (2-chain) and VH+VL+target (3-chain). Both metric sets stored.
- **Ranking auto-fallback**: `method: boltzgen` with `secondary_method: worst_metric_rank`. Falls back automatically if no boltzgen_rank data.
- **Calibration baselines**: Step 00 runs ProteinMPNN/AntiFold/Protenix on controls when `run_validation_baselines: true`.
65. **3-chain prediction via Modal** — `predict_complex_multichain()` in `modal/boltz2_app.py` and `predict_complex_3chain()` in `Boltz2Predictor`. Chain assignment: A=target, B=VH, C=VL. Interface metrics union chains B+C as binder.
66. **ComplexPredictionResult.prediction_mode** — "scfv" or "3chain" field tracks which prediction mode produced the result.
67. **Step 04a pre-filters before scoring** — Only candidates passing hard filters (binding, humanness, CDR liabilities) get expensive validation scoring. Reduces ~200 → ~30-50 candidates.
68. **Validation scores used for ranking** — `proteinmpnn_ll`, `antifold_ll`, `protenix_iptm` fields on both `CandidateScore` and `RankedCandidate`. Included in worst_metric_rank with configurable weights.
69. **RankingConfig.secondary_method** — Fallback ranking method when primary method's data unavailable. Default: "worst_metric_rank".
70. **run_calibration() returns cif_strings** — `calibration["cif_strings"]` dict maps control name → CIF content. Saved to `data/outputs/calibration_cif/`.
71. **ProteinMPNN requires CIF→PDB conversion** — `parse_PDB` only handles PDB format. `_cif_to_pdb()` in `affinity_scoring.py` uses BioPython `MMCIFParser` + `PDBIO` to convert. Temp PDB file cleaned up after parsing.
72. **AntiFold requires custom_chain_mode for single-chain** — For VHH/scFv (single binder chain B), AntiFold needs `custom_chain_mode=True` and `Lchain=None` in the DataFrame. Without this, AntiFold errors on missing Lchain column. For 3-chain Fab files, use standard mode with `Hchain="B"`, `Lchain="C"`.
73. **AntiFold residue column is pdb_res** — AntiFold logits DataFrame uses `pdb_res` column for residue identity (not `wt`). NLL computed by indexing into amino acid logit columns per residue.
74. **ProteinMPNN + AntiFold installed in conda env** — `pip install proteinmpnn antifold` in cd3-binder env. Note: this downgrades torch to 2.3.1 (from 2.8.0) — no impact on local CPU scoring.
75. **Test run validated (2026-02-09)** — 5 VHH + 5 Fab → 250 candidates → 5 passed hard filters → all 5 scored with MPNN (1.10-1.36) + AntiFold (0.25-0.65) + Protenix (0.21-0.34). Full pipeline end-to-end with all scoring tools working.
76. **NLL metrics are lower=better** — Both ProteinMPNN and AntiFold output negative log-likelihood per residue. Lower values indicate better structural fit. The ranking code has a `lower_is_better` set for `{"proteinmpnn_ll", "antifold_ll"}` that reverses sort direction.
77. **AntiFold NLL is NOT negated** — `score_antifold()` returns raw NLL (positive values ~0.3-0.5 for controls). Previously was negated to make "higher=better" but this was inconsistent with ProteinMPNN. Now both are consistently lower=better.
78. **YAML sequences have whitespace** — Starting sequence YAML files use `>-` folded scalar format, which converts newlines to spaces. `load_starting_sequence()` strips all whitespace from sequences. Without this, Protenix silently fails (Boltz-2 tolerates spaces).
79. **Protenix output path** — Protenix v1.0+ outputs to `output_dir/<name>/seed_<seed>/predictions/`. Files: `<name>_summary_confidence_sample_0.json` (metrics) and `<name>_sample_0.cif` (structure). The pLDDT key is `"plddt"` (0-100 scale), normalized to 0-1 in our parsing.
80. **Calibration baselines complete** — All 3 validation tools scored on 3 controls. Teplizumab Protenix ipTM=0.855 (strong, consistent with high-affinity known binder). Baselines in `calibration.json` under `validation_baselines` with `proteinmpnn_ll`, `antifold_ll`, `protenix_iptm`.
81. **CRITICAL: BoltzGen Fab output uses vh_sequence/vl_sequence keys** — BoltzGen outputs `vh_sequence` and `vl_sequence` for Fab designs, but the rest of the pipeline expects `vh` and `vl`. Step 04 now normalizes these keys at load time. Without this fix, Fab VL sequences are silently dropped: Fabs get predicted as VH-only (not scFv), 3-chain prediction is skipped, and no scFv-based bispecific formats are generated.
82. **Step 04 scFv construction** — When a candidate has both `vh` and `vl`, step 04 constructs scFv (VH-linker-VL) for 2-chain prediction. Previously it used the raw `sequence` field (VH-only for Fabs) which produced incorrect single-chain predictions for paired antibodies.
83. **Modal gRPC deadline on 100+ Fab designs** — `run_boltzgen_fab` with 100 designs × 12 scaffolds × 2 targets causes gRPC `Deadline exceeded` on the client side. 50+50 is the reliable config. The Modal function has 90-min server timeout but the gRPC client deadline is shorter.
84. **Step 04 processes all cumulative denovo results** — Structure prediction loads ALL `denovo_results_*.json` files and predicts/re-predicts everything. This ensures consistent metrics but is slow for large cumulative sets. Archive old result files if only new designs are needed.
85. **BoltzGen ranking now primary** — As of Feb 11 run, all designs have `boltzgen_rank` data. Config uses `method: boltzgen` with no fallback needed. BoltzGen's decision tree (ipTM, pTM, PAE, H-bonds, salt bridges, SASA) is experimentally validated at 66% nanobody hit rate.
86. **Config reduced to 50+50** — `num_vhh_designs: 50`, `num_fab_designs: 50`. Produces ~95 designs (~10% Fab RMSD filter failures), plenty for 10 final candidates.
87. **Humanness is the yield bottleneck** — 88% of new designs pass binding filters, but only 17% pass humanness ≥0.8. BoltzGen CDR redesign on fully-human scaffolds still introduces non-human CDR motifs. Many near-misses at 0.78–0.80. Lowering threshold to 0.78 would recover ~11 additional candidates per 50+50 run.
88. **Protenix ipTM 0.69 for Fabs** — Feb 11 run Fab candidates scored 0.686 and 0.693 Protenix ipTM, the highest observed across all runs (previous best: 0.570 VHH). This suggests the VL fix + scFv construction produces better-quality predictions.
89. **Filtering uses scFv metrics, not 3-chain** — Step 04a hard filter checks `structure_prediction` (scFv/2-chain mode). 3-chain predictions are supplementary scoring only. Both Fabs in final candidates pass binding on scFv alone.
90. **Feb 11 run: 3 final candidates** — 387 cumulative candidates → 3 passed all hard filters (2 Fab, 1 VHH). Bottleneck: 152 rejected for humanness, 27 for interface area, 18 for CDR deamidation.
91. **HuDiff post-hoc humanization integrated** — Step 04b uses HuDiff (AFL-3.0, Tencent AI4S) to humanize near-miss candidates (humanness 0.70-0.80). Diffusion-based framework regeneration preserves CDRs while replacing non-human framework residues. Deployed on Modal (A100). Config: `humanization.enabled: true`.
92. **HuDiff has separate Ab/Nb models** — `hudiffab.pt` for VH+VL pairs (antibody), `hudiffnb.pt` for VHH (nanobody). Both downloaded from `cloud77/HuDiff` on HuggingFace (~2 GB tar.gz). Deploy: `modal deploy modal/hudiff_app.py && modal run modal/hudiff_app.py --download`.
97. **HuDiff-Ab has no internal humanness scorer** — The antibody model is pure diffusion (relies on learned human framework distribution from OAS training data). Only HuDiff-Nb bundles AbNatiV for explicit humanness guidance. This means HuDiff-Ab outputs are *likely* human but not *guaranteed* to pass any threshold — Sapiens re-scoring in step 04b is essential.
93. **HuDiff requires ANARCI** — Uses `abnumber` for IMGT numbering of input sequences. ANARCI is installed in the Modal image (not in local conda env). All HuDiff inference runs on Modal.
94. **Step 04b loads full candidate pool** — Reads `candidates_with_structures.json` (all candidates) to find near-misses, not just `candidates_with_scores.json` (hard-filter survivors). Near-misses are candidates that failed humanness but may have passed binding.
95. **Humanized candidates tagged** — `source: "hudiff_humanized"`, `parent_design_id` linking to original, `humanization_mutations` list with position/chain/original/mutated. These go through full re-prediction and re-scoring before merging.
96. **Step 05 prefers humanized input** — If `candidates_humanized.json` exists (from 04b), step 05 loads it instead of `candidates_with_scores.json`. The humanized file contains both original scored candidates and new humanized variants.
97. **HuDiff production run complete** — 275 near-misses humanized → 1375 HuDiff variants → 378 pass humanness ≥0.8 (27.5% success rate). VHH humanization works much better than Fab (~40% vs ~15% pass rate). Total candidate pool: 381 (3 original + 378 humanized).
98. **Protenix parallel via spawn()** — Step 04b uses `protenix_fn.spawn()` + `handle.get()` for parallel Protenix cross-validation. 378 predictions complete in ~2 hours vs ~83 hours serial. Some GPU jobs timeout at 1500s (~10% failure rate) but results are collected for the rest.
99. **Humanized Protenix scores dramatically higher** — Best humanized: vhh_1XIW_0007_hudiff_2 ipTM=0.937, vhh_1XIW_0034_hudiff_0 ipTM=0.876, fab_1SY6_0007_hudiff_0 ipTM=0.812. Previous best was 0.693. Humanized framework regions may improve fold confidence.
100. **HuDiff pymol/patent_eval mocked** — HuDiff imports `pymol.cmd` (dead import) and `patent_eval` (eval-only). Both mocked with `MagicMock` in Modal functions to avoid installing unnecessary heavy dependencies.
101. **Pipeline step ordering updated** — `00 → 01 → 02 → 03 → 04 → 04a → 04b → 05 → 05b → 06 → 07`. Step 04b is now included in `run_full_pipeline.py`.
102. **Step 05 filtering: 10 final candidates** — 381 humanized pool → 130 pass all hard filters → 10 after diversity selection. 2 original candidates (vhh_1XIW_0001, fab_1XIW_0003) + 8 HuDiff-humanized variants. Pipeline goal of ~10 candidates achieved.
103. **BoltzGen rank only on originals** — Humanized variants lack `boltzgen_rank` data (re-predicted via Boltz-2, not BoltzGen). Only 2/130 passing candidates have BoltzGen rank. Ranking falls back to `secondary_method: worst_metric_rank` for the rest, then diversity selection picks final 10.
104. **Step 05b: 6/10 ipTM disagreement** — Cross-validation found 6/10 candidates have Boltz-2 vs Protenix ipTM disagreement >0.1. The 4 with agreement are more trustworthy for experimental prioritization.
105. **Step 06: 26 bispecific constructs** — 5 formats (crossmab, fab_scfv, fab_vhh, igg_scfv, igg_vhh). scFv/crossmab formats only for 2 Fab candidates; VHH formats for all 10. Tumor arm: trastuzumab (HER2).
106. **Full pipeline complete Feb 14** — Steps 00-07 all passing. 10 candidates, 26 bispecific constructs, HTML report generated. Ready for experimental validation.
107. **Boltz-2/Protenix agreement analysis** — 4/10 candidates have agreement (ipTM delta ≤0.1): fab_1XIW_0003, vhh_1XIW_0019_hudiff_0, vhh_1SY6_0013_hudiff_0, vhh_1SY6_0019_hudiff_2. These are highest-confidence for experimental prioritization.
108. **CDR preservation verified** — ANARCI IMGT numbering confirms 57 framework mutations, 0 CDR mutations across all 8 humanized variants. HuDiff `sample_method: FR` preserves CDRs as designed.
109. **Canonical disulfide restored by HuDiff** — Common mutations at IMGT positions 23 and 104 introduce Cys residues, restoring the canonical human VH intradomain disulfide bond absent in camelid VHH frameworks.
110. **fab_1XIW_0003 is top candidate** — Highest ipTM from both Boltz-2 (0.720) and Protenix (0.686), best humanness (0.830), Fab format with human scaffold. Only Fab in the agreement set.

### Future: Affinity Prediction Tools (Not Yet Integrated)
- **Boltz-2 IC50** - Enable with `--sampling_steps_affinity 200` (MIT, not antibody-validated)
- **AttABseq** - Sequence-based ΔΔG prediction (MIT, https://github.com/ruofanjin/AttABseq)
- **Graphinity** - Structure-based ΔΔG prediction (BSD-3, https://github.com/oxpig/Graphinity)
- See HANDOFF.md "Future Enhancements" section for details
