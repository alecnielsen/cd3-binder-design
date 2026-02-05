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

## Quick Reference

### Running the Pipeline
```bash
conda activate cd3-binder
export PYTHONPATH=/Users/alec/kernel/cd3-binder-design

# One-time setup
modal deploy modal/boltzgen_app.py      # Deploy BoltzGen to Modal
modal run modal/boltzgen_app.py --download  # Download model weights
python scripts/setup_fab_scaffolds.py   # Download Fab scaffolds

# Run pipeline
python scripts/00_run_calibration.py    # Run first!
python scripts/run_full_pipeline.py --config config.yaml
```

### Key Directories
- `src/` - Core library modules
- `modal/` - GPU compute deployments
- `scripts/` - Pipeline execution (00-07)
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
