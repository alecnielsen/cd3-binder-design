#!/usr/bin/env python3
"""Step 4b: Post-hoc humanization of near-miss candidates using HuDiff.

Identifies candidates that narrowly failed humanness filtering (0.70-0.80)
and generates humanized variants using HuDiff (diffusion-based framework
regeneration). CDR sequences are preserved.

Input: data/outputs/structures/candidates_with_scores.json (from 04a)
Output: data/outputs/structures/candidates_humanized.json

Usage:
    python scripts/04b_humanize_candidates.py [--config config.yaml] [--no-modal]
"""

import argparse
from pathlib import Path
from typing import Optional
import json


def extract_sequence_from_pdb(pdb_path: str, chain_id: str = "A") -> str:
    """Extract amino acid sequence from a specific chain in PDB file."""
    AA_3TO1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    sequence = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                if chain != chain_id:
                    continue
                res_name = line[17:20].strip()
                res_num = int(line[22:26])
                key = (chain, res_num, res_name)
                if key not in seen and res_name in AA_3TO1:
                    seen.add(key)
                    sequence.append((res_num, AA_3TO1[res_name]))
    sequence.sort(key=lambda x: x[0])
    return "".join(aa for _, aa in sequence)


def main():
    parser = argparse.ArgumentParser(description="Post-hoc humanization of near-miss candidates")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input candidates JSON (from step 04a)")
    parser.add_argument("--output", type=str, default="data/outputs/structures/candidates_humanized.json",
                        help="Output file")
    parser.add_argument("--no-modal", action="store_true", help="Skip Modal GPU steps")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Post-Hoc Humanization (Step 04b)")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig
    from src.analysis.humanization import (
        identify_near_misses,
        prepare_hudiff_input,
        parse_hudiff_output,
    )

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    if not config.humanization.enabled:
        print("Humanization disabled in config. Set humanization.enabled: true to enable.")
        return 0

    # Find input
    if args.input:
        input_path = Path(args.input)
    else:
        # Prefer scored candidates from 04a
        scored_path = Path("data/outputs/structures/candidates_with_scores.json")
        structures_path = Path("data/outputs/structures/candidates_with_structures.json")
        if scored_path.exists():
            input_path = scored_path
        elif structures_path.exists():
            input_path = structures_path
        else:
            print("ERROR: No input found. Run step 04a first or specify --input")
            return 1

    print(f"\nInput: {input_path}")

    # Load ALL candidates (including those that failed hard filter)
    # We need access to the full pool to find near-misses
    with open(input_path) as f:
        data = json.load(f)

    # candidates_with_scores.json only contains hard-filter survivors.
    # We need the full candidate pool. Check for the pre-filter input too.
    scored_candidates = data.get("candidates", data if isinstance(data, list) else [data])

    # Also load the full candidate pool (pre-filtering) to find near-misses
    full_pool_path = Path("data/outputs/structures/candidates_with_structures.json")
    if full_pool_path.exists() and str(input_path) != str(full_pool_path):
        with open(full_pool_path) as f:
            full_data = json.load(f)
        all_candidates = full_data.get("candidates", full_data if isinstance(full_data, list) else [full_data])
        print(f"Loaded full candidate pool: {len(all_candidates)} candidates")
    else:
        all_candidates = scored_candidates
        print(f"Using scored candidates as full pool: {len(all_candidates)} candidates")

    # Score humanness for all candidates to identify near-misses
    print(f"\n--- Scoring Humanness for Near-Miss Identification ---")
    from src.analysis.humanness import score_humanness_pair

    for c in all_candidates:
        if "_oasis_score_mean" in c:
            continue  # Already scored
        vh_seq = c.get("sequence") or c.get("vh", "")
        vl_seq = c.get("sequence_vl") or c.get("vl")
        if not vh_seq:
            continue
        try:
            humanness_report = score_humanness_pair(vh_seq, vl_seq)
            c["_oasis_score_mean"] = humanness_report.mean_score
            c["_oasis_score_vh"] = humanness_report.vh_report.oasis_score
            if humanness_report.vl_report:
                c["_oasis_score_vl"] = humanness_report.vl_report.oasis_score
        except Exception:
            pass

    # Identify near-misses
    min_h = config.humanization.min_humanness_for_humanization
    max_h = config.humanization.max_humanness_for_humanization
    near_misses = identify_near_misses(all_candidates, min_score=min_h, max_score=max_h)

    print(f"\nHumanization range: [{min_h:.2f}, {max_h:.2f})")
    print(f"Near-miss candidates: {len(near_misses)}")

    if not near_misses:
        print("No near-miss candidates found. Nothing to humanize.")
        # Save output with original scored candidates only
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump({
                "candidates": scored_candidates,
                "humanized_count": 0,
                "near_miss_count": 0,
            }, f, indent=2)
        print(f"Results saved to: {output_path}")
        return 0

    for nm in near_misses:
        nid = nm.get("design_id", nm.get("name", "unknown"))
        score = nm.get("_oasis_score_mean", "N/A")
        btype = nm.get("binder_type", "vhh")
        print(f"  {nid}: humanness={score:.3f}, type={btype}")

    # ---- 1. Call HuDiff on Modal ----
    use_modal = not args.no_modal
    num_variants = config.humanization.num_variants_per_candidate

    if not use_modal:
        print("\nSkipping HuDiff (--no-modal). Cannot humanize without GPU.")
        return 1

    print(f"\n--- Running HuDiff Humanization ---")
    print(f"  Variants per candidate: {num_variants}")

    hudiff_inputs = [prepare_hudiff_input(nm) for nm in near_misses]
    for inp in hudiff_inputs:
        inp["num_samples"] = num_variants

    try:
        import modal
        batch_humanize_fn = modal.Function.from_name("hudiff-cd3", "batch_humanize")
        hudiff_results = batch_humanize_fn.remote(candidates=hudiff_inputs)
    except Exception as e:
        print(f"ERROR: HuDiff not available: {e}")
        print("Deploy with: modal deploy modal/hudiff_app.py")
        print("Download weights: modal run modal/hudiff_app.py --download")
        return 1

    # Parse results into candidate dicts
    humanized_candidates = []
    for hudiff_result, parent in zip(hudiff_results, near_misses):
        variants = parse_hudiff_output(hudiff_result, parent)
        humanized_candidates.extend(variants)

    total_variants = sum(r.get("num_generated", 0) for r in hudiff_results)
    print(f"\nGenerated {total_variants} humanized variants from {len(near_misses)} near-misses")

    if not humanized_candidates:
        print("WARNING: No humanized variants generated. Check HuDiff deployment.")
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump({
                "candidates": scored_candidates,
                "humanized_count": 0,
                "near_miss_count": len(near_misses),
            }, f, indent=2)
        print(f"Results saved to: {output_path}")
        return 0

    # ---- 2. Re-score humanness locally ----
    print(f"\n--- Re-scoring Humanness (Sapiens) ---")
    passing_humanized = []
    for hc in humanized_candidates:
        vh_seq = hc.get("sequence") or hc.get("vh", "")
        vl_seq = hc.get("sequence_vl") or hc.get("vl")
        try:
            humanness_report = score_humanness_pair(vh_seq, vl_seq)
            hc["_oasis_score_mean"] = humanness_report.mean_score
            hc["_oasis_score_vh"] = humanness_report.vh_report.oasis_score
            if humanness_report.vl_report:
                hc["_oasis_score_vl"] = humanness_report.vl_report.oasis_score

            score = humanness_report.mean_score
            hid = hc.get("design_id", "unknown")
            if score is not None and score >= 0.8:
                passing_humanized.append(hc)
                print(f"  PASS {hid}: humanness={score:.3f}")
            else:
                print(f"  FAIL {hid}: humanness={score:.3f}" if score else f"  FAIL {hid}: no score")
        except Exception as e:
            print(f"  ERROR {hc.get('design_id', 'unknown')}: {e}")

    print(f"\nHumanized variants passing humanness: {len(passing_humanized)}/{len(humanized_candidates)}")

    if not passing_humanized:
        print("WARNING: No humanized variants passed humanness threshold.")
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump({
                "candidates": scored_candidates,
                "humanized_count": 0,
                "near_miss_count": len(near_misses),
                "variants_generated": len(humanized_candidates),
                "variants_passing_humanness": 0,
            }, f, indent=2)
        print(f"Results saved to: {output_path}")
        return 0

    # ---- 3. Re-predict structures (Boltz-2) ----
    if config.humanization.repredict_structures:
        print(f"\n--- Re-predicting Structures (Boltz-2) ---")
        from src.structure.boltz_complex import Boltz2Predictor

        default_target_pdb = None
        if config.design.target_structures:
            default_target_pdb = config.design.target_structures[0]
            if not Path(default_target_pdb).exists():
                default_target_pdb = None

        if not default_target_pdb:
            print("WARNING: No target structure found. Skipping re-prediction.")
        else:
            scfv_linker = config.formatting.scfv_linker
            predictor = Boltz2Predictor(use_modal=True)
            cif_dir = Path("data/outputs/structures/cif")
            cif_dir.mkdir(parents=True, exist_ok=True)

            for i, hc in enumerate(passing_humanized):
                vh = hc.get("vh") or hc.get("sequence", "")
                vl = hc.get("vl") or hc.get("sequence_vl")
                if vh and vl:
                    binder_sequence = vh + scfv_linker + vl
                else:
                    binder_sequence = vh

                design_id = hc.get("design_id", f"humanized_{i}")
                print(f"  Predicting {design_id}...")

                try:
                    result = predictor.predict_complex(
                        binder_sequence=binder_sequence,
                        target_pdb_path=default_target_pdb,
                        seed=config.reproducibility.sampling_seed + 10000 + i,
                    )
                    hc["structure_prediction"] = result.to_dict()
                    hc["structure_prediction"]["binder_sequence_used"] = binder_sequence
                    hc["structure_prediction"]["target_structure"] = default_target_pdb

                    # Save CIF
                    if config.output.export_cif and result.pdb_string:
                        cif_path = cif_dir / f"{design_id}.cif"
                        result.save_cif(str(cif_path))

                    sp = hc["structure_prediction"]
                    print(f"    ipTM={sp.get('iptm', 0):.3f}, "
                          f"area={sp.get('interface_area', 0):.0f}, "
                          f"contacts={sp.get('num_contacts', 0)}")

                    # 3-chain prediction for Fab designs
                    if vh and vl:
                        try:
                            result_3chain = predictor.predict_complex_3chain(
                                vh_sequence=vh,
                                vl_sequence=vl,
                                target_pdb_path=default_target_pdb,
                                seed=config.reproducibility.sampling_seed + 10000 + i,
                            )
                            hc["structure_prediction_3chain"] = result_3chain.to_dict()
                            if config.output.export_cif and result_3chain.pdb_string:
                                cif_path_3c = cif_dir / f"{design_id}_3chain.cif"
                                result_3chain.save_cif(str(cif_path_3c))
                        except Exception as e3:
                            print(f"    3-chain failed: {e3}")

                except Exception as e:
                    print(f"    Failed: {e}")
                    hc["structure_prediction"] = None

    # ---- 4. Re-score with validation tools ----
    if config.humanization.rescore_validation:
        print(f"\n--- Re-scoring with Validation Tools ---")

        # Extract target sequence for Protenix
        target_sequence = None
        if config.design.target_structures:
            primary_target = config.design.target_structures[0]
            if Path(primary_target).exists():
                target_sequence = extract_sequence_from_pdb(primary_target)

        cif_dir = Path("data/outputs/structures/cif")

        # ProteinMPNN + AntiFold (local CPU)
        if config.validation.run_proteinmpnn or config.validation.run_antifold:
            try:
                from src.analysis.affinity_scoring import batch_score_affinity
                design_ids = [hc.get("design_id", "unknown") for hc in passing_humanized]
                binder_types = [hc.get("binder_type", "vhh") for hc in passing_humanized]

                results = batch_score_affinity(
                    cif_dir=str(cif_dir),
                    design_ids=design_ids,
                    binder_types=binder_types,
                    run_proteinmpnn=config.validation.run_proteinmpnn,
                    run_antifold=config.validation.run_antifold,
                )

                for r in results:
                    for hc in passing_humanized:
                        if hc.get("design_id") == r.design_id:
                            vs = hc.setdefault("validation_scores", {})
                            if r.proteinmpnn_ll is not None:
                                vs["proteinmpnn_ll_scfv"] = r.proteinmpnn_ll
                            if r.antifold_ll is not None:
                                vs["antifold_ll_scfv"] = r.antifold_ll
                            status = []
                            if r.proteinmpnn_ll is not None:
                                status.append(f"MPNN={r.proteinmpnn_ll:.3f}")
                            if r.antifold_ll is not None:
                                status.append(f"AF={r.antifold_ll:.3f}")
                            print(f"  {r.design_id}: {', '.join(status)}")
                            break
            except Exception as e:
                print(f"  Affinity scoring failed: {e}")

        # Protenix (Modal GPU)
        if config.validation.run_protenix and target_sequence:
            try:
                import modal
                protenix_fn = modal.Function.from_name("protenix-cd3", "predict_complex")
                protenix_cif_dir = Path("data/outputs/structures/protenix_cif")
                protenix_cif_dir.mkdir(parents=True, exist_ok=True)

                for hc in passing_humanized:
                    hid = hc.get("design_id", "unknown")
                    seq = hc.get("sequence") or hc.get("vh", "")
                    seq_vl = hc.get("sequence_vl") or hc.get("vl")
                    binder_type = hc.get("binder_type", "vhh")

                    print(f"  Protenix: {hid}...")
                    try:
                        pred = protenix_fn.remote(
                            binder_sequence=seq,
                            target_sequence=target_sequence,
                            binder_type=binder_type,
                            binder_sequence_vl=seq_vl,
                            seed=config.validation.protenix_seeds[0],
                            use_msa=config.validation.protenix_use_msa,
                        )
                        vs = hc.setdefault("validation_scores", {})
                        if not pred.get("error"):
                            vs["protenix_iptm"] = pred.get("iptm")
                            vs["protenix_ptm"] = pred.get("ptm")
                            vs["protenix_ranking_score"] = pred.get("ranking_score")
                            print(f"    ipTM={pred.get('iptm')}, pTM={pred.get('ptm')}")
                        else:
                            vs["protenix_error"] = pred["error"][:200]
                            print(f"    Error: {pred['error'][:80]}")

                        if pred.get("cif_string"):
                            cif_path = protenix_cif_dir / f"{hid}_protenix.cif"
                            cif_path.write_text(pred["cif_string"])
                    except Exception as e:
                        print(f"    Error: {e}")
            except Exception as e:
                print(f"  Protenix not available: {e}")

    # ---- 5. Merge and save results ----
    print(f"\n--- Merging Results ---")

    # Strip CIF strings to keep output manageable
    for hc in passing_humanized:
        for key in ["structure_prediction", "structure_prediction_3chain"]:
            sp = hc.get(key)
            if sp and "cif_string" in sp:
                del sp["cif_string"]

    # Merge: original scored candidates + humanized variants
    merged = list(scored_candidates) + passing_humanized

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump({
            "candidates": merged,
            "count": len(merged),
            "humanization_stats": {
                "near_miss_count": len(near_misses),
                "variants_generated": len(humanized_candidates),
                "variants_passing_humanness": len(passing_humanized),
                "humanness_range": f"[{min_h:.2f}, {max_h:.2f})",
                "original_candidates": len(scored_candidates),
                "total_after_merge": len(merged),
            },
        }, f, indent=2)

    print(f"\nResults:")
    print(f"  Original scored candidates: {len(scored_candidates)}")
    print(f"  Near-misses humanized: {len(near_misses)}")
    print(f"  Humanized variants passing: {len(passing_humanized)}")
    print(f"  Total merged candidates: {len(merged)}")
    print(f"\nSaved to: {output_path}")

    print("\n" + "=" * 60)
    print("Humanization complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
