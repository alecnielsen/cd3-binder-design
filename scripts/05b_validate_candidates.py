#!/usr/bin/env python3
"""Step 5b: Validate top candidates with affinity scoring and cross-prediction.

Runs after filtering (step 05) on the final candidates:
1. Affinity proxy scoring (ProteinMPNN + AntiFold log-likelihoods) — local CPU
2. Protenix re-prediction — Modal GPU cross-validation vs Boltz-2
3. Cross-validation analysis — flag ipTM disagreements

Results are informational (included in report), not used for filtering/ranking.

Usage:
    python scripts/05b_validate_candidates.py [--config config.yaml]
"""

import argparse
from pathlib import Path
import json


AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_sequence_from_pdb(pdb_path: str, chain_id: str = "A") -> str:
    """Extract amino acid sequence from a specific chain in PDB file."""
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
    parser = argparse.ArgumentParser(description="Validate candidates")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input filtered candidates JSON")
    parser.add_argument("--output", type=str, default="data/outputs/validated", help="Output dir")
    parser.add_argument("--no-modal", action="store_true", help="Skip Protenix (Modal GPU)")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Candidate Validation (Step 05b)")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    if not config.validation.enabled:
        print("Validation disabled in config. Skipping.")
        return 0

    # Find input
    if args.input:
        input_path = Path(args.input)
    else:
        filtered_path = Path("data/outputs/filtered/filtered_candidates.json")
        if filtered_path.exists():
            input_path = filtered_path
        else:
            print("ERROR: No input found. Run step 05 first or specify --input")
            return 1

    print(f"\nInput: {input_path}")

    # Load candidates
    with open(input_path) as f:
        data = json.load(f)

    candidates_data = data.get("candidates", [])
    filter_stats = data.get("filter_stats", {})
    print(f"Loaded {len(candidates_data)} candidates")

    if not candidates_data:
        print("No candidates to validate.")
        return 0

    # CIF directory for structure-based scoring
    cif_dir = Path("data/outputs/structures/cif")

    # Extract target sequence from first target PDB
    target_sequence = None
    target_structures = config.design.target_structures
    if target_structures:
        primary_target = target_structures[0]
        if Path(primary_target).exists():
            target_sequence = extract_sequence_from_pdb(primary_target)
            print(f"Target sequence: {len(target_sequence)} residues from {primary_target}")

    # ---- 1. Affinity scoring (local, CPU) ----
    affinity_results = {}
    if config.validation.run_proteinmpnn or config.validation.run_antifold:
        print(f"\n--- Affinity Scoring ---")
        from src.analysis.affinity_scoring import batch_score_affinity

        design_ids = [c["candidate_id"] for c in candidates_data]
        binder_types = [c.get("binder_type", "vhh") for c in candidates_data]

        results = batch_score_affinity(
            cif_dir=str(cif_dir),
            design_ids=design_ids,
            binder_types=binder_types,
            run_proteinmpnn=config.validation.run_proteinmpnn,
            run_antifold=config.validation.run_antifold,
        )

        for r in results:
            affinity_results[r.design_id] = r
            status = []
            if r.proteinmpnn_ll is not None:
                status.append(f"MPNN={r.proteinmpnn_ll:.3f}")
            if r.antifold_ll is not None:
                status.append(f"AF={r.antifold_ll:.3f}")
            if r.error:
                status.append(f"err={r.error[:60]}")
            print(f"  {r.design_id}: {', '.join(status) if status else 'no scores'}")

    # ---- 2. Protenix cross-prediction (Modal GPU) ----
    protenix_results = {}
    run_protenix = config.validation.run_protenix and not args.no_modal
    if run_protenix and target_sequence:
        print(f"\n--- Protenix Cross-Validation ---")
        try:
            import modal
            protenix_fn = modal.Function.from_name("protenix-cd3", "predict_complex")

            protenix_cif_dir = Path("data/outputs/structures/protenix_cif")
            protenix_cif_dir.mkdir(parents=True, exist_ok=True)

            for c in candidates_data:
                cid = c["candidate_id"]
                seq = c.get("sequence", "")
                seq_vl = c.get("sequence_vl")
                binder_type = c.get("binder_type", "vhh")

                print(f"  Predicting {cid} with Protenix...")
                try:
                    pred = protenix_fn.remote(
                        binder_sequence=seq,
                        target_sequence=target_sequence,
                        binder_type=binder_type,
                        binder_sequence_vl=seq_vl,
                        seed=config.validation.protenix_seeds[0],
                        use_msa=config.validation.protenix_use_msa,
                    )

                    protenix_results[cid] = pred

                    if pred.get("error"):
                        print(f"    Error: {pred['error'][:80]}")
                    else:
                        iptm = pred.get("iptm")
                        ptm = pred.get("ptm")
                        rs = pred.get("ranking_score")
                        print(f"    ipTM={iptm}, pTM={ptm}, ranking_score={rs}")

                    # Save CIF
                    if pred.get("cif_string"):
                        cif_path = protenix_cif_dir / f"{cid}_protenix.cif"
                        cif_path.write_text(pred["cif_string"])

                except Exception as e:
                    print(f"    Error: {e}")
                    protenix_results[cid] = {"error": str(e)}

        except Exception as e:
            print(f"  Protenix not available: {e}")
            print("  Deploy with: modal deploy modal/protenix_app.py")
            run_protenix = False

    elif run_protenix and not target_sequence:
        print("\nSkipping Protenix: no target sequence available")
        run_protenix = False

    # ---- 3. Cross-validation analysis ----
    print(f"\n--- Cross-Validation Analysis ---")
    threshold = config.validation.iptm_disagreement_threshold

    validated_candidates = []
    disagreement_count = 0

    for c in candidates_data:
        cid = c["candidate_id"]
        validated = dict(c)  # Copy all existing fields

        # Add affinity scores
        validation = {}
        aff = affinity_results.get(cid)
        if aff:
            validation["proteinmpnn_ll"] = aff.proteinmpnn_ll
            validation["antifold_ll"] = aff.antifold_ll

        # Add Protenix scores
        prot = protenix_results.get(cid, {})
        if prot and not prot.get("error"):
            validation["protenix_iptm"] = prot.get("iptm")
            validation["protenix_ptm"] = prot.get("ptm")
            validation["protenix_ranking_score"] = prot.get("ranking_score")

            # Compute ipTM delta (Boltz-2 vs Protenix)
            boltz2_iptm = c.get("binding", {}).get("iptm")
            protenix_iptm = prot.get("iptm")

            if boltz2_iptm is not None and protenix_iptm is not None:
                delta = abs(boltz2_iptm - protenix_iptm)
                validation["boltz2_protenix_iptm_delta"] = round(delta, 4)

                if delta > threshold:
                    disagreement_count += 1
                    validation.setdefault("agreement_flags", []).append(
                        f"ipTM_disagreement: Boltz2={boltz2_iptm:.3f} vs Protenix={protenix_iptm:.3f} (delta={delta:.3f})"
                    )
        elif prot and prot.get("error"):
            validation["protenix_error"] = prot["error"][:200]

        validated["validation"] = validation
        validated_candidates.append(validated)

    # Print summary
    n_with_mpnn = sum(1 for c in validated_candidates if c.get("validation", {}).get("proteinmpnn_ll") is not None)
    n_with_af = sum(1 for c in validated_candidates if c.get("validation", {}).get("antifold_ll") is not None)
    n_with_protenix = sum(1 for c in validated_candidates if c.get("validation", {}).get("protenix_iptm") is not None)

    print(f"\n  ProteinMPNN scores: {n_with_mpnn}/{len(validated_candidates)}")
    print(f"  AntiFold scores: {n_with_af}/{len(validated_candidates)}")
    print(f"  Protenix predictions: {n_with_protenix}/{len(validated_candidates)}")
    if n_with_protenix > 0:
        print(f"  ipTM disagreements (>{threshold}): {disagreement_count}/{n_with_protenix}")

    # ---- Save results ----
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "validated_candidates.json"

    with open(output_path, "w") as f:
        json.dump({
            "candidates": validated_candidates,
            "count": len(validated_candidates),
            "filter_stats": filter_stats,
            "validation_summary": {
                "proteinmpnn_scored": n_with_mpnn,
                "antifold_scored": n_with_af,
                "protenix_predicted": n_with_protenix,
                "iptm_disagreements": disagreement_count,
                "iptm_disagreement_threshold": threshold,
            },
        }, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    print("\n" + "=" * 60)
    print("Validation complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
