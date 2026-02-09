#!/usr/bin/env python3
"""Step 4a: Pre-filter and score candidates with validation tools.

Applies hard filters to reduce ~200 candidates to ~30-50, then runs
ProteinMPNN, AntiFold, and Protenix on survivors. All scores are stored
for downstream ranking.

Input: data/outputs/structures/candidates_with_structures.json
Output: data/outputs/structures/candidates_with_scores.json

Usage:
    python scripts/04a_score_candidates.py [--config config.yaml] [--no-modal]
"""

import argparse
from pathlib import Path
from typing import Optional
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


def quick_hard_filter(candidate: dict, thresholds: dict, config) -> Optional[str]:
    """Apply cheapest hard filters. Returns failure reason or None if passes.

    Only checks binding quality, humanness, and CDR liabilities.
    Skips soft filters (those are for final ranking).
    """
    sp = candidate.get("structure_prediction") or {}

    # Binding: interface_area and contacts
    interface_area = sp.get("interface_area")
    if interface_area is not None and interface_area < thresholds["min_interface_area"]:
        return f"interface_area {interface_area:.0f} < {thresholds['min_interface_area']}"

    num_contacts = sp.get("num_contacts")
    if num_contacts is not None and num_contacts < thresholds["min_contacts"]:
        return f"num_contacts {num_contacts} < {thresholds['min_contacts']}"

    # Humanness
    vh_seq = candidate.get("sequence") or candidate.get("vh", "")
    vl_seq = candidate.get("sequence_vl") or candidate.get("vl")

    try:
        from src.analysis.humanness import score_humanness_pair
        humanness_report = score_humanness_pair(vh_seq, vl_seq)
        oasis_mean = humanness_report.mean_score
        if oasis_mean is not None and oasis_mean < 0.8:
            return f"humanness {oasis_mean:.3f} < 0.8"
        candidate["_oasis_score_mean"] = oasis_mean
        candidate["_oasis_score_vh"] = humanness_report.vh_report.oasis_score
        if humanness_report.vl_report:
            candidate["_oasis_score_vl"] = humanness_report.vl_report.oasis_score
    except Exception:
        pass  # Can't assess humanness, let it pass

    # CDR liabilities
    try:
        from src.analysis.liabilities import LiabilityScanner
        scanner = LiabilityScanner()
        vh_report = scanner.scan_with_cdr_detection(vh_seq, chain_type="H")
        cdr_deamidation = sum(1 for s in vh_report.deamidation_sites if s.in_cdr)
        cdr_isomerization = sum(1 for s in vh_report.isomerization_sites if s.in_cdr)
        cdr_glycosylation = sum(1 for s in vh_report.glycosylation_sites if s.in_cdr)

        if vl_seq:
            vl_report = scanner.scan_with_cdr_detection(vl_seq, chain_type="L")
            cdr_deamidation += sum(1 for s in vl_report.deamidation_sites if s.in_cdr)
            cdr_isomerization += sum(1 for s in vl_report.isomerization_sites if s.in_cdr)
            cdr_glycosylation += sum(1 for s in vl_report.glycosylation_sites if s.in_cdr)

        if not getattr(config.filtering, "allow_deamidation_cdr", False) and cdr_deamidation > 0:
            return f"CDR deamidation ({cdr_deamidation} sites)"
        if not getattr(config.filtering, "allow_isomerization_cdr", False) and cdr_isomerization > 0:
            return f"CDR isomerization ({cdr_isomerization} sites)"
        if not getattr(config.filtering, "allow_glycosylation_cdr", False) and cdr_glycosylation > 0:
            return f"CDR glycosylation ({cdr_glycosylation} sites)"
    except Exception:
        pass  # Can't assess liabilities, let it pass

    return None


def main():
    parser = argparse.ArgumentParser(description="Pre-filter and score candidates")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input candidates JSON")
    parser.add_argument("--output", type=str, default="data/outputs/structures/candidates_with_scores.json",
                        help="Output file")
    parser.add_argument("--no-modal", action="store_true", help="Skip Protenix (Modal GPU)")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Score Candidates (Step 04a)")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    # Find input
    if args.input:
        input_path = Path(args.input)
    else:
        default_path = Path("data/outputs/structures/candidates_with_structures.json")
        if default_path.exists():
            input_path = default_path
        else:
            print("ERROR: No input found. Run step 04 first or specify --input")
            return 1

    print(f"\nInput: {input_path}")

    # Load candidates
    with open(input_path) as f:
        data = json.load(f)

    candidates = data.get("candidates", data if isinstance(data, list) else [data])
    print(f"Loaded {len(candidates)} candidates")

    # Get thresholds
    thresholds = config.get_effective_thresholds()
    print(f"\nHard filter thresholds:")
    print(f"  min_interface_area: {thresholds['min_interface_area']}")
    print(f"  min_contacts: {thresholds['min_contacts']}")
    print(f"  min_oasis_score: 0.8")

    # ---- 1. Quick hard filter ----
    print(f"\n--- Quick Hard Filter ---")
    survivors = []
    rejected_reasons = {}
    for c in candidates:
        reason = quick_hard_filter(c, thresholds, config)
        if reason is None:
            survivors.append(c)
        else:
            cid = c.get("design_id", c.get("name", "unknown"))
            rejected_reasons[cid] = reason

    print(f"  Passed: {len(survivors)}/{len(candidates)}")
    if rejected_reasons:
        # Summarize rejection reasons
        reason_counts = {}
        for reason in rejected_reasons.values():
            key = reason.split("(")[0].strip() if "(" in reason else reason.split()[0]
            reason_counts[key] = reason_counts.get(key, 0) + 1
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"    Rejected ({reason}): {count}")

    if not survivors:
        print("ERROR: No candidates passed hard filters")
        return 1

    # CIF directory for structure-based scoring
    cif_dir = Path("data/outputs/structures/cif")

    # Extract target sequence for Protenix
    target_sequence = None
    target_structures = config.design.target_structures
    if target_structures:
        primary_target = target_structures[0]
        if Path(primary_target).exists():
            target_sequence = extract_sequence_from_pdb(primary_target)
            print(f"\nTarget sequence: {len(target_sequence)} residues from {primary_target}")

    # ---- 2. Affinity scoring (local CPU) ----
    print(f"\n--- Affinity Scoring ---")
    affinity_results = {}
    if config.validation.run_proteinmpnn or config.validation.run_antifold:
        from src.analysis.affinity_scoring import batch_score_affinity

        design_ids = [c.get("design_id", c.get("name", "unknown")) for c in survivors]
        binder_types = [c.get("binder_type", "vhh") for c in survivors]

        # Score scFv CIF files
        results = batch_score_affinity(
            cif_dir=str(cif_dir),
            design_ids=design_ids,
            binder_types=binder_types,
            run_proteinmpnn=config.validation.run_proteinmpnn,
            run_antifold=config.validation.run_antifold,
        )

        for r in results:
            affinity_results[r.design_id] = {"scfv": r}
            status = []
            if r.proteinmpnn_ll is not None:
                status.append(f"MPNN={r.proteinmpnn_ll:.3f}")
            if r.antifold_ll is not None:
                status.append(f"AF={r.antifold_ll:.3f}")
            if r.error:
                status.append(f"err={r.error[:40]}")
            print(f"  {r.design_id}: {', '.join(status) if status else 'no scores'}")

        # Score 3-chain CIF files for Fab designs
        fab_ids = [did for did, c in zip(design_ids, survivors)
                   if c.get("structure_prediction_3chain")]
        if fab_ids:
            fab_types = ["fab"] * len(fab_ids)
            # 3-chain CIF files are named {id}_3chain.cif
            cif_dir_3chain = str(cif_dir)
            results_3chain = []
            for fid in fab_ids:
                cif_3chain = cif_dir / f"{fid}_3chain.cif"
                if cif_3chain.exists():
                    from src.analysis.affinity_scoring import score_proteinmpnn, score_antifold
                    merged_3c = {"proteinmpnn_ll": None, "antifold_ll": None}
                    if config.validation.run_proteinmpnn:
                        r = score_proteinmpnn(str(cif_3chain), fid)
                        merged_3c["proteinmpnn_ll"] = r.proteinmpnn_ll
                    if config.validation.run_antifold:
                        r = score_antifold(str(cif_3chain), fid, "fab")
                        merged_3c["antifold_ll"] = r.antifold_ll
                    if fid in affinity_results:
                        affinity_results[fid]["3chain"] = merged_3c
                    print(f"  {fid} (3chain): MPNN={merged_3c['proteinmpnn_ll']}, AF={merged_3c['antifold_ll']}")

    # ---- 3. Protenix cross-validation (Modal GPU) ----
    protenix_results = {}
    run_protenix = config.validation.run_protenix and not args.no_modal
    if run_protenix and target_sequence:
        print(f"\n--- Protenix Cross-Validation ---")
        try:
            import modal
            protenix_fn = modal.Function.from_name("protenix-cd3", "predict_complex")

            protenix_cif_dir = Path("data/outputs/structures/protenix_cif")
            protenix_cif_dir.mkdir(parents=True, exist_ok=True)

            for c in survivors:
                cid = c.get("design_id", c.get("name", "unknown"))
                seq = c.get("sequence") or c.get("vh", "")
                seq_vl = c.get("sequence_vl") or c.get("vl")
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

    # ---- 4. Merge all scores into candidates ----
    print(f"\n--- Merging Scores ---")
    scored_candidates = []
    for c in survivors:
        cid = c.get("design_id", c.get("name", "unknown"))
        scored = dict(c)

        # Build validation_scores dict
        validation_scores = {}

        # Affinity scores (scFv)
        aff = affinity_results.get(cid, {})
        scfv_aff = aff.get("scfv")
        if scfv_aff:
            if scfv_aff.proteinmpnn_ll is not None:
                validation_scores["proteinmpnn_ll_scfv"] = scfv_aff.proteinmpnn_ll
            if scfv_aff.antifold_ll is not None:
                validation_scores["antifold_ll_scfv"] = scfv_aff.antifold_ll

        # Affinity scores (3-chain)
        threechain_aff = aff.get("3chain", {})
        if threechain_aff.get("proteinmpnn_ll") is not None:
            validation_scores["proteinmpnn_ll_3chain"] = threechain_aff["proteinmpnn_ll"]
        if threechain_aff.get("antifold_ll") is not None:
            validation_scores["antifold_ll_3chain"] = threechain_aff["antifold_ll"]

        # Protenix scores
        prot = protenix_results.get(cid, {})
        if prot and not prot.get("error"):
            validation_scores["protenix_iptm"] = prot.get("iptm")
            validation_scores["protenix_ptm"] = prot.get("ptm")
            validation_scores["protenix_ranking_score"] = prot.get("ranking_score")
        elif prot and prot.get("error"):
            validation_scores["protenix_error"] = prot["error"][:200]

        # Pre-computed humanness from quick filter
        if "_oasis_score_mean" in c:
            scored["_oasis_score_mean"] = c["_oasis_score_mean"]
        if "_oasis_score_vh" in c:
            scored["_oasis_score_vh"] = c["_oasis_score_vh"]
        if "_oasis_score_vl" in c:
            scored["_oasis_score_vl"] = c["_oasis_score_vl"]

        scored["validation_scores"] = validation_scores
        scored_candidates.append(scored)

    # ---- 5. Save results ----
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Strip CIF strings from structure predictions to keep file manageable
    # (CIF files are already saved to disk)
    for c in scored_candidates:
        for key in ["structure_prediction", "structure_prediction_3chain"]:
            sp = c.get(key)
            if sp and "cif_string" in sp:
                del sp["cif_string"]

    with open(output_path, "w") as f:
        json.dump({
            "candidates": scored_candidates,
            "count": len(scored_candidates),
            "pre_filter_stats": {
                "total_input": len(candidates),
                "passed_hard_filter": len(survivors),
                "rejected": len(candidates) - len(survivors),
            },
            "scoring_summary": {
                "proteinmpnn_scored": sum(1 for c in scored_candidates
                                          if c.get("validation_scores", {}).get("proteinmpnn_ll_scfv") is not None),
                "antifold_scored": sum(1 for c in scored_candidates
                                       if c.get("validation_scores", {}).get("antifold_ll_scfv") is not None),
                "protenix_predicted": sum(1 for c in scored_candidates
                                          if c.get("validation_scores", {}).get("protenix_iptm") is not None),
                "fab_3chain_scored": sum(1 for c in scored_candidates
                                         if c.get("validation_scores", {}).get("proteinmpnn_ll_3chain") is not None),
            },
        }, f, indent=2)

    print(f"\nScored {len(scored_candidates)} candidates")
    print(f"Results saved to: {output_path}")

    print("\n" + "=" * 60)
    print("Scoring complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
