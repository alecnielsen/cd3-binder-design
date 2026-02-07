#!/usr/bin/env python3
"""Step 5: Apply filtering cascade to candidates.

Filters candidates based on binding quality, humanness, liabilities, and developability.

Usage:
    python scripts/05_filter_candidates.py [--config config.yaml]
"""

import argparse
from pathlib import Path
import json


def main():
    parser = argparse.ArgumentParser(description="Filter candidates")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input candidates JSON")
    parser.add_argument("--output", type=str, default="data/outputs/filtered", help="Output dir")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Filtering")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig
    from src.pipeline.filter_cascade import CandidateScore, run_filter_cascade
    from src.pipeline.ranking import RankedCandidate, worst_metric_rank, diversity_select
    from src.analysis.liabilities import LiabilityScanner
    from src.analysis.humanness import score_humanness_pair
    from src.analysis.developability import DevelopabilityAssessor
    from src.structure.interface_analysis import InterfaceAnalyzer

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    # Find input
    if args.input:
        input_path = Path(args.input)
    else:
        structures_path = Path("data/outputs/structures/candidates_with_structures.json")
        if structures_path.exists():
            input_path = structures_path
        else:
            print("ERROR: No input found. Run step 04 first or specify --input")
            return 1

    print(f"\nInput: {input_path}")

    # Load candidates
    with open(input_path, "r") as f:
        data = json.load(f)

    candidates_data = data.get("candidates", data if isinstance(data, list) else [data])
    print(f"Loaded {len(candidates_data)} candidates")

    # Print threshold info
    thresholds = config.get_effective_thresholds()
    print(f"\nFilter thresholds:")
    print(f"  min_pdockq: {thresholds['min_pdockq']}")
    print(f"  min_interface_area: {thresholds['min_interface_area']}")
    print(f"  min_contacts: {thresholds['min_contacts']}")
    if config.calibrated_min_pdockq:
        print("  (Using calibrated thresholds)")

    # Convert to CandidateScore objects
    scanner = LiabilityScanner()
    dev_assessor = DevelopabilityAssessor()
    interface_analyzer = InterfaceAnalyzer(
        okt3_epitope_residues=config.epitope.okt3_epitope_residues
    )

    scored_candidates = []
    for c in candidates_data:
        # Extract sequences - handle both single sequence and vh/vl pairs
        vh_seq = c.get("sequence") or c.get("vh", "")
        vl_seq = c.get("sequence_vl") or c.get("vl")

        # Determine binder type - check if single sequence is actually an scFv
        binder_type = c.get("binder_type")
        if binder_type is None:
            if vl_seq is not None:
                binder_type = "scfv"
            else:
                # Check if the single sequence is a concatenated scFv
                from src.utils.constants import parse_scfv, is_likely_scfv
                if is_likely_scfv(vh_seq):
                    parsed = parse_scfv(vh_seq)
                    if parsed:
                        vh_seq, vl_seq = parsed
                        binder_type = "scfv"
                    else:
                        binder_type = "vhh"
                else:
                    binder_type = "vhh"

        score = CandidateScore(
            candidate_id=c.get("design_id", c.get("name", "unknown")),
            sequence=vh_seq,
            sequence_vl=vl_seq,
            binder_type=binder_type,
            source=c.get("source", "unknown"),
        )

        # Add structure prediction metrics
        sp = c.get("structure_prediction", {})
        if sp:
            score.iptm = sp.get("iptm") or sp.get("ipTM") or 0.0
            score.pdockq = sp.get("pdockq")
            score.interface_area = sp.get("interface_area")
            score.num_contacts = sp.get("num_contacts")

            # Epitope annotation
            if sp.get("interface_residues_target"):
                epitope_class, overlap = interface_analyzer.annotate_epitope_class(
                    sp["interface_residues_target"]
                )
                score.epitope_class = epitope_class
                score.okt3_overlap = overlap

        # Liability analysis - scan with CDR detection for accurate filtering
        try:
            # Scan VH with CDR detection
            vh_report = scanner.scan_with_cdr_detection(vh_seq, chain_type="H")

            # Scan VL if present
            if vl_seq:
                vl_report = scanner.scan_with_cdr_detection(vl_seq, chain_type="L")
                # Combine reports - offset VL positions by VH length
                vh_len = len(vh_seq)
                from src.analysis.liabilities import LiabilitySite
                all_deamidation = vh_report.deamidation_sites + [
                    LiabilitySite(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                    for s in vl_report.deamidation_sites
                ]
                all_isomerization = vh_report.isomerization_sites + [
                    LiabilitySite(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                    for s in vl_report.isomerization_sites
                ]
                all_glycosylation = vh_report.glycosylation_sites + [
                    LiabilitySite(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                    for s in vl_report.glycosylation_sites
                ]
                all_oxidation = vh_report.oxidation_sites + [
                    LiabilitySite(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                    for s in vl_report.oxidation_sites
                ]
            else:
                all_deamidation = vh_report.deamidation_sites
                all_isomerization = vh_report.isomerization_sites
                all_glycosylation = vh_report.glycosylation_sites
                all_oxidation = vh_report.oxidation_sites

            # Extract positions as integers for JSON serialization
            score.deamidation_sites = [s.position for s in all_deamidation]
            score.isomerization_sites = [s.position for s in all_isomerization]
            score.glycosylation_sites = [s.position for s in all_glycosylation]
            score.oxidation_sites = [s.position for s in all_oxidation]
            score.unpaired_cys = vh_report.unpaired_cysteines + (vl_report.unpaired_cysteines if vl_seq else 0)

            # Count CDR-specific liabilities for filtering
            score.cdr_deamidation_count = sum(1 for s in all_deamidation if s.in_cdr)
            score.cdr_isomerization_count = sum(1 for s in all_isomerization if s.in_cdr)
            score.cdr_glycosylation_count = sum(1 for s in all_glycosylation if s.in_cdr)
            score.cdr_oxidation_count = sum(1 for s in all_oxidation if s.in_cdr)
        except Exception:
            pass

        # Humanness scoring
        try:
            humanness_report = score_humanness_pair(vh_seq, vl_seq)
            score.oasis_score_vh = humanness_report.vh_report.oasis_score
            score.oasis_score_vl = humanness_report.vl_report.oasis_score if humanness_report.vl_report else None
            score.oasis_score_mean = humanness_report.mean_score
        except Exception:
            pass

        # Developability
        try:
            dev_report = dev_assessor.assess(vh_seq, vl_seq, include_humanness=False)
            score.cdr_h3_length = dev_report.cdr_h3_length
            score.net_charge = dev_report.physicochemical.net_charge
            score.isoelectric_point = dev_report.physicochemical.isoelectric_point
            score.hydrophobic_patches = dev_report.aggregation.hydrophobic_patches
        except Exception:
            pass

        scored_candidates.append(score)

    # Run filtering
    print(f"\nRunning filter cascade...")
    filtered, stats = run_filter_cascade(
        candidates=scored_candidates,
        config=config,
        min_candidates=config.filtering.min_candidates,
    )

    print(f"\nFilter statistics:")
    print(f"  Input: {stats['total_input']}")
    print(f"  Passing (first pass): {stats['passing_first_pass']}")
    print(f"  Final passing: {stats['final_passing']}")
    if stats.get("used_fallback"):
        print(f"  Fallback relaxations: {len(stats['relaxations_applied'])}")

    # --- Ranking ---
    ranking_method = config.ranking.method
    print(f"\nRanking method: {ranking_method}")

    # Build a lookup from candidate_id to raw data for ptm/plddt extraction
    raw_lookup = {}
    for c in candidates_data:
        cid = c.get("design_id", c.get("name", "unknown"))
        raw_lookup[cid] = c

    if ranking_method == "worst_metric_rank" and filtered:
        # Build RankedCandidate objects
        ranked_candidates = []
        for score in filtered:
            raw = raw_lookup.get(score.candidate_id, {})
            sp = raw.get("structure_prediction", {})

            rc = RankedCandidate(
                candidate_id=score.candidate_id,
                sequence=score.sequence,
                sequence_vl=score.sequence_vl,
                iptm=score.iptm or 0.0,
                ptm=sp.get("ptm", 0.0),
                plddt=sp.get("plddt_mean", 0.0),
                interface_area=score.interface_area or 0.0,
                num_contacts=score.num_contacts or 0,
                humanness=score.oasis_score_mean or score.oasis_score_vh or 0.0,
                _score_ref=score,
            )
            ranked_candidates.append(rc)

        # Apply worst-metric-rank
        worst_metric_rank(ranked_candidates, config.ranking.metric_weights)

        # Apply diversity selection if enabled
        n_final = config.output.num_final_candidates
        if config.ranking.use_diversity_selection:
            selected = diversity_select(
                ranked_candidates, n_final, alpha=config.ranking.diversity_alpha
            )
            print(f"  Diversity selection: {len(ranked_candidates)} -> {len(selected)}")
        else:
            selected = ranked_candidates[:n_final]

        # Map back to CandidateScore, updating ranks
        filtered = []
        for i, rc in enumerate(selected):
            score = rc._score_ref
            score.rank = i + 1
            score.composite_score = rc.quality_key  # Store quality_key as score
            filtered.append(score)

        stats["ranking_method"] = "worst_metric_rank"
        stats["diversity_selection"] = config.ranking.use_diversity_selection
    else:
        # Legacy composite score ranking (already sorted by run_filter_cascade)
        n_final = config.output.num_final_candidates
        filtered = filtered[:n_final]
        stats["ranking_method"] = "composite"

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "filtered_candidates.json"

    with open(output_path, "w") as f:
        json.dump({
            "candidates": [c.to_dict() for c in filtered],
            "count": len(filtered),
            "filter_stats": stats,
        }, f, indent=2)

    print(f"\nFiltered to {len(filtered)} candidates")
    print(f"Results saved to: {output_path}")

    # Print top candidates
    print(f"\nTop candidates:")
    for c in filtered[:5]:
        iptm_str = f"{c.iptm:.3f}" if c.iptm else "N/A"
        area_str = f"{c.interface_area:.0f}" if c.interface_area else "N/A"
        humanness_str = f"{c.oasis_score_mean:.3f}" if c.oasis_score_mean else "N/A"
        if ranking_method == "worst_metric_rank":
            print(f"  {c.rank}. {c.candidate_id}: quality_key={c.composite_score:.1f}, "
                  f"ipTM={iptm_str}, area={area_str}, humanness={humanness_str}")
        else:
            print(f"  {c.rank}. {c.candidate_id}: score={c.composite_score:.3f}, "
                  f"ipTM={iptm_str}")

    print("\n" + "=" * 60)
    print("Filtering complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
