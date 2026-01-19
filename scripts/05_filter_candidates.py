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
    from src.analysis.liabilities import LiabilityScanner
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
    interface_analyzer = InterfaceAnalyzer()

    scored_candidates = []
    for c in candidates_data:
        # Extract sequences - handle both single sequence and vh/vl pairs
        vh_seq = c.get("sequence") or c.get("vh", "")
        vl_seq = c.get("sequence_vl") or c.get("vl")

        # Determine binder type
        binder_type = c.get("binder_type")
        if binder_type is None:
            binder_type = "vhh" if vl_seq is None else "scfv"

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

        # Liability analysis - scan full sequence
        try:
            full_sequence = vh_seq + (vl_seq or "")
            liability_report = scanner.scan(full_sequence)
            score.deamidation_sites = liability_report.deamidation_sites
            score.isomerization_sites = liability_report.isomerization_sites
            score.glycosylation_sites = liability_report.glycosylation_sites
            score.oxidation_sites = liability_report.oxidation_sites
        except Exception:
            pass

        # Developability
        try:
            dev_report = dev_assessor.assess(vh_seq, vl_seq, include_humanness=False)
            score.cdr_h3_length = dev_report.cdr_h3_length
            score.net_charge = dev_report.physicochemical.net_charge
            score.isoelectric_point = dev_report.physicochemical.isoelectric_point
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
    print(f"\nTop 5 candidates:")
    for c in filtered[:5]:
        pdockq_str = f"{c.pdockq:.3f}" if c.pdockq is not None else "N/A"
        print(f"  {c.rank}. {c.candidate_id}: score={c.composite_score:.3f}, pDockQ={pdockq_str}")

    print("\n" + "=" * 60)
    print("Filtering complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
