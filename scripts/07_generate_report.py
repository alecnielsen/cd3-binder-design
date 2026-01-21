#!/usr/bin/env python3
"""Step 7: Generate final report and scorecards.

Generates HTML and JSON reports for filtered candidates.

Usage:
    python scripts/07_generate_report.py [--config config.yaml]
"""

import argparse
from pathlib import Path
import json


def main():
    parser = argparse.ArgumentParser(description="Generate report")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input filtered candidates JSON")
    parser.add_argument("--output", type=str, default="data/outputs/reports", help="Output dir")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Report Generation")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig, get_provenance
    from src.pipeline.filter_cascade import CandidateScore, FilterResult
    from src.pipeline.report_generator import generate_report

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

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
    with open(input_path, "r") as f:
        data = json.load(f)

    candidates_data = data.get("candidates", [])
    filter_stats = data.get("filter_stats", {})

    print(f"Loaded {len(candidates_data)} candidates")

    # Convert to CandidateScore objects
    candidates = []
    for c in candidates_data:
        score = CandidateScore(
            candidate_id=c.get("candidate_id", "unknown"),
            sequence=c.get("sequence", ""),
            binder_type=c.get("binder_type", "vhh"),
            source=c.get("source", "unknown"),
            composite_score=c.get("composite_score", 0.0),
            rank=c.get("rank", 0),
        )

        # Copy metrics
        if "binding" in c:
            score.pdockq = c["binding"].get("pdockq")
            score.interface_area = c["binding"].get("interface_area")
            score.num_contacts = c["binding"].get("num_contacts")

        if "humanness" in c:
            score.oasis_score_vh = c["humanness"].get("oasis_score_vh")
            score.oasis_score_mean = c["humanness"].get("oasis_score_mean")

        if "epitope" in c:
            score.epitope_class = c["epitope"].get("epitope_class", "unknown")
            score.okt3_overlap = c["epitope"].get("okt3_overlap", 0.0)

        if "liabilities" in c:
            score.deamidation_sites = c["liabilities"].get("deamidation_sites", [])
            score.glycosylation_sites = c["liabilities"].get("glycosylation_sites", [])
            score.oxidation_sites = c["liabilities"].get("oxidation_sites", [])

        # Restore filter results for report display
        if "filter_results" in c:
            for filter_name, result_value in c["filter_results"].items():
                score.filter_results[filter_name] = FilterResult(result_value)

        score.risk_flags = c.get("risk_flags", [])

        candidates.append(score)

    # Select final candidates
    num_final = config.output.num_final_candidates
    final_candidates = candidates[:num_final]
    print(f"Selecting top {num_final} for report")

    # Get provenance
    provenance = get_provenance()
    provenance["config_hash"] = config.config_hash()

    # Generate reports
    print(f"\nGenerating reports...")
    saved_files = generate_report(
        candidates=final_candidates,
        filter_stats=filter_stats,
        output_dir=args.output,
        provenance=provenance,
    )

    print(f"\n" + "=" * 60)
    print("Report generation complete!")
    print(f"Reports saved to: {args.output}")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
