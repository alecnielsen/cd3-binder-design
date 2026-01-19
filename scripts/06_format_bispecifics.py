#!/usr/bin/env python3
"""Step 6: Format candidates into bispecific antibody constructs.

Converts filtered candidates into multiple bispecific formats:
- CrossMab (Fab x Fab)
- Fab + scFv
- Fab + VHH
- IgG-scFv Morrison
- IgG-VHH Morrison

Usage:
    python scripts/06_format_bispecifics.py [--config config.yaml]
"""

import argparse
from pathlib import Path
import json


def main():
    parser = argparse.ArgumentParser(description="Format bispecifics")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input filtered candidates JSON")
    parser.add_argument("--output", type=str, default="data/outputs/formatted", help="Output dir")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Bispecific Formatting")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig
    from src.formatting import format_all, load_target_sequences, FORMATTERS

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

    candidates = data.get("candidates", data if isinstance(data, list) else [data])
    print(f"Loaded {len(candidates)} candidates")

    # Load target arm sequences from placeholder_targets.yaml based on config
    tumor_target_name = config.formatting.tumor_target
    try:
        target_vh, target_vl, target_display = load_target_sequences(tumor_target_name)
        print(f"\nTarget arm: {tumor_target_name} ({target_display})")
        print(f"  VH: {len(target_vh)} aa, VL: {len(target_vl)} aa")
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: Failed to load target sequences: {e}")
        return 1
    print(f"Formats: {config.formatting.formats}")

    # Format each candidate
    all_formatted = {}
    for fmt_name in config.formatting.formats:
        all_formatted[fmt_name] = []

    for c in candidates:
        candidate_id = c.get("candidate_id", "unknown")
        cd3_binder = c.get("sequence", "")
        cd3_binder_vl = c.get("sequence_vl")

        if not cd3_binder:
            continue

        try:
            constructs = format_all(
                target_vh=target_vh,
                target_vl=target_vl,
                cd3_binder=cd3_binder,
                cd3_binder_vl=cd3_binder_vl,
                name_prefix=candidate_id,
                target_name=target_display,
                formats=config.formatting.formats,
            )

            for fmt_name, construct in constructs.items():
                all_formatted[fmt_name].append({
                    "candidate_id": candidate_id,
                    "format": fmt_name,
                    "construct": construct.to_dict() if hasattr(construct, "to_dict") else {
                        "name": construct.name,
                        "chains": [{"name": ch.name, "sequence": ch.sequence} for ch in construct.chains],
                    },
                })

        except Exception as e:
            print(f"  Warning: Failed to format {candidate_id}: {e}")

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save all formats together
    all_output = output_dir / "all_bispecifics.json"
    with open(all_output, "w") as f:
        json.dump(all_formatted, f, indent=2)
    print(f"\nAll formats saved to: {all_output}")

    # Save individual format files
    for fmt_name, constructs in all_formatted.items():
        fmt_path = output_dir / f"{fmt_name}_constructs.json"
        with open(fmt_path, "w") as f:
            json.dump(constructs, f, indent=2)
        print(f"  {fmt_name}: {len(constructs)} constructs -> {fmt_path}")

    # Summary
    print(f"\nFormat summary:")
    for fmt_name, constructs in all_formatted.items():
        print(f"  {fmt_name}: {len(constructs)} constructs")

    print("\n" + "=" * 60)
    print("Bispecific formatting complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
