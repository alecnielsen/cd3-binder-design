#!/usr/bin/env python3
"""Step 4: Run structure prediction on candidates.

Runs Boltz-2 to predict binder-CD3 complex structures.

Usage:
    python scripts/04_predict_structures.py [--config config.yaml] [--no-modal]
"""

import argparse
from pathlib import Path
import json


def main():
    parser = argparse.ArgumentParser(description="Run structure prediction")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--input", type=str, help="Input candidates JSON (auto-detect if not specified)")
    parser.add_argument("--no-modal", action="store_true", help="Skip Modal (requires local GPU)")
    parser.add_argument("--output", type=str, default="data/outputs/structures", help="Output dir")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Structure Prediction")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig
    from src.structure.boltz_complex import Boltz2Predictor

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    # Find input candidates
    if args.input:
        input_paths = [Path(args.input)]
    else:
        # Auto-detect from outputs - aggregate ALL files from both directories
        output_dir = Path("data/outputs")
        denovo_files = list(output_dir.glob("denovo/*.json"))
        optimized_files = list(output_dir.glob("optimized/*.json"))
        input_paths = denovo_files + optimized_files

        if not input_paths:
            print("ERROR: No input candidate files found")
            print("Run steps 02 and 03 first, or specify --input")
            return 1

    print(f"\nFound {len(input_paths)} input file(s):")
    for p in input_paths:
        print(f"  - {p}")

    # Load and aggregate candidates from all input files
    candidates = []
    for input_path in input_paths:
        with open(input_path, "r") as f:
            data = json.load(f)

        if "designs" in data:
            file_candidates = data["designs"]
        elif "variants" in data:
            file_candidates = data["variants"]
        elif "vhh_designs" in data or "fab_designs" in data:
            # De novo design output: extract individual designs from both tracks
            file_candidates = []
            file_candidates.extend(data.get("vhh_designs", []))
            file_candidates.extend(data.get("fab_designs", []))
        else:
            file_candidates = data if isinstance(data, list) else [data]

        # Tag source file for traceability
        for c in file_candidates:
            if "source_file" not in c:
                c["source_file"] = str(input_path)

        candidates.extend(file_candidates)
        print(f"  Loaded {len(file_candidates)} candidates from {input_path.name}")

    print(f"\nTotal candidates: {len(candidates)}")

    # Get default target structure (fallback for candidates without specific target)
    default_target_pdb = None
    if config.design.target_structures:
        default_target_pdb = config.design.target_structures[0]
        if not Path(default_target_pdb).exists():
            print(f"WARNING: Default target not found: {default_target_pdb}")
            default_target_pdb = None

    if default_target_pdb:
        print(f"Default target: {default_target_pdb}")
    else:
        print("No default target - candidates must have target_structure field")

    # Get scFv linker from config
    scfv_linker = config.formatting.scfv_linker

    # Run predictions
    use_modal = not args.no_modal
    predictor = Boltz2Predictor(use_modal=use_modal)

    print(f"\nRunning Boltz-2 predictions (Modal={use_modal})...")
    print("This may take a while...")

    results = []
    for i, candidate in enumerate(candidates):
        # Get binder sequence - handle vh/vl pairs by creating scFv
        if "sequence" in candidate and candidate["sequence"]:
            binder_sequence = candidate["sequence"]
        elif "vh" in candidate:
            vh = candidate["vh"]
            vl = candidate.get("vl")
            if vl:
                # Create scFv: VH-linker-VL
                binder_sequence = vh + scfv_linker + vl
            else:
                binder_sequence = vh
        else:
            print(f"  Warning: No sequence found for candidate {i}")
            candidate["structure_prediction"] = None
            results.append(candidate)
            continue

        # Use per-candidate target_structure if available, otherwise fall back to default
        target_pdb = candidate.get("target_structure", default_target_pdb)
        if target_pdb and not Path(target_pdb).exists():
            print(f"  Warning: Target not found for candidate {i}: {target_pdb}")
            candidate["structure_prediction"] = None
            results.append(candidate)
            continue
        if not target_pdb:
            print(f"  Warning: No target structure for candidate {i}")
            candidate["structure_prediction"] = None
            results.append(candidate)
            continue

        try:
            result = predictor.predict_complex(
                binder_sequence=binder_sequence,
                target_pdb_path=target_pdb,
                seed=config.reproducibility.sampling_seed + i,
            )

            candidate["structure_prediction"] = result.to_dict()
            candidate["structure_prediction"]["binder_sequence_used"] = binder_sequence
            candidate["structure_prediction"]["target_structure"] = target_pdb
            results.append(candidate)

            if (i + 1) % 10 == 0:
                print(f"  Predicted {i + 1}/{len(candidates)}")

        except Exception as e:
            print(f"  Warning: Failed for candidate {i}: {e}")
            candidate["structure_prediction"] = None
            results.append(candidate)

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "candidates_with_structures.json"

    with open(output_path, "w") as f:
        json.dump({"candidates": results, "count": len(results)}, f, indent=2)

    print(f"\nPredicted structures for {len(results)} candidates")
    print(f"Results saved to: {output_path}")

    print("\n" + "=" * 60)
    print("Structure prediction complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
