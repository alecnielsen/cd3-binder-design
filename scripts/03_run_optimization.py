#!/usr/bin/env python3
"""Step 3: Run optimization of existing binders.

Generates humanized and affinity variants of known CD3 binders.

Usage:
    python scripts/03_run_optimization.py [--config config.yaml]
"""

import argparse
from pathlib import Path
import json


def main():
    parser = argparse.ArgumentParser(description="Run binder optimization")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--output", type=str, default="data/outputs/optimized", help="Output dir")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Optimization")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig
    from src.design.optimization import optimize_existing_binders
    from src.design.affinity_variants import AffinityMutationLibrary, AffinityVariantGenerator

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    print(f"\nStarting sequences: {config.design.starting_sequences}")
    print(f"Affinity variants: {config.design.affinity_variants}")

    # Generate optimized variants
    print("\nGenerating optimized variants...")
    variants = optimize_existing_binders(
        binder_names=config.design.starting_sequences,
        include_humanization=True,
        include_back_mutations=config.filtering.generate_back_mutations,
    )

    # Generate affinity variants
    print("\nGenerating affinity variants...")
    library = AffinityMutationLibrary()
    generator = AffinityVariantGenerator(library)

    all_variants = []
    for v in variants:
        all_variants.append(v.to_dict())

        affinity_panel = generator.generate_affinity_panel(
            parent_name=v.name,
            vh=v.vh,
            vl=v.vl,
            target_classes=config.design.affinity_variants,
        )
        for av in affinity_panel:
            all_variants.append(av.to_dict())

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "optimized_variants.json"

    with open(output_path, "w") as f:
        json.dump({"variants": all_variants, "count": len(all_variants)}, f, indent=2)

    print(f"\nGenerated {len(all_variants)} total variants")
    print(f"Results saved to: {output_path}")

    print("\n" + "=" * 60)
    print("Optimization complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
