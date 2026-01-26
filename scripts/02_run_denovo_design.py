#!/usr/bin/env python3
"""Step 2: Run de novo design generation.

Runs BoltzGen on Modal to generate VHH and Fab binders.

VHH: De novo nanobody design (~120 aa)
Fab: CDR redesign using human antibody scaffolds (VH ~120 aa + VL ~107 aa)

Usage:
    python scripts/02_run_denovo_design.py [--config config.yaml] [--no-modal]

Prerequisites:
    1. Run scripts/setup_fab_scaffolds.py to download scaffold files
    2. Deploy Modal app: modal deploy modal/boltzgen_app.py
"""

import argparse
from pathlib import Path
import json


def main():
    parser = argparse.ArgumentParser(description="Run de novo design generation")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--no-modal", action="store_true", help="Run mock locally")
    parser.add_argument("--output", type=str, default="data/outputs/denovo", help="Output dir")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - De Novo Design")
    print("=" * 60)

    from src.pipeline.config import PipelineConfig
    from src.design.denovo_design import run_denovo_design

    # Load config
    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
    else:
        config = PipelineConfig()

    if not config.design.target_structures:
        print("ERROR: No target structures in config")
        return 1

    print(f"\nTarget structures: {config.design.target_structures}")
    print(f"VHH designs: {config.design.num_vhh_designs}")
    print(f"Fab designs: {config.design.num_fab_designs}")
    print(f"Fab scaffolds: {config.design.fab_scaffolds}")
    print(f"Seed: {config.reproducibility.boltzgen_seed}")

    # Check scaffold files exist
    scaffold_dir = Path(config.design.fab_scaffold_dir)
    if config.design.num_fab_designs > 0:
        if not scaffold_dir.exists():
            print(f"\nWARNING: Scaffold directory not found: {scaffold_dir}")
            print("Run 'python scripts/setup_fab_scaffolds.py' first.")
            if not args.no_modal:
                return 1

    # Run design
    use_modal = not args.no_modal
    print(f"\nRunning BoltzGen (Modal={use_modal})...")

    result = run_denovo_design(
        target_structures=config.design.target_structures,
        num_vhh=config.design.num_vhh_designs,
        num_fab=config.design.num_fab_designs,
        fab_scaffolds=config.design.fab_scaffolds,
        fab_scaffold_dir=config.design.fab_scaffold_dir,
        seed=config.reproducibility.boltzgen_seed,
        output_dir=args.output,
        use_modal=use_modal,
    )

    print(f"\nGenerated {result.total_designs} designs")
    print(f"  VHH: {len(result.vhh_designs)}")
    print(f"  Fab: {len(result.fab_designs)}")

    # Save results
    output_path = result.save()
    print(f"\nResults saved to: {output_path}")

    print("\n" + "=" * 60)
    print("De novo design complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
