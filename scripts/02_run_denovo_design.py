#!/usr/bin/env python3
"""Step 2: Run de novo design generation.

Runs BoltzGen on Modal to generate VHH and scFv binders.

Usage:
    python scripts/02_run_denovo_design.py [--config config.yaml] [--no-modal]
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
    print(f"scFv designs: {config.design.num_scfv_designs}")
    print(f"Seed: {config.reproducibility.boltzgen_seed}")

    # Run design
    use_modal = not args.no_modal
    print(f"\nRunning BoltzGen (Modal={use_modal})...")

    result = run_denovo_design(
        target_structures=config.design.target_structures,
        num_vhh=config.design.num_vhh_designs,
        num_scfv=config.design.num_scfv_designs,
        seed=config.reproducibility.boltzgen_seed,
        output_dir=args.output,
        use_modal=use_modal,
    )

    print(f"\nGenerated {result.total_designs} designs")
    print(f"  VHH: {len(result.vhh_designs)}")
    print(f"  scFv: {len(result.scfv_designs)}")

    # Save results
    output_path = result.save()
    print(f"\nResults saved to: {output_path}")

    print("\n" + "=" * 60)
    print("De novo design complete!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
