#!/usr/bin/env python3
"""Run the complete CD3 binder design pipeline.

This script orchestrates all pipeline steps:
0. Calibration (optional)
1. Setup targets
2. De novo design
3. Optimization
4. Structure prediction
4a. Score candidates (ProteinMPNN + AntiFold + Protenix)
4b. Post-hoc humanization (HuDiff)
5. Filtering
5b. Cross-validation
6. Bispecific formatting
7. Report generation

Usage:
    python scripts/run_full_pipeline.py [--config config.yaml] [--no-modal] [--skip-calibration]
"""

import argparse
from pathlib import Path
import subprocess
import sys


def run_step(script_name: str, args: list[str] = None) -> int:
    """Run a pipeline step script.

    Args:
        script_name: Name of script in scripts/ directory.
        args: Additional command line arguments.

    Returns:
        Return code from script.
    """
    script_path = Path("scripts") / script_name
    if not script_path.exists():
        print(f"ERROR: Script not found: {script_path}")
        return 1

    cmd = [sys.executable, str(script_path)]
    if args:
        cmd.extend(args)

    result = subprocess.run(cmd)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description="Run full CD3 binder design pipeline")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--no-modal", action="store_true", help="Run without Modal (local/mock)")
    parser.add_argument("--skip-calibration", action="store_true", help="Skip calibration step")
    parser.add_argument("--skip-setup", action="store_true", help="Skip target setup (if already done)")
    parser.add_argument("--start-from", type=int, default=0, help="Start from step N (0-9)")
    parser.add_argument("--stop-after", type=int, default=9, help="Stop after step N (0-9)")
    args = parser.parse_args()

    print("=" * 70)
    print("CD3 Binder Design Pipeline - Full Run")
    print("=" * 70)
    print(f"Config: {args.config}")
    print(f"Modal: {'disabled' if args.no_modal else 'enabled'}")
    print(f"Steps: {args.start_from} to {args.stop_after}")
    print("=" * 70)

    common_args = ["--config", args.config]
    if args.no_modal:
        common_args.append("--no-modal")

    steps = [
        ("00_run_calibration.py", "Calibration"),
        ("01_setup_targets.py", "Setup targets"),
        ("02_run_denovo_design.py", "De novo design"),
        ("03_run_optimization.py", "Optimization"),
        ("04_predict_structures.py", "Structure prediction"),
        ("04a_score_candidates.py", "Score candidates (ProteinMPNN + AntiFold + Protenix)"),
        ("04b_humanize_candidates.py", "Post-hoc humanization (HuDiff)"),
        ("05_filter_candidates.py", "Filtering"),
        ("05b_validate_candidates.py", "Cross-validation"),
        ("06_format_bispecifics.py", "Bispecific formatting"),
        ("07_generate_report.py", "Report generation"),
    ]

    for i, (script, description) in enumerate(steps):
        if i < args.start_from:
            continue
        if i > args.stop_after:
            break

        # Handle skip options
        if i == 0 and args.skip_calibration:
            print(f"\n[Step {i}] Skipping {description}")
            continue
        if i == 1 and args.skip_setup:
            print(f"\n[Step {i}] Skipping {description}")
            continue

        print(f"\n{'=' * 70}")
        print(f"[Step {i}] {description}")
        print("=" * 70)

        step_args = common_args.copy()

        # Add step-specific args
        if i == 0:  # Calibration
            if args.no_modal:
                step_args = ["--config", args.config, "--no-modal"]
            else:
                step_args = ["--config", args.config]

        ret = run_step(script, step_args)

        if ret != 0:
            print(f"\nERROR: Step {i} ({description}) failed with code {ret}")
            print("Pipeline aborted.")
            return ret

        print(f"\n[Step {i}] {description} - COMPLETE")

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE!")
    print("=" * 70)
    print("\nOutputs:")
    print("  - data/outputs/denovo/       : De novo designs")
    print("  - data/outputs/optimized/    : Optimized variants")
    print("  - data/outputs/structures/   : Structure predictions + scored candidates")
    print("  - data/outputs/filtered/     : Filtered candidates")
    print("  - data/outputs/validated/    : Cross-validated candidates")
    print("  - data/outputs/formatted/    : Bispecific sequences")
    print("  - data/outputs/reports/      : Final reports")
    print("\nNext steps:")
    print("  1. Review HTML report in data/outputs/reports/")
    print("  2. Select candidates for experimental validation")
    print("  3. Order gene synthesis for bispecific sequences")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    exit(main())
