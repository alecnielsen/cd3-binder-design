#!/usr/bin/env python3
"""Step 0: Run calibration to set filter thresholds.

This script runs Boltz-2 on known binders (teplizumab, SP34, UCHT1)
to establish appropriate filter thresholds before running the main pipeline.

IMPORTANT: Run this BEFORE filtering to ensure thresholds are appropriate.

Usage:
    python scripts/00_run_calibration.py [--config config.yaml] [--no-modal]
"""

import argparse
from pathlib import Path
import json
import yaml


def main():
    parser = argparse.ArgumentParser(description="Run calibration for filter thresholds")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file path")
    parser.add_argument("--no-modal", action="store_true", help="Run locally without Modal")
    parser.add_argument("--output", type=str, default="data/outputs/calibration.json", help="Output file")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Calibration")
    print("=" * 60)

    # Load config
    from src.pipeline.config import PipelineConfig

    if Path(args.config).exists():
        config = PipelineConfig.load(args.config)
        print(f"Loaded config from {args.config}")
    else:
        config = PipelineConfig()
        print("Using default config")

    # Check target structures exist
    if not config.design.target_structures:
        print("ERROR: No target structures specified in config")
        print("Add target_structures to design section of config.yaml")
        return 1

    target_pdb = config.design.target_structures[0]
    if not Path(target_pdb).exists():
        print(f"ERROR: Target structure not found: {target_pdb}")
        return 1

    print(f"\nTarget structure: {target_pdb}")
    print(f"Known binders: {config.calibration.positive_controls}")

    # Load known binder sequences
    from src.design.optimization import SequenceOptimizer

    optimizer = SequenceOptimizer()
    known_sequences = []
    scfv_linker = config.formatting.scfv_linker

    for name in config.calibration.positive_controls:
        try:
            seq = optimizer.load_starting_sequence(name)
            # For paired antibodies (VH+VL), construct scFv for accurate calibration
            if seq.vl:
                binder_seq = seq.vh + scfv_linker + seq.vl
                print(f"  Loaded {name}: {len(seq.vh)} aa VH + {len(seq.vl)} aa VL = {len(binder_seq)} aa scFv")
            else:
                binder_seq = seq.vh
                print(f"  Loaded {name}: {len(binder_seq)} aa VHH")
            known_sequences.append(binder_seq)
        except Exception as e:
            print(f"  Warning: Could not load {name}: {e}")

    if not known_sequences:
        print("ERROR: No known binder sequences could be loaded")
        return 1

    # Run calibration
    print(f"\nRunning Boltz-2 on {len(known_sequences)} known binders...")
    print("  (This may take several minutes on Modal)")

    use_modal = not args.no_modal

    from src.structure.boltz_complex import run_calibration

    try:
        calibration_results = run_calibration(
            known_binder_sequences=known_sequences,
            target_pdb_path=target_pdb,
            use_modal=use_modal,
            control_names=config.calibration.positive_controls,
        )
    except Exception as e:
        print(f"ERROR: Calibration failed: {e}")
        return 1

    # Save CIF files for controls
    cif_strings = calibration_results.pop("cif_strings", {})
    if cif_strings:
        cif_dir = Path("data/outputs/calibration_cif")
        cif_dir.mkdir(parents=True, exist_ok=True)
        for name, cif_content in cif_strings.items():
            cif_path = cif_dir / f"{name}.cif"
            cif_path.write_text(cif_content)
            print(f"  Saved CIF: {cif_path}")

    # Print results
    print("\n" + "=" * 60)
    print("Calibration Results")
    print("=" * 60)

    stats = calibration_results["known_binder_stats"]
    print(f"\nKnown binder statistics:")
    print(f"  pDockQ: min={stats['pdockq']['min']:.3f}, max={stats['pdockq']['max']:.3f}, mean={stats['pdockq']['mean']:.3f}")
    print(f"  Interface area: min={stats['interface_area']['min']:.1f}, max={stats['interface_area']['max']:.1f}")
    print(f"  Contacts: min={stats['contacts']['min']}, max={stats['contacts']['max']}")

    thresholds = calibration_results["calibrated_thresholds"]
    print(f"\nCalibrated thresholds (with margin):")
    print(f"  min_pdockq: {thresholds['min_pdockq']:.3f}")
    print(f"  min_interface_area: {thresholds['min_interface_area']:.1f}")
    print(f"  min_contacts: {thresholds['min_contacts']}")

    # --- Validation baselines (ProteinMPNN, AntiFold, Protenix on controls) ---
    if config.calibration.run_validation_baselines and cif_strings:
        print(f"\n--- Validation Baselines ---")
        validation_baselines = {}

        # Affinity scoring on controls
        try:
            from src.analysis.affinity_scoring import batch_score_affinity
            cif_dir_str = str(Path("data/outputs/calibration_cif"))
            control_names = list(cif_strings.keys())
            control_types = ["scfv"] * len(control_names)  # Controls are scFvs

            aff_results = batch_score_affinity(
                cif_dir=cif_dir_str,
                design_ids=control_names,
                binder_types=control_types,
                run_proteinmpnn=config.validation.run_proteinmpnn,
                run_antifold=config.validation.run_antifold,
            )

            mpnn_scores = {}
            af_scores = {}
            for r in aff_results:
                if r.proteinmpnn_ll is not None:
                    mpnn_scores[r.design_id] = r.proteinmpnn_ll
                    print(f"  {r.design_id} ProteinMPNN: {r.proteinmpnn_ll:.3f}")
                if r.antifold_ll is not None:
                    af_scores[r.design_id] = r.antifold_ll
                    print(f"  {r.design_id} AntiFold: {r.antifold_ll:.3f}")

            if mpnn_scores:
                vals = list(mpnn_scores.values())
                validation_baselines["proteinmpnn_ll"] = {
                    "min": min(vals), "max": max(vals),
                    "mean": sum(vals) / len(vals),
                    "by_control": mpnn_scores,
                }
            if af_scores:
                vals = list(af_scores.values())
                validation_baselines["antifold_ll"] = {
                    "min": min(vals), "max": max(vals),
                    "mean": sum(vals) / len(vals),
                    "by_control": af_scores,
                }
        except Exception as e:
            print(f"  Affinity scoring failed: {e}")

        # Protenix on controls
        if use_modal and config.validation.run_protenix:
            try:
                import modal
                from src.structure.pdb_utils import extract_sequence_from_pdb
                target_sequence = extract_sequence_from_pdb(target_pdb)
                protenix_fn = modal.Function.from_name("protenix-cd3", "predict_complex")

                protenix_scores = {}
                for i, name in enumerate(config.calibration.positive_controls):
                    seq = known_sequences[i]
                    print(f"  Running Protenix on {name}...")
                    try:
                        pred = protenix_fn.remote(
                            binder_sequence=seq,
                            target_sequence=target_sequence,
                            binder_type="scfv",
                            seed=config.validation.protenix_seeds[0],
                        )
                        if pred.get("iptm") is not None:
                            protenix_scores[name] = pred["iptm"]
                            print(f"    ipTM={pred['iptm']}")
                    except Exception as e:
                        print(f"    Error: {e}")

                if protenix_scores:
                    vals = list(protenix_scores.values())
                    validation_baselines["protenix_iptm"] = {
                        "min": min(vals), "max": max(vals),
                        "mean": sum(vals) / len(vals),
                        "by_control": protenix_scores,
                    }
            except Exception as e:
                print(f"  Protenix scoring failed: {e}")

        if validation_baselines:
            calibration_results["validation_baselines"] = validation_baselines
            print(f"\n  Validation baselines computed for {len(validation_baselines)} metrics")

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(calibration_results, f, indent=2)
    print(f"\nCalibration results saved to: {output_path}")

    # Update config with calibrated thresholds
    config.calibrated_min_pdockq = thresholds["min_pdockq"]
    config.calibrated_min_interface_area = thresholds["min_interface_area"]
    config.calibrated_min_contacts = thresholds["min_contacts"]
    config.save(args.config)
    print(f"Config updated with calibrated thresholds: {args.config}")

    print("\n" + "=" * 60)
    print("Calibration complete!")
    print("You can now run the design pipeline with calibrated thresholds.")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
