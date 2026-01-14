"""End-to-end pipeline orchestration for CD3 binder design.

This module coordinates all pipeline stages:
1. Calibration (optional, run once)
2. Design generation (de novo + optimization)
3. Structure prediction
4. Filtering cascade
5. Candidate selection and ranking
6. Bispecific formatting
7. Report generation
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import json
import datetime

from src.pipeline.config import PipelineConfig, get_provenance
from src.pipeline.filter_cascade import FilterCascade, CandidateScore, run_filter_cascade


@dataclass
class PipelineResult:
    """Complete result from pipeline execution."""

    # Candidates at each stage
    denovo_candidates: list[dict] = field(default_factory=list)
    optimized_candidates: list[dict] = field(default_factory=list)
    all_candidates: list[CandidateScore] = field(default_factory=list)
    filtered_candidates: list[CandidateScore] = field(default_factory=list)
    final_candidates: list[CandidateScore] = field(default_factory=list)

    # Statistics
    filter_stats: dict = field(default_factory=dict)
    calibration_results: Optional[dict] = None

    # Metadata
    config: Optional[PipelineConfig] = None
    provenance: dict = field(default_factory=dict)
    timestamp: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "num_denovo_candidates": len(self.denovo_candidates),
            "num_optimized_candidates": len(self.optimized_candidates),
            "num_total_candidates": len(self.all_candidates),
            "num_filtered_candidates": len(self.filtered_candidates),
            "num_final_candidates": len(self.final_candidates),
            "filter_stats": self.filter_stats,
            "calibration_results": self.calibration_results,
            "final_candidates": [c.to_dict() for c in self.final_candidates],
            "provenance": self.provenance,
            "timestamp": self.timestamp,
        }

    def save(self, output_path: str) -> str:
        """Save results to JSON file."""
        with open(output_path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        return output_path


class DesignPipeline:
    """End-to-end pipeline for CD3 binder design.

    This class orchestrates the complete workflow from
    design generation through final candidate selection.
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        """Initialize pipeline.

        Args:
            config: Pipeline configuration.
        """
        self.config = config or PipelineConfig()
        self.result = PipelineResult(config=self.config)

    def run_calibration(self, use_modal: bool = True) -> dict:
        """Run calibration to set filter thresholds.

        Uses known binders to establish appropriate thresholds.

        Args:
            use_modal: If True, use Modal for GPU compute.

        Returns:
            Calibration results dictionary.
        """
        print("Running calibration with known binders...")

        # Load known binder sequences
        from src.design.optimization import SequenceOptimizer

        optimizer = SequenceOptimizer()
        known_sequences = []

        for name in self.config.calibration.positive_controls:
            try:
                seq = optimizer.load_starting_sequence(name)
                known_sequences.append(seq.vh)
                print(f"  Loaded {name}")
            except Exception as e:
                print(f"  Warning: Could not load {name}: {e}")

        if not known_sequences:
            raise RuntimeError("No known binder sequences could be loaded for calibration")

        # Run Boltz-2 on known binders
        from src.structure.boltz_complex import run_calibration

        target_pdb = self.config.design.target_structures[0] if self.config.design.target_structures else None

        if target_pdb is None:
            raise RuntimeError("No target structure specified for calibration")

        calibration_results = run_calibration(
            known_binder_sequences=known_sequences,
            target_pdb_path=target_pdb,
            use_modal=use_modal,
        )

        # Update config with calibrated thresholds
        thresholds = calibration_results["calibrated_thresholds"]
        self.config.calibrated_min_pdockq = thresholds["min_pdockq"]
        self.config.calibrated_min_interface_area = thresholds["min_interface_area"]
        self.config.calibrated_min_contacts = thresholds["min_contacts"]

        self.result.calibration_results = calibration_results

        print(f"Calibration complete. Thresholds:")
        print(f"  min_pdockq: {thresholds['min_pdockq']:.3f}")
        print(f"  min_interface_area: {thresholds['min_interface_area']:.1f}")
        print(f"  min_contacts: {thresholds['min_contacts']}")

        return calibration_results

    def run_denovo_design(self, use_modal: bool = True) -> list[dict]:
        """Run de novo design generation.

        Args:
            use_modal: If True, use Modal for GPU compute.

        Returns:
            List of design dictionaries.
        """
        print("Running de novo design generation...")

        from src.design.denovo_design import run_denovo_design

        result = run_denovo_design(
            target_structures=self.config.design.target_structures,
            num_vhh=self.config.design.num_vhh_designs,
            num_scfv=self.config.design.num_scfv_designs,
            seed=self.config.reproducibility.boltzgen_seed,
            output_dir=self.config.output.output_dir,
            use_modal=use_modal,
        )

        # Convert to dicts
        designs = []
        for d in result.vhh_designs + result.scfv_designs:
            designs.append(d.to_dict())

        self.result.denovo_candidates = designs
        print(f"Generated {len(designs)} de novo designs")

        return designs

    def run_optimization(self) -> list[dict]:
        """Run optimization of existing binders.

        Returns:
            List of optimized variant dictionaries.
        """
        print("Running optimization of existing binders...")

        from src.design.optimization import optimize_existing_binders
        from src.design.affinity_variants import generate_affinity_variants, AffinityMutationLibrary, AffinityVariantGenerator

        # Generate humanized variants
        variants = optimize_existing_binders(
            binder_names=self.config.design.starting_sequences,
            include_humanization=True,
            include_back_mutations=self.config.filtering.generate_back_mutations,
        )

        # Generate affinity variants for each
        all_variants = []
        for v in variants:
            all_variants.append(v.to_dict())

            # Generate affinity panel
            library = AffinityMutationLibrary()
            generator = AffinityVariantGenerator(library)
            affinity_panel = generator.generate_affinity_panel(
                parent_name=v.name,
                vh=v.vh,
                vl=v.vl,
                target_classes=self.config.design.affinity_variants,
            )
            for av in affinity_panel:
                all_variants.append(av.to_dict())

        self.result.optimized_candidates = all_variants
        print(f"Generated {len(all_variants)} optimized variants")

        return all_variants

    def run_structure_prediction(
        self,
        candidates: list[dict],
        use_modal: bool = True,
    ) -> list[dict]:
        """Run structure prediction on candidates.

        Args:
            candidates: List of candidate dictionaries.
            use_modal: If True, use Modal for GPU compute.

        Returns:
            Updated candidates with structure predictions.
        """
        print(f"Running structure prediction on {len(candidates)} candidates...")

        from src.structure.boltz_complex import Boltz2Predictor

        predictor = Boltz2Predictor(use_modal=use_modal)
        target_pdb = self.config.design.target_structures[0]

        for i, candidate in enumerate(candidates):
            try:
                result = predictor.predict_complex(
                    binder_sequence=candidate["sequence"],
                    target_pdb_path=target_pdb,
                    seed=self.config.reproducibility.sampling_seed + i,
                )

                candidate["structure_prediction"] = {
                    "pdockq": result.pdockq,
                    "ptm": result.ptm,
                    "plddt_mean": result.plddt_mean,
                    "interface_area": result.interface_area,
                    "num_contacts": result.num_contacts,
                    "interface_residues_target": result.interface_residues_target,
                }

                if (i + 1) % 10 == 0:
                    print(f"  Predicted {i + 1}/{len(candidates)}")

            except Exception as e:
                print(f"  Warning: Structure prediction failed for candidate {i}: {e}")
                candidate["structure_prediction"] = None

        return candidates

    def run_analysis(self, candidates: list[dict]) -> list[CandidateScore]:
        """Run analysis on candidates to populate scores.

        Args:
            candidates: List of candidate dictionaries.

        Returns:
            List of CandidateScore objects.
        """
        print(f"Running analysis on {len(candidates)} candidates...")

        from src.analysis.liabilities import LiabilityScanner
        from src.analysis.humanness import HumannessScorer
        from src.analysis.developability import DevelopabilityScorer
        from src.structure.interface_analysis import InterfaceAnalyzer

        liability_scanner = LiabilityScanner()
        humanness_scorer = HumannessScorer()
        developability_scorer = DevelopabilityScorer()
        interface_analyzer = InterfaceAnalyzer()

        scored_candidates = []

        for candidate in candidates:
            score = CandidateScore(
                candidate_id=candidate.get("design_id", candidate.get("name", "unknown")),
                sequence=candidate.get("sequence", candidate.get("vh", "")),
                sequence_vl=candidate.get("sequence_vl", candidate.get("vl")),
                binder_type=candidate.get("binder_type", "vhh"),
                source=candidate.get("source", "unknown"),
            )

            # Structure prediction metrics
            if candidate.get("structure_prediction"):
                sp = candidate["structure_prediction"]
                score.pdockq = sp.get("pdockq")
                score.interface_area = sp.get("interface_area")
                score.num_contacts = sp.get("num_contacts")

                # Epitope annotation
                if sp.get("interface_residues_target"):
                    epitope_class, overlap = interface_analyzer.annotate_epitope_class(
                        sp["interface_residues_target"],
                        self.config.epitope.overlap_threshold,
                    )
                    score.epitope_class = epitope_class
                    score.okt3_overlap = overlap

            # Liability analysis
            try:
                liabilities = liability_scanner.scan_sequence(score.sequence)
                score.deamidation_sites = liabilities.get("deamidation", [])
                score.isomerization_sites = liabilities.get("isomerization", [])
                score.glycosylation_sites = liabilities.get("glycosylation", [])
                score.oxidation_sites = liabilities.get("oxidation", [])
                score.unpaired_cys = liabilities.get("unpaired_cysteine", 0)
            except Exception as e:
                print(f"  Warning: Liability analysis failed: {e}")

            # Humanness scoring
            try:
                humanness = humanness_scorer.score_oasis(score.sequence, score.sequence_vl)
                score.oasis_score_vh = humanness.get("vh_oasis")
                score.oasis_score_vl = humanness.get("vl_oasis")
                score.oasis_score_mean = humanness.get("mean_oasis")
            except Exception as e:
                print(f"  Warning: Humanness scoring failed: {e}")

            # Developability scoring
            try:
                dev = developability_scorer.score(score.sequence)
                score.cdr_h3_length = dev.get("cdr_h3_length")
                score.net_charge = dev.get("net_charge")
                score.isoelectric_point = dev.get("isoelectric_point")
                score.hydrophobic_patches = dev.get("hydrophobic_patches", 0)
            except Exception as e:
                print(f"  Warning: Developability scoring failed: {e}")

            scored_candidates.append(score)

        return scored_candidates

    def run_filtering(
        self,
        candidates: list[CandidateScore],
    ) -> tuple[list[CandidateScore], dict]:
        """Run filtering cascade.

        Args:
            candidates: List of CandidateScore objects.

        Returns:
            Tuple of (filtered_candidates, filter_stats).
        """
        print(f"Running filter cascade on {len(candidates)} candidates...")

        filtered, stats = run_filter_cascade(
            candidates=candidates,
            config=self.config,
            min_candidates=self.config.filtering.min_candidates,
        )

        print(f"Filter results:")
        print(f"  Input: {stats['total_input']}")
        print(f"  Passing (first pass): {stats['passing_first_pass']}")
        print(f"  Final passing: {stats['final_passing']}")
        if stats.get("used_fallback"):
            print(f"  Used fallback: {len(stats['relaxations_applied'])} relaxations")

        self.result.filter_stats = stats
        self.result.filtered_candidates = filtered

        return filtered, stats

    def select_final_candidates(
        self,
        candidates: list[CandidateScore],
    ) -> list[CandidateScore]:
        """Select final candidates with diversity.

        Args:
            candidates: Filtered candidates.

        Returns:
            Final selected candidates.
        """
        print(f"Selecting {self.config.output.num_final_candidates} final candidates...")

        # Already sorted by composite score
        # Implement diversity selection (simplified version)

        final = []
        seen_epitopes = {"OKT3-like": 0, "novel_epitope": 0}
        seen_types = {"vhh": 0, "scfv": 0}
        seen_sources = {"denovo": 0, "optimized": 0}

        for candidate in candidates:
            if len(final) >= self.config.output.num_final_candidates:
                break

            # Ensure diversity (simplified - could use clustering)
            # For now, just take top candidates
            final.append(candidate)

            seen_epitopes[candidate.epitope_class] = seen_epitopes.get(candidate.epitope_class, 0) + 1
            seen_types[candidate.binder_type] = seen_types.get(candidate.binder_type, 0) + 1
            seen_sources[candidate.source] = seen_sources.get(candidate.source, 0) + 1

        print(f"Selected {len(final)} candidates:")
        print(f"  By epitope: {seen_epitopes}")
        print(f"  By type: {seen_types}")
        print(f"  By source: {seen_sources}")

        self.result.final_candidates = final
        return final

    def run_formatting(
        self,
        candidates: list[CandidateScore],
    ) -> dict[str, list[dict]]:
        """Convert candidates to bispecific formats.

        Args:
            candidates: Final candidates.

        Returns:
            Dictionary mapping format names to formatted sequences.
        """
        print(f"Converting {len(candidates)} candidates to bispecific formats...")

        from src.formatting import format_all

        formatted = {}

        for candidate in candidates:
            # Need placeholder target sequences
            # In real use, these would come from config or user input
            target_vh = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
            target_vl = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

            try:
                constructs = format_all(
                    target_vh=target_vh,
                    target_vl=target_vl,
                    cd3_binder=candidate.sequence,
                    cd3_binder_vl=candidate.sequence_vl,
                    name_prefix=candidate.candidate_id,
                    formats=self.config.formatting.formats,
                )

                for fmt_name, construct in constructs.items():
                    if fmt_name not in formatted:
                        formatted[fmt_name] = []
                    formatted[fmt_name].append({
                        "candidate_id": candidate.candidate_id,
                        "construct": construct.to_dict() if hasattr(construct, "to_dict") else str(construct),
                    })

            except Exception as e:
                print(f"  Warning: Formatting failed for {candidate.candidate_id}: {e}")

        return formatted

    def run(
        self,
        run_calibration: bool = True,
        use_modal: bool = True,
    ) -> PipelineResult:
        """Run complete pipeline.

        Args:
            run_calibration: If True, run calibration first.
            use_modal: If True, use Modal for GPU compute.

        Returns:
            PipelineResult with all results.
        """
        self.result.timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.result.provenance = get_provenance()

        print("=" * 60)
        print("CD3 Binder Design Pipeline")
        print("=" * 60)

        # Step 0: Calibration
        if run_calibration:
            self.run_calibration(use_modal=use_modal)

        # Step 1: Design generation
        denovo = self.run_denovo_design(use_modal=use_modal)
        optimized = self.run_optimization()

        all_candidates = denovo + optimized
        print(f"Total candidates: {len(all_candidates)}")

        # Step 2: Structure prediction
        all_candidates = self.run_structure_prediction(all_candidates, use_modal=use_modal)

        # Step 3: Analysis
        scored = self.run_analysis(all_candidates)
        self.result.all_candidates = scored

        # Step 4: Filtering
        filtered, stats = self.run_filtering(scored)

        # Step 5: Final selection
        final = self.select_final_candidates(filtered)

        # Step 6: Formatting
        formatted = self.run_formatting(final)

        # Save results
        output_dir = Path(self.config.output.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_path = output_dir / f"pipeline_results_{self.result.timestamp}.json"
        self.result.save(str(output_path))
        print(f"\nResults saved to: {output_path}")

        print("\n" + "=" * 60)
        print("Pipeline complete!")
        print(f"Final candidates: {len(final)}")
        print("=" * 60)

        return self.result


def run_full_pipeline(
    config_path: Optional[str] = None,
    run_calibration: bool = True,
    use_modal: bool = True,
) -> PipelineResult:
    """Convenience function to run full pipeline.

    Args:
        config_path: Path to config YAML (optional).
        run_calibration: If True, run calibration first.
        use_modal: If True, use Modal for GPU compute.

    Returns:
        PipelineResult with all results.
    """
    if config_path:
        config = PipelineConfig.load(config_path)
    else:
        config = PipelineConfig()

    pipeline = DesignPipeline(config)
    return pipeline.run(run_calibration=run_calibration, use_modal=use_modal)
