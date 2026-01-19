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
        scfv_linker = self.config.formatting.scfv_linker

        for name in self.config.calibration.positive_controls:
            try:
                seq = optimizer.load_starting_sequence(name)
                # For paired antibodies (VH+VL), construct scFv for accurate calibration
                if seq.vl:
                    binder_seq = seq.vh + scfv_linker + seq.vl
                    print(f"  Loaded {name} as scFv ({len(binder_seq)} aa)")
                else:
                    binder_seq = seq.vh
                    print(f"  Loaded {name} as VHH ({len(binder_seq)} aa)")
                known_sequences.append(binder_seq)
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
            pdockq_margin=self.config.calibration.pdockq_margin,
            interface_area_margin=self.config.calibration.interface_area_margin,
            contacts_margin=self.config.calibration.contacts_margin,
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

        # Default to first target structure
        default_target = self.config.design.target_structures[0] if self.config.design.target_structures else None
        scfv_linker = self.config.formatting.scfv_linker

        for i, candidate in enumerate(candidates):
            try:
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
                    continue

                # Use target structure from candidate if available, otherwise default
                target_pdb = candidate.get("target_structure", default_target)
                if target_pdb is None:
                    print(f"  Warning: No target structure for candidate {i}")
                    candidate["structure_prediction"] = None
                    continue

                result = predictor.predict_complex(
                    binder_sequence=binder_sequence,
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
                    "target_structure": target_pdb,
                    "binder_sequence_used": binder_sequence,
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
        from src.analysis.humanness import score_humanness_pair
        from src.analysis.developability import DevelopabilityAssessor
        from src.structure.interface_analysis import InterfaceAnalyzer

        liability_scanner = LiabilityScanner()
        developability_assessor = DevelopabilityAssessor()
        interface_analyzer = InterfaceAnalyzer(
            okt3_epitope_residues=self.config.epitope.okt3_epitope_residues
        )

        scored_candidates = []

        for candidate in candidates:
            # Extract sequences - handle both single sequence and vh/vl pairs
            vh_seq = candidate.get("sequence") or candidate.get("vh", "")
            vl_seq = candidate.get("sequence_vl") or candidate.get("vl")

            # Determine binder type - check if single sequence is actually an scFv
            binder_type = candidate.get("binder_type")
            if binder_type is None:
                if vl_seq is not None:
                    binder_type = "scfv"
                else:
                    # Check if the single sequence is a concatenated scFv
                    from src.utils.constants import parse_scfv, is_likely_scfv
                    if is_likely_scfv(vh_seq):
                        parsed = parse_scfv(vh_seq)
                        if parsed:
                            vh_seq, vl_seq = parsed
                            binder_type = "scfv"
                        else:
                            binder_type = "vhh"
                    else:
                        binder_type = "vhh"

            score = CandidateScore(
                candidate_id=candidate.get("design_id", candidate.get("name", "unknown")),
                sequence=vh_seq,
                sequence_vl=vl_seq,
                binder_type=binder_type,
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

            # Liability analysis - scan with CDR detection for accurate filtering
            try:
                # Scan VH with CDR detection
                vh_report = liability_scanner.scan_with_cdr_detection(vh_seq, chain_type="H")

                # Scan VL if present
                if vl_seq:
                    vl_report = liability_scanner.scan_with_cdr_detection(vl_seq, chain_type="L")
                    # Combine reports - offset VL positions by VH length
                    vh_len = len(vh_seq)
                    all_deamidation = vh_report.deamidation_sites + [
                        type(s)(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                        for s in vl_report.deamidation_sites
                    ]
                    all_isomerization = vh_report.isomerization_sites + [
                        type(s)(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                        for s in vl_report.isomerization_sites
                    ]
                    all_glycosylation = vh_report.glycosylation_sites + [
                        type(s)(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                        for s in vl_report.glycosylation_sites
                    ]
                    all_oxidation = vh_report.oxidation_sites + [
                        type(s)(s.motif, s.position + vh_len, s.liability_type, s.in_cdr, s.cdr_name, s.severity)
                        for s in vl_report.oxidation_sites
                    ]
                else:
                    all_deamidation = vh_report.deamidation_sites
                    all_isomerization = vh_report.isomerization_sites
                    all_glycosylation = vh_report.glycosylation_sites
                    all_oxidation = vh_report.oxidation_sites

                # Extract positions as integers for JSON serialization
                score.deamidation_sites = [s.position for s in all_deamidation]
                score.isomerization_sites = [s.position for s in all_isomerization]
                score.glycosylation_sites = [s.position for s in all_glycosylation]
                score.oxidation_sites = [s.position for s in all_oxidation]
                score.unpaired_cys = vh_report.unpaired_cysteines + (vl_report.unpaired_cysteines if vl_seq else 0)

                # Count CDR-specific liabilities for filtering
                score.cdr_deamidation_count = sum(1 for s in all_deamidation if s.in_cdr)
                score.cdr_isomerization_count = sum(1 for s in all_isomerization if s.in_cdr)
                score.cdr_glycosylation_count = sum(1 for s in all_glycosylation if s.in_cdr)
                score.cdr_oxidation_count = sum(1 for s in all_oxidation if s.in_cdr)
            except Exception as e:
                print(f"  Warning: Liability analysis failed: {e}")

            # Humanness scoring
            try:
                humanness_report = score_humanness_pair(vh_seq, vl_seq)
                score.oasis_score_vh = humanness_report.vh_report.oasis_score
                score.oasis_score_vl = humanness_report.vl_report.oasis_score if humanness_report.vl_report else None
                score.oasis_score_mean = humanness_report.mean_score
            except Exception as e:
                print(f"  Warning: Humanness scoring failed: {e}")

            # Developability scoring
            try:
                dev_report = developability_assessor.assess(vh_seq, vl_seq, include_humanness=False)
                score.cdr_h3_length = dev_report.cdr_h3_length
                score.net_charge = dev_report.physicochemical.net_charge
                score.isoelectric_point = dev_report.physicochemical.isoelectric_point
                score.hydrophobic_patches = dev_report.aggregation.hydrophobic_patches
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

        from src.formatting import format_all, load_target_sequences

        # Load target sequences from placeholder_targets.yaml based on config
        tumor_target_name = self.config.formatting.tumor_target
        try:
            target_vh, target_vl, target_display = load_target_sequences(tumor_target_name)
            print(f"  Target arm: {tumor_target_name} ({target_display})")
        except (FileNotFoundError, ValueError) as e:
            print(f"  ERROR: Failed to load target sequences: {e}")
            return {}

        formatted = {}

        for candidate in candidates:
            try:
                constructs = format_all(
                    target_vh=target_vh,
                    target_vl=target_vl,
                    cd3_binder=candidate.sequence,
                    cd3_binder_vl=candidate.sequence_vl,
                    name_prefix=candidate.candidate_id,
                    target_name=target_display,
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
