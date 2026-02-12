"""Pipeline configuration management."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import yaml
import json
import hashlib
import datetime


@dataclass
class DesignConfig:
    """Configuration for design generation."""

    # De novo design
    num_vhh_designs: int = 200
    num_fab_designs: int = 200  # Fab CDR redesign (replaces broken scFv)
    target_structures: list[str] = field(default_factory=list)

    # Fab scaffold configuration
    fab_scaffolds: list[str] = field(default_factory=lambda: ["adalimumab", "belimumab", "dupilumab"])
    fab_scaffold_dir: str = "data/fab_scaffolds"

    # Optimization
    starting_sequences: list[str] = field(default_factory=lambda: ["teplizumab", "sp34", "ucht1"])
    affinity_variants: list[str] = field(default_factory=lambda: ["wild_type", "10x_weaker", "100x_weaker"])


@dataclass
class CalibrationConfig:
    """Configuration for threshold calibration."""

    positive_controls: list[str] = field(default_factory=lambda: ["teplizumab", "sp34", "ucht1"])

    # Margins to subtract from known binder minimums
    pdockq_margin: float = 0.05
    interface_area_margin: float = 100.0
    contacts_margin: int = 2

    # Whether to run ProteinMPNN, AntiFold, Protenix on controls for baselines
    run_validation_baselines: bool = True


@dataclass
class FilteringConfig:
    """Configuration for filtering cascade."""

    # Binding quality (defaults - use calibrated if available)
    min_pdockq: float = 0.5
    min_interface_area: float = 800.0
    min_contacts: int = 10
    use_calibrated: bool = True

    # Humanness
    min_oasis_score: float = 0.8
    generate_back_mutations: bool = True

    # Liabilities
    allow_deamidation_cdr: bool = False
    allow_isomerization_cdr: bool = False
    allow_glycosylation_cdr: bool = False
    max_oxidation_sites: int = 2  # Soft filter

    # Developability
    cdr_h3_length_range: tuple[int, int] = (8, 20)
    net_charge_range: tuple[int, int] = (-2, 4)
    pi_range: tuple[float, float] = (6.0, 9.0)
    max_hydrophobic_patches: int = 2

    # Fallback
    min_candidates: int = 10
    relax_soft_filters_first: bool = True
    max_threshold_relaxation: float = 0.1  # 10%


@dataclass
class FormattingConfig:
    """Configuration for bispecific formatting."""

    tumor_target: str = "trastuzumab"
    formats: list[str] = field(default_factory=lambda: [
        "crossmab", "fab_scfv", "fab_vhh", "igg_scfv", "igg_vhh"
    ])
    scfv_linker: str = "GGGGSGGGGSGGGGS"
    fc_fusion_linker: str = "GGGGSGGGGS"


@dataclass
class RankingConfig:
    """Configuration for candidate ranking."""

    # Ranking method: "boltzgen" (default, uses BoltzGen's internal ranking),
    # "worst_metric_rank" (custom re-ranking), or "composite" (legacy)
    method: str = "boltzgen"

    # Fallback ranking method when primary method's data is unavailable
    secondary_method: str = "worst_metric_rank"

    # Metric weights for worst_metric_rank (only used if method="worst_metric_rank")
    metric_weights: dict[str, int] = field(default_factory=lambda: {
        "iptm": 1,
        "ptm": 1,
        "interface_area": 1,
        "humanness": 1,
        "num_contacts": 2,
        "plddt": 2,
        "proteinmpnn_ll": 1,
        "antifold_ll": 2,
        "protenix_iptm": 1,
    })

    # Diversity selection
    use_diversity_selection: bool = True
    diversity_alpha: float = 0.001  # Greedy maximin alpha


@dataclass
class OutputConfig:
    """Configuration for pipeline output."""

    num_final_candidates: int = 10
    include_structures: bool = True
    generate_report: bool = True
    include_provenance: bool = True
    output_dir: str = "data/outputs"
    export_cif: bool = True


@dataclass
class ReproducibilityConfig:
    """Configuration for reproducibility."""

    boltzgen_seed: int = 42
    sampling_seed: int = 12345
    clustering_seed: int = 0


@dataclass
class EpitopeConfig:
    """Configuration for epitope annotation."""

    # OKT3 epitope residues on CD3ε
    # Default is None, which triggers dynamic extraction from 1SY6 crystal structure.
    # Users can override with explicit residue list if needed.
    okt3_epitope_residues: Optional[list[int]] = None
    overlap_threshold: float = 0.5


@dataclass
class ValidationConfig:
    """Configuration for candidate validation (step 05b)."""

    enabled: bool = True
    run_protenix: bool = True
    run_proteinmpnn: bool = True
    run_antifold: bool = True
    protenix_model: str = "protenix_base_default_v1.0.0"
    protenix_use_msa: bool = False
    protenix_seeds: list[int] = field(default_factory=lambda: [101])
    iptm_disagreement_threshold: float = 0.1


@dataclass
class HumanizationConfig:
    """Configuration for post-hoc humanization (step 04b)."""

    enabled: bool = False
    min_humanness_for_humanization: float = 0.70
    max_humanness_for_humanization: float = 0.80
    num_variants_per_candidate: int = 5
    sample_method: str = "FR"  # Framework-only regeneration (preserves CDRs)
    repredict_structures: bool = True
    rescore_validation: bool = True


@dataclass
class PipelineConfig:
    """Complete pipeline configuration."""

    design: DesignConfig = field(default_factory=DesignConfig)
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    filtering: FilteringConfig = field(default_factory=FilteringConfig)
    formatting: FormattingConfig = field(default_factory=FormattingConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    reproducibility: ReproducibilityConfig = field(default_factory=ReproducibilityConfig)
    epitope: EpitopeConfig = field(default_factory=EpitopeConfig)
    ranking: RankingConfig = field(default_factory=RankingConfig)
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    humanization: HumanizationConfig = field(default_factory=HumanizationConfig)

    # Calibrated thresholds (set after calibration)
    calibrated_min_pdockq: Optional[float] = None
    calibrated_min_interface_area: Optional[float] = None
    calibrated_min_contacts: Optional[int] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "design": {
                "num_vhh_designs": self.design.num_vhh_designs,
                "num_fab_designs": self.design.num_fab_designs,
                "target_structures": self.design.target_structures,
                "fab_scaffolds": self.design.fab_scaffolds,
                "fab_scaffold_dir": self.design.fab_scaffold_dir,
                "starting_sequences": self.design.starting_sequences,
                "affinity_variants": self.design.affinity_variants,
            },
            "calibration": {
                "positive_controls": self.calibration.positive_controls,
                "pdockq_margin": self.calibration.pdockq_margin,
                "interface_area_margin": self.calibration.interface_area_margin,
                "contacts_margin": self.calibration.contacts_margin,
                "run_validation_baselines": self.calibration.run_validation_baselines,
            },
            "filtering": {
                "min_pdockq": self.filtering.min_pdockq,
                "min_interface_area": self.filtering.min_interface_area,
                "min_contacts": self.filtering.min_contacts,
                "use_calibrated": self.filtering.use_calibrated,
                "min_oasis_score": self.filtering.min_oasis_score,
                "generate_back_mutations": self.filtering.generate_back_mutations,
                "allow_deamidation_cdr": self.filtering.allow_deamidation_cdr,
                "allow_isomerization_cdr": self.filtering.allow_isomerization_cdr,
                "allow_glycosylation_cdr": self.filtering.allow_glycosylation_cdr,
                "max_oxidation_sites": self.filtering.max_oxidation_sites,
                "cdr_h3_length_range": self.filtering.cdr_h3_length_range,
                "net_charge_range": self.filtering.net_charge_range,
                "pi_range": self.filtering.pi_range,
                "max_hydrophobic_patches": self.filtering.max_hydrophobic_patches,
                "min_candidates": self.filtering.min_candidates,
                "relax_soft_filters_first": self.filtering.relax_soft_filters_first,
                "max_threshold_relaxation": self.filtering.max_threshold_relaxation,
            },
            "formatting": {
                "tumor_target": self.formatting.tumor_target,
                "formats": self.formatting.formats,
                "scfv_linker": self.formatting.scfv_linker,
                "fc_fusion_linker": self.formatting.fc_fusion_linker,
            },
            "output": {
                "num_final_candidates": self.output.num_final_candidates,
                "include_structures": self.output.include_structures,
                "generate_report": self.output.generate_report,
                "include_provenance": self.output.include_provenance,
                "output_dir": self.output.output_dir,
                "export_cif": self.output.export_cif,
            },
            "reproducibility": {
                "boltzgen_seed": self.reproducibility.boltzgen_seed,
                "sampling_seed": self.reproducibility.sampling_seed,
                "clustering_seed": self.reproducibility.clustering_seed,
            },
            "epitope": {
                "okt3_epitope_residues": self.epitope.okt3_epitope_residues,
                "overlap_threshold": self.epitope.overlap_threshold,
            },
            "ranking": {
                "method": self.ranking.method,
                "secondary_method": self.ranking.secondary_method,
                "metric_weights": self.ranking.metric_weights,
                "use_diversity_selection": self.ranking.use_diversity_selection,
                "diversity_alpha": self.ranking.diversity_alpha,
            },
            "validation": {
                "enabled": self.validation.enabled,
                "run_protenix": self.validation.run_protenix,
                "run_proteinmpnn": self.validation.run_proteinmpnn,
                "run_antifold": self.validation.run_antifold,
                "protenix_model": self.validation.protenix_model,
                "protenix_use_msa": self.validation.protenix_use_msa,
                "protenix_seeds": self.validation.protenix_seeds,
                "iptm_disagreement_threshold": self.validation.iptm_disagreement_threshold,
            },
            "humanization": {
                "enabled": self.humanization.enabled,
                "min_humanness_for_humanization": self.humanization.min_humanness_for_humanization,
                "max_humanness_for_humanization": self.humanization.max_humanness_for_humanization,
                "num_variants_per_candidate": self.humanization.num_variants_per_candidate,
                "sample_method": self.humanization.sample_method,
                "repredict_structures": self.humanization.repredict_structures,
                "rescore_validation": self.humanization.rescore_validation,
            },
            "calibrated_thresholds": {
                "min_pdockq": self.calibrated_min_pdockq,
                "min_interface_area": self.calibrated_min_interface_area,
                "min_contacts": self.calibrated_min_contacts,
            },
        }

    def config_hash(self) -> str:
        """Generate hash of configuration."""
        config_str = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.sha256(config_str.encode()).hexdigest()[:12]

    def get_effective_thresholds(self) -> dict:
        """Get effective filter thresholds (calibrated or default)."""
        if self.filtering.use_calibrated and self.calibrated_min_pdockq is not None:
            return {
                "min_pdockq": self.calibrated_min_pdockq,
                "min_interface_area": self.calibrated_min_interface_area,
                "min_contacts": self.calibrated_min_contacts,
            }
        else:
            return {
                "min_pdockq": self.filtering.min_pdockq,
                "min_interface_area": self.filtering.min_interface_area,
                "min_contacts": self.filtering.min_contacts,
            }

    def save(self, output_path: str) -> str:
        """Save configuration to YAML file."""
        with open(output_path, "w") as f:
            # Use safe_dump to avoid Python-specific types like !!python/tuple
            yaml.safe_dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)
        return output_path

    @classmethod
    def load(cls, config_path: str) -> "PipelineConfig":
        """Load configuration from YAML file.

        Supports both flat and nested schemas for backwards compatibility.
        Nested schema (from README):
            design:
              denovo:
                num_vhh_designs: 200
            filtering:
              binding:
                min_pdockq: 0.5

        Flat schema:
            design:
              num_vhh_designs: 200
            filtering:
              min_pdockq: 0.5
        """
        with open(config_path, "r") as f:
            data = yaml.safe_load(f)

        config = cls()

        # Design - support nested (design.denovo) and flat (design) schemas
        if "design" in data:
            d = data["design"]
            # Check for nested 'denovo' key
            denovo = d.get("denovo", {})
            optimization = d.get("optimization", {})

            config.design.num_vhh_designs = denovo.get("num_vhh_designs", d.get("num_vhh_designs", 200))
            # Support both old (num_scfv_designs) and new (num_fab_designs) keys
            config.design.num_fab_designs = denovo.get("num_fab_designs", d.get("num_fab_designs",
                denovo.get("num_scfv_designs", d.get("num_scfv_designs", 200))))
            config.design.target_structures = denovo.get("target_structures", d.get("target_structures", []))
            config.design.fab_scaffolds = d.get("fab_scaffolds", ["adalimumab", "belimumab", "dupilumab"])
            config.design.fab_scaffold_dir = d.get("fab_scaffold_dir", "data/fab_scaffolds")
            config.design.starting_sequences = optimization.get("starting_sequences", d.get("starting_sequences", ["teplizumab", "sp34", "ucht1"]))
            config.design.affinity_variants = optimization.get("affinity_variants", d.get("affinity_variants", ["wild_type", "10x_weaker", "100x_weaker"]))

        # Calibration - support nested calibration_margin and flat
        if "calibration" in data:
            c = data["calibration"]
            margins = c.get("calibration_margin", {})

            config.calibration.positive_controls = c.get("positive_controls", ["teplizumab", "sp34", "ucht1"])
            config.calibration.pdockq_margin = margins.get("pdockq", c.get("pdockq_margin", 0.05))
            config.calibration.interface_area_margin = margins.get("interface_area", c.get("interface_area_margin", 100.0))
            config.calibration.contacts_margin = margins.get("contacts", c.get("contacts_margin", 2))
            config.calibration.run_validation_baselines = c.get("run_validation_baselines", True)

        # Filtering - support nested (filtering.binding, filtering.humanness, etc.) and flat
        if "filtering" in data:
            f = data["filtering"]
            binding = f.get("binding", {})
            humanness = f.get("humanness", {})
            liabilities = f.get("liabilities", {})
            developability = f.get("developability", {})
            fallback = f.get("fallback", {})

            # Binding thresholds
            config.filtering.min_pdockq = binding.get("min_pdockq", f.get("min_pdockq", 0.5))
            config.filtering.min_interface_area = binding.get("min_interface_area", f.get("min_interface_area", 800.0))
            config.filtering.min_contacts = binding.get("min_contacts", f.get("min_contacts", 10))
            config.filtering.use_calibrated = binding.get("use_calibrated", f.get("use_calibrated", True))

            # Humanness
            config.filtering.min_oasis_score = humanness.get("min_oasis_score", f.get("min_oasis_score", 0.8))
            config.filtering.generate_back_mutations = humanness.get("generate_back_mutations", f.get("generate_back_mutations", True))

            # Liabilities
            config.filtering.allow_deamidation_cdr = liabilities.get("allow_deamidation_cdr", f.get("allow_deamidation_cdr", False))
            config.filtering.allow_isomerization_cdr = liabilities.get("allow_isomerization_cdr", f.get("allow_isomerization_cdr", False))
            config.filtering.allow_glycosylation_cdr = liabilities.get("allow_glycosylation_cdr", f.get("allow_glycosylation_cdr", False))
            config.filtering.max_oxidation_sites = liabilities.get("max_oxidation_sites", f.get("max_oxidation_sites", 2))

            # Developability
            cdr_h3_range = developability.get("cdr_h3_length_range", f.get("cdr_h3_length_range", [8, 20]))
            if isinstance(cdr_h3_range, list):
                config.filtering.cdr_h3_length_range = tuple(cdr_h3_range)
            charge_range = developability.get("net_charge_range", f.get("net_charge_range", [-2, 4]))
            if isinstance(charge_range, list):
                config.filtering.net_charge_range = tuple(charge_range)
            pi_range = developability.get("pi_range", f.get("pi_range", [6.0, 9.0]))
            if isinstance(pi_range, list):
                config.filtering.pi_range = tuple(pi_range)
            config.filtering.max_hydrophobic_patches = developability.get("max_hydrophobic_patches", f.get("max_hydrophobic_patches", 2))

            # Fallback
            config.filtering.min_candidates = fallback.get("min_candidates", f.get("min_candidates", 10))
            config.filtering.relax_soft_filters_first = fallback.get("relax_soft_filters_first", f.get("relax_soft_filters_first", True))
            config.filtering.max_threshold_relaxation = fallback.get("max_threshold_relaxation", f.get("max_threshold_relaxation", 0.1))

        # Formatting
        if "formatting" in data:
            fmt = data["formatting"]
            config.formatting.tumor_target = fmt.get("tumor_target", "trastuzumab")
            config.formatting.formats = fmt.get("formats", ["crossmab", "fab_scfv", "fab_vhh", "igg_scfv", "igg_vhh"])
            linkers = fmt.get("linkers", {})
            config.formatting.scfv_linker = linkers.get("scfv", fmt.get("scfv_linker", "GGGGSGGGGSGGGGS"))
            config.formatting.fc_fusion_linker = linkers.get("fc_fusion", fmt.get("fc_fusion_linker", "GGGGSGGGGS"))

        # Output
        if "output" in data:
            o = data["output"]
            config.output.num_final_candidates = o.get("num_final_candidates", 10)
            config.output.include_structures = o.get("include_structures", True)
            config.output.generate_report = o.get("generate_report", True)
            config.output.include_provenance = o.get("include_provenance", True)
            config.output.output_dir = o.get("output_dir", "data/outputs")
            config.output.export_cif = o.get("export_cif", True)

        # Reproducibility
        if "reproducibility" in data:
            r = data["reproducibility"]
            config.reproducibility.boltzgen_seed = r.get("boltzgen_seed", 42)
            config.reproducibility.sampling_seed = r.get("sampling_seed", 12345)
            config.reproducibility.clustering_seed = r.get("clustering_seed", 0)

        # Epitope - support nested epitope_annotation and flat epitope
        if "epitope_annotation" in data:
            e = data["epitope_annotation"]
            config.epitope.okt3_epitope_residues = e.get("okt3_epitope_residues", config.epitope.okt3_epitope_residues)
            config.epitope.overlap_threshold = e.get("overlap_threshold", 0.5)
        elif "epitope" in data:
            e = data["epitope"]
            config.epitope.okt3_epitope_residues = e.get("okt3_epitope_residues", config.epitope.okt3_epitope_residues)
            config.epitope.overlap_threshold = e.get("overlap_threshold", 0.5)

        # Ranking
        if "ranking" in data:
            r = data["ranking"]
            config.ranking.method = r.get("method", "worst_metric_rank")
            config.ranking.secondary_method = r.get("secondary_method", "worst_metric_rank")
            config.ranking.metric_weights = r.get("metric_weights", config.ranking.metric_weights)
            config.ranking.use_diversity_selection = r.get("use_diversity_selection", True)
            config.ranking.diversity_alpha = r.get("diversity_alpha", 0.001)

        # Validation
        if "validation" in data:
            v = data["validation"]
            config.validation.enabled = v.get("enabled", True)
            config.validation.run_protenix = v.get("run_protenix", True)
            config.validation.run_proteinmpnn = v.get("run_proteinmpnn", True)
            config.validation.run_antifold = v.get("run_antifold", True)
            config.validation.protenix_model = v.get("protenix_model", "protenix_base_default_v1.0.0")
            config.validation.protenix_use_msa = v.get("protenix_use_msa", False)
            config.validation.protenix_seeds = v.get("protenix_seeds", [101])
            config.validation.iptm_disagreement_threshold = v.get("iptm_disagreement_threshold", 0.1)

        # Humanization
        if "humanization" in data:
            h = data["humanization"]
            config.humanization.enabled = h.get("enabled", False)
            config.humanization.min_humanness_for_humanization = h.get("min_humanness_for_humanization", 0.70)
            config.humanization.max_humanness_for_humanization = h.get("max_humanness_for_humanization", 0.80)
            config.humanization.num_variants_per_candidate = h.get("num_variants_per_candidate", 5)
            config.humanization.sample_method = h.get("sample_method", "FR")
            config.humanization.repredict_structures = h.get("repredict_structures", True)
            config.humanization.rescore_validation = h.get("rescore_validation", True)

        # Calibrated thresholds
        if "calibrated_thresholds" in data:
            ct = data["calibrated_thresholds"]
            config.calibrated_min_pdockq = ct.get("min_pdockq")
            config.calibrated_min_interface_area = ct.get("min_interface_area")
            config.calibrated_min_contacts = ct.get("min_contacts")

        return config


def get_provenance() -> dict:
    """Get provenance information for output files."""
    import subprocess

    provenance = {
        "pipeline_version": "1.0.0",
        "run_timestamp": datetime.datetime.now().isoformat(),
    }

    # Try to get git commit
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            provenance["git_commit"] = result.stdout.strip()[:12]
    except Exception:
        pass

    return provenance


def create_default_config(output_path: Optional[str] = None) -> PipelineConfig:
    """Create a default configuration.

    Args:
        output_path: If provided, save config to this path.

    Returns:
        Default PipelineConfig.
    """
    config = PipelineConfig()

    # Set default target structures
    config.design.target_structures = [
        "data/targets/cd3_epsilon_delta_1XIW.pdb",
        "data/targets/cd3_epsilon_gamma_1SY6.pdb",
    ]

    if output_path:
        config.save(output_path)
        print(f"Default config saved to: {output_path}")

    return config
