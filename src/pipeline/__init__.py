"""Pipeline orchestration modules.

This package provides:
- Configuration management
- Multi-stage filtering cascade
- End-to-end pipeline orchestration
- Report generation
"""

from src.pipeline.config import (
    DesignConfig,
    CalibrationConfig,
    FilteringConfig,
    FormattingConfig,
    OutputConfig,
    ReproducibilityConfig,
    EpitopeConfig,
    PipelineConfig,
    get_provenance,
    create_default_config,
)
from src.pipeline.filter_cascade import (
    FilterResult,
    CandidateScore,
    FilterCascade,
    run_filter_cascade,
)
from src.pipeline.design_pipeline import (
    PipelineResult,
    DesignPipeline,
    run_full_pipeline,
)
from src.pipeline.report_generator import (
    ReportConfig,
    ReportGenerator,
    generate_report,
)

__all__ = [
    # Config
    "DesignConfig",
    "CalibrationConfig",
    "FilteringConfig",
    "FormattingConfig",
    "OutputConfig",
    "ReproducibilityConfig",
    "EpitopeConfig",
    "PipelineConfig",
    "get_provenance",
    "create_default_config",
    # Filter cascade
    "FilterResult",
    "CandidateScore",
    "FilterCascade",
    "run_filter_cascade",
    # Pipeline
    "PipelineResult",
    "DesignPipeline",
    "run_full_pipeline",
    # Report
    "ReportConfig",
    "ReportGenerator",
    "generate_report",
]
