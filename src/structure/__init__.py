"""Structure prediction and analysis modules.

This package provides:
- ABodyBuilder2 wrapper for antibody structure prediction
- Boltz-2 wrapper for complex structure prediction
- PDB utilities for parsing and manipulation
- Interface analysis for binding assessment
"""

from src.structure.abodybuilder import (
    ABodyBuilderResult,
    ABodyBuilder,
    predict_antibody_structure,
)
from src.structure.boltz_complex import (
    ComplexPredictionResult,
    Boltz2Predictor,
    predict_binder_complex,
    run_calibration,
)
from src.structure.pdb_utils import (
    extract_sequence_from_pdb,
    extract_chain,
    get_chain_ids,
    get_residue_positions,
    calculate_interface_residues,
    count_contacts,
    estimate_interface_area,
    renumber_chain,
    combine_pdbs,
)
from src.structure.interface_analysis import (
    InterfaceMetrics,
    EpitopeComparison,
    InterfaceAnalyzer,
    analyze_complex_interface,
    batch_epitope_annotation,
)

__all__ = [
    # ABodyBuilder
    "ABodyBuilderResult",
    "ABodyBuilder",
    "predict_antibody_structure",
    # Boltz-2
    "ComplexPredictionResult",
    "Boltz2Predictor",
    "predict_binder_complex",
    "run_calibration",
    # PDB utilities
    "extract_sequence_from_pdb",
    "extract_chain",
    "get_chain_ids",
    "get_residue_positions",
    "calculate_interface_residues",
    "count_contacts",
    "estimate_interface_area",
    "renumber_chain",
    "combine_pdbs",
    # Interface analysis
    "InterfaceMetrics",
    "EpitopeComparison",
    "InterfaceAnalyzer",
    "analyze_complex_interface",
    "batch_epitope_annotation",
]
