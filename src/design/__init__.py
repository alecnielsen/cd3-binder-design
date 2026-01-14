"""Design modules for CD3 binder generation.

This package provides:
- De novo design using BoltzGen
- Optimization of existing binders
- Affinity variant generation
"""

from src.design.boltzgen_runner import (
    BoltzGenConfig,
    BoltzGenDesign,
    BoltzGenRunner,
    run_boltzgen_batch,
)
from src.design.denovo_design import (
    DeNovoDesignConfig,
    DeNovoDesignResult,
    DeNovoDesigner,
    run_denovo_design,
)
from src.design.optimization import (
    AntibodySequences,
    OptimizedVariant,
    SequenceOptimizer,
    optimize_existing_binders,
)
from src.design.affinity_variants import (
    AffinityMutation,
    AffinityVariant,
    AffinityMutationLibrary,
    AffinityVariantGenerator,
    generate_affinity_variants,
)

__all__ = [
    # BoltzGen runner
    "BoltzGenConfig",
    "BoltzGenDesign",
    "BoltzGenRunner",
    "run_boltzgen_batch",
    # De novo design
    "DeNovoDesignConfig",
    "DeNovoDesignResult",
    "DeNovoDesigner",
    "run_denovo_design",
    # Optimization
    "AntibodySequences",
    "OptimizedVariant",
    "SequenceOptimizer",
    "optimize_existing_binders",
    # Affinity variants
    "AffinityMutation",
    "AffinityVariant",
    "AffinityMutationLibrary",
    "AffinityVariantGenerator",
    "generate_affinity_variants",
]
