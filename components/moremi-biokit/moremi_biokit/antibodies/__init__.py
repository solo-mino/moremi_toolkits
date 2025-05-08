"""
Moremi BioKit - Antibodies Subpackage
=====================================

This subpackage provides a comprehensive suite of tools for antibody
analysis, validation, ranking, and reporting.

Key Modules:
- `analysis_tools`: Core bioinformatics tools for property prediction.
- `antibody_validator`: Validates antibody sequences and collects metrics.
- `antibody_ranker`: Ranks antibodies based on weighted scores.
- `batch_antibody_processor`: Processes batches of antibodies.
- `reports`: Generates PDF and CSV reports.
- `utils`: Utility functions, including PDB fetching.
"""

# Import key components from within the antibodies subpackage

from .antibody_validator import (
    AntibodyValidator,
    AntibodyMetrics,
    MetricCategory,
    ProcessingResult
)
from .antibody_ranker import (
    AntibodyRanker,
    ScoringConfig
)
from .batch_antibody_processor import BatchAntibodyProcessor

# Re-export all from analysis_tools
from .analysis_tools import * # noqa: F403

# Re-export key components from utils
from .utils.pdb_fetcher import (
    list_internal_pdb_ids,
    get_internal_pdb_path,
    read_internal_pdb,
    fetch_pdb,
    download_pdb_from_rcsb
    # download_structure_from_ncbi, # Not fully implemented
    # download_structure_from_uniprot # Not fully implemented
)

__all__ = [
    # From antibody_validator
    "AntibodyValidator",
    "AntibodyMetrics",
    "MetricCategory",
    "ProcessingResult",

    # From antibody_ranker
    "AntibodyRanker",
    "ScoringConfig",

    # From batch_antibody_processor
    "BatchAntibodyProcessor",

    # From utils.pdb_fetcher
    "list_internal_pdb_ids",
    "get_internal_pdb_path",
    "read_internal_pdb",
    "fetch_pdb",
    "download_pdb_from_rcsb",
    # "download_structure_from_ncbi",
    # "download_structure_from_uniprot",
]

# Add all exports from analysis_tools to __all__
# Assuming analysis_tools.__all__ is correctly defined
try:
    from .analysis_tools import __all__ as analysis_tools_all
    __all__.extend(analysis_tools_all)
except ImportError:
    # Fallback or log if necessary, though analysis_tools should be present
    pass

# Clean up to avoid duplicates if any module also exports something explicitly named above
__all__ = sorted(list(set(__all__)))