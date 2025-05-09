"""
Moremi BioKit - Proteins Subpackage
=====================================

This subpackage provides a comprehensive suite of tools for protein
analysis, validation, ranking, and reporting.

Key Modules:
- `analysis_tools`: Core bioinformatics tools for property prediction.
- `protein_validator`: Validates protein sequences and collects metrics.
- `protein_ranker`: Ranks proteins based on weighted scores.
- `batch_protein_processor`: Processes batches of proteins.
- `reports`: Generates PDF and CSV reports.
- `utils`: Utility functions, including PDB fetching.
"""

# Import key components from within the proteins subpackage
from .batch_protein_processor import BatchAntibodyProcessor
from .protein_validator import (
    ProteinValidator,
    ProteinMetrics,
    MetricCategory,
    ProcessingResult
)
from .protein_ranker import (
    ProteinRanker,
    ScoringConfig
)

# Expose all submodules for easier access if needed by users
from . import analysis_tools
from . import utils
from . import reports

__all__ = [
    # From protein_validator
    "ProteinValidator",
    "ProteinMetrics",
    "MetricCategory",
    "ProcessingResult",

    # From protein_ranker
    "ProteinRanker",
    "ScoringConfig",

    # From batch_protein_processor
    "BatchAntibodyProcessor",
    
    "utils",
    "reports",
    "analysis_tools"
]