"""
Tools for Small Molecule (SMILES) Processing and Analysis.

This subpackage provides a comprehensive suite of tools for handling SMILES strings,
including batch processing, validation against various molecular property criteria,
calculation of diverse physicochemical and ADMET properties, ranking based on
customizable scoring profiles, and generation of detailed reports.

Key Components:
- `BatchMoleculeProcessor`: Orchestrates the processing of multiple SMILES strings from input files.
- `SmallMoleculeValidator`: Validates individual SMILES strings and gathers a wide array of metrics.
- `SmallMoleculeRankerV4`: Ranks molecules based on calculated metrics and a scoring configuration.
- `MoleculeMetrics`: A data class holding all metrics for a single molecule.
- `ProcessingResult`: A data class representing the outcome of processing a single SMILES string.
- `MetricCategory`: Enum for different categories of molecular metrics.
- `ScoringConfig`: Data class for configuring the scoring behavior of the ranker.

Submodules:
- `property_calculators`: Contains individual modules for calculating specific sets of properties 
  (e.g., physicochemical, lipophilicity, ADMET).
- `reports`: Handles the generation of various reports (e.g., PDF, CSV).
- `utils`: Provides utility functions, such as SMILES conversion tools.

Example Usage:
    To process a batch of SMILES from a file:
    >>> from moremi_biokit.smiles import BatchMoleculeProcessor
    >>> processor = BatchMoleculeProcessor(input_file="molecules.smi", output_dir="analysis_results")
    >>> processor.process_batch()

    To validate a single SMILES string:
    >>> from moremi_biokit.smiles import SmallMoleculeValidator, MoleculeMetrics
    >>> validator = SmallMoleculeValidator()
    >>> result = validator.process_molecule("CCO")
    >>> if result.success and isinstance(result.metrics, MoleculeMetrics):
    ...     print(f"Ethanol MW: {result.metrics.molecular_weight}")
"""

# Core components
from .batch_molecule_processor import BatchMoleculeProcessor
from .small_molecule_validator_v3 import (
    SmallMoleculeValidator,
    MoleculeMetrics,
    ProcessingResult,
    MetricCategory,
)
from .small_molecule_ranker_v4 import (
    SmallMoleculeRankerV4,
    ScoringConfig
)

# Expose submodules for easier access if needed by users
from . import property_calculators
from . import reports
from . import utils

__all__ = [
    "BatchMoleculeProcessor",
    "SmallMoleculeValidator",
    "SmallMoleculeRankerV4",
    "MoleculeMetrics",
    "ProcessingResult",
    "MetricCategory",
    "ScoringConfig",
    "property_calculators",
    "reports",
    "utils",
] 