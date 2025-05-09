"""
Moremi BioKit (moremi-biokit)

This is the main package for the Moremi BioKit, a foundational Python toolkit 
designed to provide robust processing, validation, analysis, and reporting 
capabilities for common biological entities.

Currently, it primarily includes the `smiles` subpackage for handling small 
molecules represented by SMILES strings.

Subpackages:
  - smiles: Tools for SMILES processing, validation, ranking, and reporting.
  - proteins: Tools for protein sequence analysis.

Example Usage:
    To use the SMILES validation tools:
    >>> from moremi_biokit import SmallMoleculeValidator
    >>> validator = SmallMoleculeValidator()
    >>> result = validator.process_molecule("CCO") # Ethanol
    >>> if result.success:
    ...     print(f"SMILES: {result.metrics.smiles}, MW: {result.metrics.molecular_weight}")

    To access the smiles subpackage directly:
    >>> from moremi_biokit import smiles
    >>> processor = smiles.BatchMoleculeProcessor("input.smi", "output_dir")
    >>> # processor.process_batch()
"""

# Import and re-export key components from the smiles subpackage
from .smiles import (
    BatchMoleculeProcessor,
    SmallMoleculeValidator,
    SmallMoleculeRankerV4,
    MoleculeMetrics,
    ProcessingResult,
    MetricCategory,
    ScoringConfig
)

# Also make the subpackages themselves available
from . import smiles
# from . import proteins # When proteins module is ready

__all__ = [
    # From smiles subpackage
    "BatchMoleculeProcessor",
    "SmallMoleculeValidator",
    "SmallMoleculeRankerV4",
    "MoleculeMetrics",
    "ProcessingResult",
    "MetricCategory",
    "ScoringConfig",
    
    # Subpackages
    "smiles",
    # "proteins", # When proteins module is ready
]

# Define package version (consider centralizing this, e.g., in pyproject.toml and reading here)
__version__ = "0.1.0"