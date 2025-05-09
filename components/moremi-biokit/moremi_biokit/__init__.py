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
    
    To access the proteins subpackage directly:
    >>> from moremi_biokit import proteins
    >>> processor = proteins.BatchProteinProcessor("input.txt", "output_dir")
    >>> # processor.process_btach()
"""
from .pdb_fetcher import (
    list_internal_pdb_ids,
    get_internal_pdb_path,
    read_internal_pdb,
    fetch_pdb
)

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

from .proteins import (
    BatchAntibodyProcessor,
    ProteinValidator,
    ProteinRanker,
    MetricCategory,
    ProcessingResult,
    ProteinMetrics,
    ScoringConfig,
)

# Also make the subpackages themselves available
from . import smiles
from . import proteins
from . import connectors

__all__ = [
    # From smiles subpackage
    "BatchMoleculeProcessor",
    "SmallMoleculeValidator",
    "SmallMoleculeRankerV4",
    "MoleculeMetrics",
    "ProcessingResult",
    "MetricCategory",
    "ScoringConfig",
    
    "BatchAntibodyProcessor",
    "ProteinValidator",
    "ProteinRanker",
    "MetricCategory",
    "ProcessingResult",
    "ProteinMetrics",
    "ScoringConfig",
    
    # Subpackages
    "smiles",
    "proteins",
    "connectors",
    
    # PDB Fetcher functions
    "list_internal_pdb_ids",
    "get_internal_pdb_path",
    "read_internal_pdb",
    "fetch_pdb",

]

__version__ = "0.1.0"