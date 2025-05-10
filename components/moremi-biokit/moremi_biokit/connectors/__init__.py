"""
Connectors Sub-package for Moremi Biokit.

This sub-package provides modules for interacting with various external
biological databases and APIs. Each module within `connectors` is responsible
for handling communication, data retrieval, and potentially data submission
for a specific external resource.

Modules included:
- rcsb_pdb: For interacting with the RCSB Protein Data Bank.
- ncbi: For interacting with NCBI services like Entrez and BLAST.
- uniprot: For interacting with the UniProt knowledgebase.
- _utils: Internal helper utilities for connector modules (not for public use).

These connectors are designed to be used by higher-level modules within
`moremi_biokit` (e.g., `data_retrieval`, `proteins.analysis_tools`) to
abstract away the complexities of direct API interactions.

Example Usage (conceptual, depends on functions exposed here):
    from moremi_biokit.connectors import fetch_rcsb_pdb_entry # If re-exported
    # or
    from moremi_biokit.connectors.rcsb import fetch_entry_summary
"""

from .rcsb import (
    download_pdb_from_rcsb,
    fetch_pdb_verification_summary,
    search_rcsb_by_sequence,
    search_rcsb_by_text
   
)
from .uniprot import(
    download_pdb_from_uniprot
)
from .ncbi import(
    download_pdb_from_ncbi
)

__all__ = [
    "download_pdb_from_rcsb",
    "fetch_pdb_verification_summary",
    "search_rcsb_by_sequence",
    "search_rcsb_by_text",
    
    "download_pdb_from_ncbi",
    
    "download_pdb_from_uniprot",
]