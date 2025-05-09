"""
UniProt Connector for Moremi Biokit.

This module provides functions for accessing data from the UniProt knowledgebase,
a comprehensive resource for protein sequence and functional information.
It primarily interacts with the UniProt REST API.

Key functionalities include:
- Retrieving protein entries by UniProt accession numbers.
- Searching UniProt for proteins based on various criteria (keywords, gene names, organism, etc.).
- Fetching specific information fields from UniProt entries (e.g., sequence, function,
  PDB cross-references, pathway information).
- Batch retrieval of data for multiple accession numbers.

This module abstracts the UniProt API calls and data parsing, making it easier
to integrate UniProt data into `moremi_biokit` workflows.
"""

from typing import Optional

def download_pdb_from_uniprot(uniprot_id: str, output_dir: str, source_db: str = 'pdb') -> Optional[str]:
    """(Placeholder) Download associated structure files for a UniProt ID.

    Args:
        uniprot_id (str): UniProt accession number.
        output_dir (str): Directory to save the downloaded file(s).
        source_db (str, optional): Database to get structures from ('pdb', 'alphafold'). 
                                Defaults to 'pdb'. AlphaFold results might be large.

    Returns:
        Optional[str]: Path to the downloaded file (or perhaps a list? TBD).
                      Returns None on failure.
    """
    raise NotImplementedError("External structure download via UniProt not yet implemented.")