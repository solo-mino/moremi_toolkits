"""
RCSB PDB Connector for Moremi Biokit.

This module provides functions to interact with the RCSB Protein Data Bank (PDB)
APIs (both REST and GraphQL where appropriate). It allows for searching,
downloading PDB files, and fetching metadata associated with PDB entries.

Key functionalities include:
- Downloading PDB files in various formats (PDB, mmCIF, etc).
- Fetching curated summary metadata for PDB entries (for verification or display).
- Performing text-based searches of the PDB.
- Performing sequence-based searches (e.g., using RCSB's sequence search capabilities).
- (Potentially) Performing chemical/ligand-based searches.

This module aims to encapsulate the details of RCSB API endpoints, query construction,
and response parsing, providing a simpler interface for other parts of `moremi_biokit`.
"""

import os
from typing import Optional
from ._utils import make_api_request, APIRequestError

def download_pdb_from_rcsb(pdb_id: str, output_dir: str, file_format: str = 'pdb') -> Optional[str]:
    """Download a structure file from RCSB PDB.

    Args:
        pdb_id (str): The 4-character PDB ID (case-insensitive).
        output_dir (str): Directory to save the downloaded file.
        file_format (str, optional): Desired file format ('pdb', 'cif', 'xml', 'fasta', 'pdb.gz', 'cif.gz'). 
                                    Defaults to 'pdb'.

    Returns:
        Optional[str]: The absolute path to the downloaded file on success, None otherwise.
    """
    if not isinstance(pdb_id, str) or len(pdb_id) != 4:
        print(f"Error: Invalid PDB ID '{pdb_id}'. Must be a 4-character string.")
        return None
        
    pdb_id_lower = pdb_id.lower()
    # Construct the download URL
    # Example: https://files.rcsb.org/download/1xyz.pdb
    download_url = f"https://files.rcsb.org/download/{pdb_id_lower}.{file_format}"
    
    # Ensure output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error creating output directory '{output_dir}': {e}")
        return None
        
    output_filename = f"{pdb_id_lower}.{file_format}"
    output_path = os.path.abspath(os.path.join(output_dir, output_filename))

    print(f"Attempting to download {pdb_id} ({file_format}) from RCSB to {output_path}...")
    
    try:
        # Use the utility function for the request
        response = make_api_request(
            url=download_url,
            method="GET",
            stream=True,
            timeout=120, # Can adjust timeout here or rely on default from _utils
            retries=2    # Can adjust retries here or rely on default from _utils
        )

        # Write the content to the file
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        print(f"Successfully downloaded {output_filename}.")
        return output_path

    except APIRequestError as api_err:
        # The make_api_request function will print warnings during retries.
        # Here, we print a final error message.
        if api_err.status_code == 404:
            print(f"Error: PDB ID '{pdb_id}' in format '{file_format}' not found on RCSB PDB (URL: {api_err.url}). Final error: {api_err}")
        else:
            print(f"Error downloading from RCSB after retries: {api_err}")
        return None
    except IOError as io_err:
        print(f"File error saving downloaded PDB to {output_path}: {io_err}")
        # Clean up partially downloaded file if it exists
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except OSError:
                pass # Ignore error during cleanup
        return None
    except Exception as e:
        print(f"An unexpected error occurred during RCSB download: {e}")
        return None