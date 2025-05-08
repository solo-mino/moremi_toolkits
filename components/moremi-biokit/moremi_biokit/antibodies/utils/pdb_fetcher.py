"""Utility for fetching Protein Data Bank (PDB) files.

Provides functions to:
- List available internal PDB files bundled with the package.
- Retrieve paths to or content of internal PDB files.
- (Placeholder) Download PDB files from external databases (PDB, NCBI, UniProt).
"""

import os
import importlib.resources
from typing import List, Optional, Dict, Any
import requests
import argparse

# Define the package path for internal PDBs
_INTERNAL_PDB_PKG_PATH = "moremi_biokit.antibodies.pdb_targets"

def list_internal_pdb_ids() -> List[str]:
    """List the PDB IDs of the internal PDB files bundled with the package.

    Reads the contents of the `moremi_biokit.antibodies.pdb_targets` directory.

    Returns:
        List[str]: A list of PDB IDs (filenames without the .pdb extension).
                   Returns an empty list if the directory doesn't exist or is empty.
    """
    pdb_ids: List[str] = []
    try:
        # Access the package resource directory
        pkg_path = importlib.resources.files(_INTERNAL_PDB_PKG_PATH)
        if pkg_path.is_dir():
            for item in pkg_path.iterdir():
                # Check if it's a file and ends with .pdb (case-insensitive)
                if item.is_file() and item.name.lower().endswith('.pdb'):
                    pdb_ids.append(item.stem) # item.stem gets filename without extension
    except ModuleNotFoundError:
        print(f"Warning: Internal PDB directory '{_INTERNAL_PDB_PKG_PATH}' not found.")
    except Exception as e:
        print(f"Error listing internal PDB files: {e}")
    return sorted(pdb_ids)

def get_internal_pdb_path(pdb_id: str) -> Optional[str]:
    """Get the filesystem path to an internal PDB file.

    Note: This function provides a temporary path managed by importlib.resources.
    It's often better to read the content directly using `read_internal_pdb` 
    unless the path is specifically needed by an external tool.

    Args:
        pdb_id (str): The PDB ID (case-insensitive, without .pdb extension) to retrieve.

    Returns:
        Optional[str]: The absolute path to the PDB file if found, otherwise None.
                       The path is temporary and valid only within a specific context 
                       if used with older importlib.resources versions or certain backends.
                       For modern usage (Python 3.9+), `files()` API gives more stable paths.
    """
    target_filename = f"{pdb_id}.pdb"
    try:
        pkg_path = importlib.resources.files(_INTERNAL_PDB_PKG_PATH)
        pdb_resource = pkg_path.joinpath(target_filename)

        if pdb_resource.is_file():
            # importlib.resources.path is deprecated, use as_file for older compatibility if needed
            # but files().joinpath() should give a usable path object directly in most cases.
            # For simplicity and modern practice, we return the resolved path string.
            # If compatibility with older path-expecting libraries is needed, 
            # consider using `importlib.resources.as_file()` in a `with` block where the path is needed.
            return str(pdb_resource.resolve()) # Resolve to an absolute path string
        else:
            # Try case-insensitive match if exact match failed
            pdb_id_lower = pdb_id.lower()
            for item in pkg_path.iterdir():
                if item.is_file() and item.stem.lower() == pdb_id_lower and item.name.lower().endswith('.pdb'):
                    return str(item.resolve())
            print(f"Warning: Internal PDB file for ID '{pdb_id}' not found in {_INTERNAL_PDB_PKG_PATH}.")
            return None
    except ModuleNotFoundError:
        print(f"Warning: Internal PDB directory '{_INTERNAL_PDB_PKG_PATH}' not found.")
        return None
    except Exception as e:
        print(f"Error getting internal PDB path for '{pdb_id}': {e}")
        return None

def read_internal_pdb(pdb_id: str) -> Optional[str]:
    """Read the content of an internal PDB file.

    Args:
        pdb_id (str): The PDB ID (case-insensitive, without .pdb extension) to read.

    Returns:
        Optional[str]: The content of the PDB file as a string, or None if not found or unreadable.
    """
    target_filename = f"{pdb_id}.pdb"
    try:
        pkg_path = importlib.resources.files(_INTERNAL_PDB_PKG_PATH)
        pdb_resource = pkg_path.joinpath(target_filename)

        if pdb_resource.is_file():
            return pdb_resource.read_text(encoding='utf-8')
        else:
            # Try case-insensitive match
            pdb_id_lower = pdb_id.lower()
            for item in pkg_path.iterdir():
                if item.is_file() and item.stem.lower() == pdb_id_lower and item.name.lower().endswith('.pdb'):
                    return item.read_text(encoding='utf-8')
            print(f"Warning: Internal PDB file for ID '{pdb_id}' not found in {_INTERNAL_PDB_PKG_PATH}.")
            return None
    except ModuleNotFoundError:
        print(f"Warning: Internal PDB directory '{_INTERNAL_PDB_PKG_PATH}' not found.")
        return None
    except Exception as e:
        print(f"Error reading internal PDB content for '{pdb_id}': {e}")
        return None

# --- Placeholder Functions for External PDB Fetching ---

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
        response = requests.get(download_url, stream=True, timeout=120) # Use stream=True for potentially large files, increase timeout
        response.raise_for_status() # Check for HTTP errors (like 404)

        # Write the content to the file
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        print(f"Successfully downloaded {output_filename}.")
        return output_path

    except requests.exceptions.HTTPError as http_err:
        if http_err.response.status_code == 404:
            print(f"Error: PDB ID '{pdb_id}' in format '{file_format}' not found on RCSB PDB (URL: {download_url}).")
        else:
            print(f"HTTP error downloading from RCSB: {http_err} (URL: {download_url})")
        return None
    except requests.exceptions.RequestException as req_err:
        print(f"Request error downloading from RCSB: {req_err} (URL: {download_url})")
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

def download_structure_from_ncbi(ncbi_id: str, output_dir: str) -> Optional[str]:
    """(Placeholder) Download a structure file from NCBI Structure.

    Args:
        ncbi_id (str): NCBI Structure identifier (e.g., MMDB ID).
        output_dir (str): Directory to save the downloaded file.

    Returns:
        Optional[str]: The path to the downloaded file on success, None otherwise.
    """
    raise NotImplementedError("External structure download from NCBI not yet implemented.")

def download_structure_from_uniprot(uniprot_id: str, output_dir: str, source_db: str = 'pdb') -> Optional[str]:
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

# --- Main Fetching Function ---

def fetch_pdb(
    identifier: str,
    source: str = 'rcsb',
    use_internal: bool = False,
    output_dir: Optional[str] = None,
    file_format: str = 'pdb' # Applicable mainly to RCSB download
) -> Dict[str, Any]:
    """Fetch a PDB file, either from internal package data or external sources.

    Args:
        identifier (str): The identifier (e.g., PDB ID, NCBI ID, UniProt ID) 
                          or internal PDB ID.
        source (str, optional): The source database ('rcsb', 'ncbi', 'uniprot') 
                                if `use_internal` is False. Defaults to 'rcsb'.
        use_internal (bool, optional): If True, attempts to retrieve the PDB file 
                                     from the package's internal data using 
                                     `identifier` as the internal ID. Defaults to False.
        output_dir (Optional[str], optional): The directory to save downloaded files. 
                                           Required if `use_internal` is False. 
                                           Defaults to None.
        file_format (str, optional): File format for RCSB download ('pdb', 'cif', etc.).
                                  Defaults to 'pdb'.

    Returns:
        Dict[str, Any]: A dictionary containing:
            - "status" (str): 'success' or 'error'.
            - "source" (str): Where the file was retrieved from ('internal', 'rcsb', 'ncbi', 'uniprot').
            - "identifier" (str): The identifier used.
            - "pdb_path" (Optional[str]): Path to the retrieved/downloaded PDB file (if applicable).
            - "pdb_content" (Optional[str]): Content of the internal PDB file (if retrieved).
            - "message" (str): A descriptive message about the outcome.
    """
    if use_internal:
        pdb_content = read_internal_pdb(identifier)
        if pdb_content:
            # For internal files, providing content is often more useful than a temp path
            return {
                "status": "success",
                "source": "internal",
                "identifier": identifier,
                "pdb_path": None, # Path not typically returned for internal content reads
                "pdb_content": pdb_content,
                "message": f"Successfully read internal PDB: {identifier}"
            }
        else:
             # Attempt to get path as fallback if content read failed but file might exist
             pdb_path = get_internal_pdb_path(identifier)
             if pdb_path:
                 return {
                    "status": "success",
                    "source": "internal",
                    "identifier": identifier,
                    "pdb_path": pdb_path,
                    "pdb_content": None,
                    "message": f"Successfully found internal PDB path: {pdb_path}. Content could not be read or was requested as path."
                 }
             else:
                return {
                    "status": "error",
                    "source": "internal",
                    "identifier": identifier,
                    "pdb_path": None,
                    "pdb_content": None,
                    "message": f"Internal PDB ID '{identifier}' not found."
                }
    else:
        # External download (currently placeholders)
        if not output_dir:
            return {"status": "error", "source": source, "identifier": identifier, "message": "output_dir is required for external downloads."}
        
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            if source.lower() == 'rcsb':
                file_path = download_pdb_from_rcsb(identifier, output_dir, file_format)
            elif source.lower() == 'ncbi':
                file_path = download_structure_from_ncbi(identifier, output_dir)
            elif source.lower() == 'uniprot':
                file_path = download_structure_from_uniprot(identifier, output_dir)
            else:
                return {"status": "error", "source": source, "identifier": identifier, "message": f"Unsupported external source: {source}. Use 'rcsb', 'ncbi', or 'uniprot'."}
            
            # This part will only be reached if the placeholder functions are implemented
            if file_path:
                 return {
                    "status": "success",
                    "source": source,
                    "identifier": identifier,
                    "pdb_path": file_path,
                    "pdb_content": None, # Content not typically read after download
                    "message": f"Successfully downloaded from {source} to {file_path}"
                 }
            else:
                # Should be handled by NotImplementedError now, but good practice
                return {"status": "error", "source": source, "identifier": identifier, "message": f"Download failed from {source} for {identifier}."}

        except NotImplementedError as nie:
            return {"status": "error", "source": source, "identifier": identifier, "message": f"Download function for {source} is not implemented or failed."}
        except Exception as e:
            return {"status": "error", "source": source, "identifier": identifier, "message": f"An error occurred during download attempt: {e}"}

def main():
    """Command-line interface for fetching PDB files."""
    parser = argparse.ArgumentParser(
        description="Fetch PDB files from internal package data or external sources (RCSB implemented).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "identifier", 
        type=str, 
        help="PDB ID (e.g., 1xyz, 6mpv) or other identifier based on source."
    )
    parser.add_argument(
        "--use-internal",
        action="store_true", # Makes it a flag, True if present
        help="Fetch the PDB file from internal package data instead of external sources."
    )
    parser.add_argument(
        "--source",
        type=str,
        default='rcsb',
        choices=['rcsb', 'ncbi', 'uniprot', 'internal'], # Added 'internal' for clarity
        help="Source database if not using internal. 'rcsb' is implemented, others are placeholders."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./pdb_downloads",
        help="Output directory for downloaded files (required if not using --use-internal)."
    )
    parser.add_argument(
        "--format",
        type=str,
        default='pdb',
        help="File format for download (e.g., 'pdb', 'cif', 'pdb.gz'). Primarily affects RCSB."
    )
    parser.add_argument(
        "--list-internal",
        action="store_true",
        help="List available internal PDB IDs and exit."
    )

    args = parser.parse_args()

    if args.list_internal:
        print("Available internal PDB IDs:")
        ids = list_internal_pdb_ids()
        if ids:
            for pdb_id in ids:
                print(f"- {pdb_id}")
        else:
            print("No internal PDB files found.")
        return # Exit after listing

    # Determine the actual source based on flags
    actual_source = 'internal' if args.use_internal else args.source.lower()
    if actual_source == 'internal':
        use_internal_flag = True
    else:
        use_internal_flag = False

    if not use_internal_flag and not args.output_dir:
        parser.error("--output-dir is required when fetching from external sources (not using --use-internal).")

    print(f"Fetching PDB for identifier: {args.identifier} from source: {actual_source}")

    results = fetch_pdb(
        identifier=args.identifier,
        source=actual_source,
        use_internal=use_internal_flag,
        output_dir=args.output_dir,
        file_format=args.format
    )

    print("--- Fetch Results ---")
    print(f"Status: {results.get('status', 'Unknown').upper()}")
    print(f"Message: {results.get('message', 'No message.')}")
    if results.get("pdb_path"):
        print(f"Output Path: {results['pdb_path']}")
    if results.get("pdb_content"):
        print(f"Content Retrieved (first 100 chars): {results['pdb_content'][:100]}...")

# Example usage block is removed, main() handles CLI.
if __name__ == '__main__':
    main() 