"""Utility for fetching Protein Data Bank (PDB) files.

Provides functions to:
- List available internal PDB files bundled with the package.
- Retrieve paths to or content of internal PDB files.
- (Placeholder) Download PDB files from external databases (PDB, NCBI, UniProt).
"""

import os
import importlib.resources
from typing import List, Optional, Dict, Any, Union
import requests
import argparse

# Import for structure prediction
from ..analysis_tools import structure_predictor

# Define the package path for internal PDBs
_INTERNAL_PDB_PKG_PATH = "moremi_biokit.proteins.pdb_targets"

def list_internal_pdb_ids() -> List[str]:
    """List the PDB IDs of the internal PDB files bundled with the package.

    Reads the contents of the `moremi_biokit.proteins.pdb_targets` directory.

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

# --- Helper function to read sequences from file ---
def _read_sequences_from_file(file_path: str) -> List[str]:
    """Reads sequences from a file (plain text or FASTA).

    Args:
        file_path (str): Path to the sequence file.

    Returns:
        List[str]: A list of sequences found in the file. Returns empty if error or no sequences.
    """
    sequences: List[str] = []
    current_seq_lines: List[str] = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq_lines: # Save previous sequence
                        sequences.append("".join(current_seq_lines))
                        current_seq_lines = []
                    # Header line, ignore for now, just start new sequence
                else:
                    # Validate if it looks like an amino acid sequence (simple check)
                    if all(c.upper() in "ACDEFGHIKLMNPQRSTVWY" for c in line.replace("*", "")):
                        current_seq_lines.append(line.upper())
                    else:
                        print(f"Warning: Line in sequence file does not appear to be a valid protein sequence, skipping: {line[:30]}...")
                        if current_seq_lines: # If we were accumulating a sequence, save it and reset
                            sequences.append("".join(current_seq_lines))
                            current_seq_lines = []
            if current_seq_lines: # Add the last sequence in the file
                sequences.append("".join(current_seq_lines))
    except FileNotFoundError:
        print(f"Error: Sequence file not found at {file_path}")
        return []
    except Exception as e:
        print(f"Error reading sequence file {file_path}: {e}")
        return []
    return sequences

# --- Helper function for predicting structure from a single sequence ---
def _fetch_pdb_from_sequence(
    sequence: str,
    output_dir: str,
    api_token: Optional[str],
    project_title_prefix: str,
    filename_prefix: str,
    polling_interval: int,
    max_attempts: int
) -> Dict[str, Any]:
    """Helper to call structure_predictor for a single sequence."""
    if not output_dir:
        return {"status": "error", "source": "prediction", "identifier": sequence[:30]+"...", "message": "output_dir is required for structure prediction."}
    
    print(f"Attempting structure prediction for sequence: {sequence[:30]}... into {output_dir}")
    prediction_result = structure_predictor.predict_structure(
        sequence=sequence,
        api_token=api_token,
        project_title_prefix=project_title_prefix,
        output_directory=output_dir,
        output_pdb_filename_prefix=filename_prefix, # predict_structure appends timestamp and .pdb
        polling_interval=polling_interval,
        max_polling_attempts=max_attempts
    )

    if "error" in prediction_result:
        return {
            "status": "error", 
            "source": "prediction", 
            "identifier": sequence[:30]+"...", 
            "message": prediction_result["error"],
            "details": prediction_result.get("details")
        }
    else:
        return {
            "status": "success",
            "source": "prediction",
            "identifier": sequence[:30]+"...",
            "pdb_path": prediction_result.get("pdb_file_path"),
            "pdb_content": prediction_result.get("pdb_content"),
            "gmqe": prediction_result.get("gmqe"),
            "model_details": prediction_result.get("model_details"),
            "message": prediction_result.get("message", "Structure prediction completed.")
        }

# --- Helper function to handle sequence input (single or list) ---
def _handle_sequence_input(
    target_sequence_input: Union[str, List[str]],
    output_dir: str,
    api_token: Optional[str],
    project_title_prefix: str,
    filename_prefix: str,
    polling_interval: int,
    max_attempts: int
) -> Dict[str, Any]:
    sequence_to_predict = ""
    if isinstance(target_sequence_input, list):
        if not target_sequence_input:
            return {"status": "error", "source": "prediction", "identifier": "empty_list", "message": "Provided sequence list is empty."}
        print(f"Warning: Multiple sequences provided. Predicting structure for the first sequence only: {target_sequence_input[0][:30]}...")
        sequence_to_predict = target_sequence_input[0]
    elif isinstance(target_sequence_input, str):
        sequence_to_predict = target_sequence_input
    else:
        return {"status": "error", "source": "prediction", "identifier": "invalid_type", "message": "Invalid type for sequence input."}

    if not sequence_to_predict.strip():
        return {"status": "error", "source": "prediction", "identifier": "empty_sequence", "message": "Provided sequence is empty."}

    return _fetch_pdb_from_sequence(sequence_to_predict, output_dir, api_token, project_title_prefix, filename_prefix, polling_interval, max_attempts)

# --- Helper function to handle sequence file input ---
def _handle_sequence_file_input(
    file_path: str,
    output_dir: str,
    api_token: Optional[str],
    project_title_prefix: str,
    filename_prefix: str,
    polling_interval: int,
    max_attempts: int
) -> Dict[str, Any]:
    sequences = _read_sequences_from_file(file_path)
    if not sequences:
        return {"status": "error", "source": "prediction", "identifier": file_path, "message": f"No valid sequences found in file: {file_path}."}
    
    print(f"Found {len(sequences)} sequence(s) in {file_path}. Predicting structure for the first valid sequence.")
    # Using the first sequence for prediction
    return _fetch_pdb_from_sequence(sequences[0], output_dir, api_token, project_title_prefix, filename_prefix, polling_interval, max_attempts)

# --- Helper function to handle local PDB file path input ---
def _handle_local_pdb_file(file_path: str) -> Dict[str, Any]:
    """Handles fetching a PDB from a local file path."""
    abs_file_path = os.path.abspath(file_path)
    if not os.path.exists(abs_file_path):
        return {"status": "error", "source": "local", "identifier": file_path, "message": f"Local PDB file not found: {abs_file_path}"}
    
    # Basic check for PDB extension (can be expanded)
    if not file_path.lower().endswith(('.pdb', '.ent', '.cif')):
        print(f"Warning: File '{file_path}' does not have a typical PDB extension (.pdb, .ent, .cif). Attempting to treat as PDB.")

    # Optionally, read content if needed, or just return path
    # For consistency, let's try to read a snippet or return None for pdb_content for local files by default
    # For now, we will not read the content to avoid loading large files into memory by default.
    # The caller can read it if needed using the returned path.
    return {
        "status": "success",
        "source": "local",
        "identifier": file_path,
        "pdb_path": abs_file_path,
        "pdb_content": None, # Content not read by default for local files
        "message": f"Successfully located local PDB file: {abs_file_path}"
    }

# --- Main Fetching Function (Refactored) ---
def fetch_pdb(
    target: Union[str, List[str]], # Can be PDB ID, sequence, list of sequences, path to sequence file, path to local PDB
    input_type: str = 'identifier', # 'identifier', 'sequence', 'sequence_list', 'sequence_file', 'local_pdb_file'
    # Parameters for 'identifier' input_type
    source_db: str = 'rcsb',      # For 'identifier' type: 'rcsb', 'ncbi', 'uniprot'
    use_internal_db: bool = False, # For 'identifier' type
    db_file_format: str = 'pdb',   # Format for DB download (e.g. rcsb)
    # Parameters for 'sequence', 'sequence_list' or 'sequence_file' input_types (passed to structure_predictor)
    swiss_model_api_token: Optional[str] = None,
    prediction_project_title_prefix: str = "pdb_fetcher_prediction",
    prediction_filename_prefix: str = "predicted_pdb_from_seq",
    prediction_polling_interval: int = structure_predictor.DEFAULT_POLLING_INTERVAL_SECONDS,
    prediction_max_attempts: int = structure_predictor.DEFAULT_MAX_POLLING_ATTEMPTS,
    # General parameters
    output_dir: Optional[str] = None # Required for downloads AND predictions. Not for local_pdb_file unless it is to be copied.
) -> Dict[str, Any]:
    """Fetch a PDB file from various sources or predict from sequence.

    Args:
        target (Union[str, List[str]]): The main input. Its meaning depends on `input_type`.
            - If `input_type` is 'identifier': A PDB ID (e.g., "1xyz") or other database ID.
            - If `input_type` is 'sequence': An amino acid sequence string.
            - If `input_type` is 'sequence_list': A list of amino acid sequence strings (only the first is used for prediction).
            - If `input_type` is 'sequence_file': Path to a file containing sequence(s) (FASTA or plain, first sequence used).
            - If `input_type` is 'local_pdb_file': Path to an existing PDB file on the local system.
        input_type (str, optional): Specifies how to interpret `target`.
            Choices: 'identifier', 'sequence', 'sequence_list', 'sequence_file', 'local_pdb_file'. 
            Defaults to 'identifier'.
        source_db (str, optional): Database for 'identifier' type ('rcsb', 'ncbi', 'uniprot'). 
            Defaults to 'rcsb'.
        use_internal_db (bool, optional): If True and `input_type` is 'identifier', use internal PDBs. 
            Defaults to False.
        db_file_format (str, optional): File format for DB downloads (e.g., 'pdb', 'cif'). 
            Defaults to 'pdb'.
        swiss_model_api_token (Optional[str], optional): API token for SWISS-MODEL (for sequence prediction).
        prediction_project_title_prefix (str, optional): Project title prefix for SWISS-MODEL.
        prediction_filename_prefix (str, optional): Filename prefix for predicted PDBs.
        prediction_polling_interval (int, optional): Polling interval for SWISS-MODEL.
        prediction_max_attempts (int, optional): Max polling attempts for SWISS-MODEL.
        output_dir (Optional[str], optional): Directory for downloaded or predicted PDBs. 
            Required if `input_type` leads to a download or prediction. Defaults to None.

    Returns:
        Dict[str, Any]: A dictionary with fetch status, source, path/content, and messages.
            Includes keys like "status" ('success'/'error'), "source", "identifier", 
            "pdb_path", "pdb_content", "message", and prediction-specific details if applicable.
    """
    if input_type == 'identifier':
        if not isinstance(target, str):
            return {"status": "error", "source": source_db or "internal", "identifier": str(target)[:50], "message": "Target for input_type 'identifier' must be a string ID."}
        identifier_str = target
        if use_internal_db:
            pdb_content = read_internal_pdb(identifier_str)
            if pdb_content:
                return {
                    "status": "success", "source": "internal", "identifier": identifier_str,
                    "pdb_path": None, "pdb_content": pdb_content,
                    "message": f"Successfully read internal PDB: {identifier_str}"
                }
            else:
                 pdb_path = get_internal_pdb_path(identifier_str)
                 if pdb_path:
                     return {
                        "status": "success", "source": "internal", "identifier": identifier_str,
                        "pdb_path": pdb_path, "pdb_content": None,
                        "message": f"Found internal PDB path: {pdb_path}."
                     }
                 else:
                    return {"status": "error", "source": "internal", "identifier": identifier_str, "message": f"Internal PDB ID '{identifier_str}' not found."}
        else: # External database download
            if not output_dir:
                return {"status": "error", "source": source_db, "identifier": identifier_str, "message": "output_dir is required for external database downloads."}
            os.makedirs(output_dir, exist_ok=True)
            file_path: Optional[str] = None
            try:
                if source_db.lower() == 'rcsb':
                    file_path = download_pdb_from_rcsb(identifier_str, output_dir, db_file_format)
                elif source_db.lower() == 'ncbi':
                    file_path = download_structure_from_ncbi(identifier_str, output_dir) # Placeholder
                elif source_db.lower() == 'uniprot':
                    file_path = download_structure_from_uniprot(identifier_str, output_dir) # Placeholder
                else:
                    return {"status": "error", "source": source_db, "identifier": identifier_str, "message": f"Unsupported external source_db: {source_db}."}
                
                if file_path:
                     return {"status": "success", "source": source_db, "identifier": identifier_str, "pdb_path": file_path, "message": f"Downloaded from {source_db} to {file_path}"}
                else:
                    # Error message should come from download_xxx functions if they fail properly
                    return {"status": "error", "source": source_db, "identifier": identifier_str, "message": f"Download failed from {source_db} for {identifier_str}."}
            except NotImplementedError as nie:
                return {"status": "error", "source": source_db, "identifier": identifier_str, "message": f"Download function for {source_db} is not implemented: {nie}"}
            except Exception as e:
                return {"status": "error", "source": source_db, "identifier": identifier_str, "message": f"Error during download from {source_db}: {e}"}

    elif input_type == 'sequence':
        if not isinstance(target, str):
            return {"status": "error", "source": "prediction", "identifier": str(target)[:50], "message": "Target for input_type 'sequence' must be a string."}
        if not output_dir:
            return {"status": "error", "source": "prediction", "identifier": target[:30]+"...", "message": "output_dir is required for structure prediction from sequence."}
        return _handle_sequence_input(target, output_dir, swiss_model_api_token, prediction_project_title_prefix, prediction_filename_prefix, prediction_polling_interval, prediction_max_attempts)
    
    elif input_type == 'sequence_list':
        if not isinstance(target, list):
            return {"status": "error", "source": "prediction", "identifier": "invalid_list_input", "message": "Target for input_type 'sequence_list' must be a list of strings."}
        if not output_dir:
            return {"status": "error", "source": "prediction", "identifier": "list_input", "message": "output_dir is required for structure prediction from sequence list."}
        return _handle_sequence_input(target, output_dir, swiss_model_api_token, prediction_project_title_prefix, prediction_filename_prefix, prediction_polling_interval, prediction_max_attempts)

    elif input_type == 'sequence_file':
        if not isinstance(target, str):
             return {"status": "error", "source": "prediction", "identifier": str(target)[:50], "message": "Target for input_type 'sequence_file' must be a file path string."}
        if not output_dir:
            return {"status": "error", "source": "prediction", "identifier": target, "message": "output_dir is required for structure prediction from sequence file."}
        return _handle_sequence_file_input(target, output_dir, swiss_model_api_token, prediction_project_title_prefix, prediction_filename_prefix, prediction_polling_interval, prediction_max_attempts)

    elif input_type == 'local_pdb_file':
        if not isinstance(target, str):
            return {"status": "error", "source": "local", "identifier": str(target)[:50], "message": "Target for input_type 'local_pdb_file' must be a file path string."}
        return _handle_local_pdb_file(target)

    else:
        return {"status": "error", "source": "unknown", "identifier": str(target)[:50], "message": f"Invalid input_type: '{input_type}'."}

def main():
    """Command-line interface for fetching PDB files or predicting from sequence."""
    parser = argparse.ArgumentParser(
        description="Fetch PDB files from various sources or predict structure from sequence using SWISS-MODEL.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "target", 
        help="The target input: PDB ID, sequence, path to sequence file, or path to local PDB file. Its interpretation depends on --input-type."
    )
    parser.add_argument(
        "--input-type",
        type=str,
        default='identifier',
        choices=['identifier', 'sequence', 'sequence_file', 'local_pdb_file'], # sequence_list not practical for CLI target, handled if API user passes list
        help="Specify how to interpret the 'target' argument."
    )
    parser.add_argument(
        "--use-internal-db",
        action="store_true", 
        help="For 'identifier' input_type: fetch from internal package data."
    )
    parser.add_argument(
        "--source-db",
        type=str,
        default='rcsb',
        choices=['rcsb', 'ncbi', 'uniprot'], 
        help="For 'identifier' input_type (and not --use-internal-db): external database to use."
    )
    parser.add_argument(
        "--db-file-format",
        type=str,
        default='pdb',
        help="For 'identifier' input_type (RCSB): file format for download (e.g., 'pdb', 'cif')."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./pdb_fetched_or_predicted",
        help="Output directory for downloaded or predicted files. Always created if it doesn't exist."
    )
    # Prediction specific arguments
    parser.add_argument(
        "--sm-token", 
        type=str, 
        default=os.environ.get("SWISS_MODEL_API_TOKEN"), 
        help="SWISS-MODEL API token for sequence prediction. Reads from SWISS_MODEL_API_TOKEN env var if not set."
    )
    parser.add_argument(
        "--sm-title-prefix", type=str, default="cli_pdb_fetcher_prediction", 
        help="Project title prefix for SWISS-MODEL submission."
    )
    parser.add_argument(
        "--sm-filename-prefix", type=str, default="predicted_from_seq",
        help="Filename prefix for PDB files predicted by SWISS-MODEL."
    )
    parser.add_argument(
        "--sm-polling-interval", type=int, default=structure_predictor.DEFAULT_POLLING_INTERVAL_SECONDS,
        help="Polling interval for SWISS-MODEL job status."
    )
    parser.add_argument(
        "--sm-max-attempts", type=int, default=structure_predictor.DEFAULT_MAX_POLLING_ATTEMPTS,
        help="Max polling attempts for SWISS-MODEL job."
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
        return

    # Create output_dir if it doesn't exist, as it's used by most branches
    if args.input_type in ['identifier', 'sequence', 'sequence_file'] and not args.use_internal_db:
        if not args.output_dir:
            parser.error("--output-dir is required for this operation unless using --use-internal-db or --input-type local_pdb_file.")
        os.makedirs(args.output_dir, exist_ok=True)
    elif args.input_type == 'local_pdb_file':
        # output_dir is not strictly needed for local_pdb_file unless we plan to copy it, but let's be consistent if provided
        if args.output_dir:
            os.makedirs(args.output_dir, exist_ok=True)
    
    # For sequence prediction, check token
    if args.input_type in ['sequence', 'sequence_file'] and not args.sm_token and not structure_predictor.DEFAULT_SWISS_MODEL_API_TOKENS:
         parser.error("SWISS-MODEL API token (--sm-token or SWISS_MODEL_API_TOKEN env var) is required for structure prediction if no default tokens are configured in structure_predictor.py.")
    elif args.input_type in ['sequence', 'sequence_file'] and not args.sm_token:
        print("Warning: SWISS-MODEL API token not provided. Attempting to use default tokens from structure_predictor.py (if available). This is not recommended for production.")


    print(f"Fetching PDB for target: {args.target[:60]}{'...' if len(args.target) > 60 else ''}")
    print(f"Input Type: {args.input_type}, Output Directory: {args.output_dir}")

    results = fetch_pdb(
        target=args.target,
        input_type=args.input_type,
        source_db=args.source_db,
        use_internal_db=args.use_internal_db,
        db_file_format=args.db_file_format,
        swiss_model_api_token=args.sm_token,
        prediction_project_title_prefix=args.sm_title_prefix,
        prediction_filename_prefix=args.sm_filename_prefix,
        prediction_polling_interval=args.sm_polling_interval,
        prediction_max_attempts=args.sm_max_attempts,
        output_dir=args.output_dir
    )

    print("--- Fetch Results ---")
    print(f"Status: {results.get('status', 'Unknown').upper()}")
    print(f"Source: {results.get('source', 'N/A')}")
    print(f"Identifier/Target processed: {str(results.get('identifier', 'N/A'))[:100]}")
    print(f"Message: {results.get('message', 'No message.')}")
    if results.get("pdb_path"):
        print(f"Output PDB Path: {results['pdb_path']}")
    if results.get("pdb_content") and results.get("source") == "internal": # Only show for internal for brevity
        print(f"Internal PDB Content (first 100 chars): {results['pdb_content'][:100]}...")
    if results.get("status") == "success" and results.get("source") == "prediction":
        print("Prediction Details:")
        if "gmqe" in results: print(f"  GMQE: {results['gmqe']}")
        if "model_details" in results:
            for k,v in results["model_details"].items():
                print(f"  {k.replace('_', ' ').capitalize()}: {v}")
    elif results.get("status") == "error" and "details" in results:
        print(f"Error Details: {results['details']}")

if __name__ == '__main__':
    main() 