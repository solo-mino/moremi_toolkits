"""
Local PDB FASTA Database Parser for Moremi Biokit.

This module provides functions to interact with a local copy of the RCSB PDB
sequence database stored in FASTA format (rcsb_pdb_sequence_db.txt).
It uses `pyfaidx` for efficient indexing and retrieval of sequences.

The main functionalities are:
- Retrieving a specific sequence by PDB ID and chain ID.
- Retrieving all sequences associated with a given PDB ID.

The module expects 'rcsb_pdb_sequence_db.txt' to be available as a package resource
within the 'moremi_biokit.connectors' path.
"""

# TODO: CHANGE THIS FUNCTIONALITY TO USE NEW DOCKERIZED MOREMI-SERVICE API
import importlib.resources
from pyfaidx import Fasta, FastaNotFoundError, FastaIndexingError
from typing import Optional, Dict, Any, List
import os # For checking if .fai index needs rebuild potentially, though pyfaidx handles it.
import traceback

# --- Module-level Configuration ---
FASTA_DB_FILENAME = "rcsb_pdb_sequence_db.txt"
_fasta_db_instance: Optional[Fasta] = None
_fasta_db_path_str: Optional[str] = None


def _get_and_initialize_fasta_db() -> Optional[Fasta]:
    """
    Initializes and returns the pyfaidx.Fasta object for the local PDB sequence database.
    Uses a module-level cache to avoid re-initializing on every call.
    The FASTA DB file is located using importlib.resources.
    """
    global _fasta_db_instance
    global _fasta_db_path_str

    if _fasta_db_instance is not None:
        return _fasta_db_instance

    try:
        # importlib.resources.files() returns a Traversable object.
        # We need its string path for pyfaidx.
        package_path = importlib.resources.files('moremi_biokit.connectors')
        db_file_traversable = package_path.joinpath(FASTA_DB_FILENAME)
        
        if not db_file_traversable.is_file():
            print(f"Error: Local PDB FASTA DB file '{FASTA_DB_FILENAME}' not found at {db_file_traversable}.")
            print("Please ensure the database file is included in the package and accessible via importlib.resources.")
            return None
        
        _fasta_db_path_str = str(db_file_traversable)
        print(f"Initializing local PDB FASTA DB from: {_fasta_db_path_str}")
        
        # pyfaidx will automatically create/use an index file (e.g., rcsb_pdb_sequence_db.txt.fai)
        # in the same directory if it has write permissions, or in a cache directory.
        # read_long_names=True is crucial for getting the full FASTA header.
        _fasta_db_instance = Fasta(_fasta_db_path_str, read_long_names=True)
        
        # It's good practice to check if keys were loaded, indicating successful indexing.
        if not len(_fasta_db_instance.keys()):
             print(f"Warning: Local PDB FASTA DB initialized from '{_fasta_db_path_str}' but contains no sequences. Check DB integrity or .fai index.")
             # _fasta_db_instance will be an empty Fasta object, subsequent calls will find no keys.
        else:
            print(f"Local PDB FASTA DB initialized successfully. Sequences available: {len(_fasta_db_instance.keys())}")
        
        return _fasta_db_instance

    except FastaNotFoundError:
        print(f"Error: pyfaidx could not find the FASTA DB file at '{_fasta_db_path_str}'. This should have been caught by is_file().")
        _fasta_db_instance = None 
        return None
    except FastaIndexingError as e:
        print(f"Error: pyfaidx failed to index the FASTA DB file at '{_fasta_db_path_str}': {e}")
        print("This might be due to a corrupted FASTA file or index, or permissions issues for creating/reading the .fai index.")
        _fasta_db_instance = None
        return None
    except ImportError as e: # Handles if 'moremi_biokit.connectors' is not a valid package path
        print(f"Error resolving FASTA DB path with importlib.resources: {e}. Is 'moremi_biokit.connectors' a proper package containing the DB?")
        _fasta_db_instance = None
        return None
    except Exception as e:
        print(f"An unexpected error occurred while initializing FASTA DB from '{_fasta_db_path_str}':")
        print(traceback.format_exc())
        _fasta_db_instance = None
        return None


def _parse_fasta_record(record_key: str, faidx_record: Any) -> Dict[str, Any]:
    """
    Helper function to parse a pyfaidx FastaRecord into the desired dictionary format.
    
    Args:
        record_key (str): The key used to retrieve the record (e.g., "1xyz_A").
                          This key should already be in the desired PDBID_CHAIN format.
        faidx_record (pyfaidx.FastaRecord): The record object from pyfaidx.

    Returns:
        Dict[str, Any]: A dictionary containing structured sequence information.
    """
    # record_key is assumed to be already normalized (e.g. 1xyz_A) by the calling function.
    pdb_id_part, chain_id_part = record_key.split('_', 1) if '_' in record_key else (record_key, None)
    
    if chain_id_part is None: # Handle cases where key might be just PDB ID (though less expected for this DB)
        print(f"Warning: FASTA record key '{record_key}' does not contain an underscore for chain ID separation. Processing as PDB ID only.")
        pdb_id_part = record_key # The whole key is the PDB ID
        # chain_id_part remains None

    return {
        "id": record_key,  # The normalized identifier (e.g., "1xyz_A")
        "pdb_id": pdb_id_part, # PDB ID part (already lowercase from normalization)
        "chain_id": chain_id_part, # Chain ID part (already uppercase from normalization)
        "sequence": str(faidx_record),
        "length": len(faidx_record),
        "description": faidx_record.long_name # Full FASTA header line after '>'
    }


def get_sequence_by_pdb_chain_id(pdb_id: str, chain_id: str) -> Optional[Dict[str, Any]]:
    """
    Retrieves sequence details for a specific PDB ID and chain ID from the local FASTA DB.

    Args:
        pdb_id (str): The PDB ID (case-insensitive).
        chain_id (str): The chain identifier (case-insensitive).

    Returns:
        Optional[Dict[str, Any]]: A dictionary with sequence details if found, 
                                  otherwise None. Structure:
                                  {"id", "pdb_id", "chain_id", "sequence", "length", "description"}
    """
    if not pdb_id or not chain_id:
        print("Error: PDB ID and chain ID must be provided.")
        return None

    norm_pdb_id = pdb_id.lower()
    norm_chain_id = chain_id.upper()
    # The keys in the FASTA DB are expected to be like '100d_A'
    query_key = f"{norm_pdb_id}_{norm_chain_id}"

    fa_db = _get_and_initialize_fasta_db()
    if not fa_db:
        print("Local FASTA DB not available for get_sequence_by_pdb_chain_id.")
        return None

    try:
        if query_key in fa_db:
            record = fa_db[query_key]
            # Pass the already normalized query_key for parsing, as it matches the desired output 'id' format.
            return _parse_fasta_record(query_key, record)
        else:
            # The provided FASTA example shows keys like '100d_A' (PDB lowercase).
            # If the database could have PDB IDs in uppercase (e.g. '1XYZ_A'), pyfaidx keys would reflect that.
            # However, our normalization standard is lowercase PDB, uppercase Chain for query_key.
            print(f"Sequence for ID '{query_key}' not found in local FASTA DB.")
            return None
    except KeyError: # Should be caught by 'in fa_db' but as a safeguard.
        print(f"Sequence for ID '{query_key}' not found in local FASTA DB (KeyError). This is unexpected after 'in' check.")
        return None
    except Exception as e:
        print(f"Error fetching sequence for '{query_key}' from local FASTA DB:")
        print(traceback.format_exc())
        return None


def get_all_sequences_by_pdb_id(pdb_id: str) -> List[Dict[str, Any]]:
    """
    Retrieves all sequence details for a given PDB ID from the local FASTA DB.

    Args:
        pdb_id (str): The PDB ID (case-insensitive).

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, each containing sequence details.
                              Returns an empty list if the PDB ID is not found or an error occurs.
    """
    if not pdb_id:
        print("Error: PDB ID must be provided.")
        return []
        
    norm_pdb_id = pdb_id.lower()
    results: List[Dict[str, Any]] = []
    
    fa_db = _get_and_initialize_fasta_db()
    if not fa_db:
        print("Local FASTA DB not available for get_all_sequences_by_pdb_id.")
        return results

    # Iterate through all keys in the Fasta object.
    # Keys are typically like 'pdbid_CHAIN' from the FASTA headers.
    # We need to match the 'pdbid' part case-insensitively (though we store it lowercase).
    prefix_to_match = f"{norm_pdb_id}_"

    for record_key_from_db in fa_db.keys():
        # record_key_from_db could be '1xyz_A', '1XyZ_b', etc.
        # We want to match norm_pdb_id (e.g., '1xyz') with the start of record_key_from_db.
        if record_key_from_db.lower().startswith(prefix_to_match):
            # Further check: ensure the part before the first '_' actually IS the PDB ID.
            key_pdb_part = record_key_from_db.split('_', 1)[0]
            if key_pdb_part.lower() == norm_pdb_id:
                try:
                    record = fa_db[record_key_from_db]
                    # For _parse_fasta_record, we need to provide a key that represents the
                    # desired output "id" format: lowercase PDB, uppercase Chain.
                    # The chain ID comes from record_key_from_db after the underscore.
                    _, chain_part_from_key = record_key_from_db.split('_',1)
                    normalized_output_id = f"{norm_pdb_id}_{chain_part_from_key.upper()}"
                    
                    parsed_data = _parse_fasta_record(normalized_output_id, record)
                    results.append(parsed_data)
                except KeyError: # Should not happen if key is from fa_db.keys()
                    print(f"Internal Warning: Key '{record_key_from_db}' from fa_db.keys() not found during record fetching.")
                except Exception as e:
                    print(f"Error processing record for key '{record_key_from_db}':")
                    print(traceback.format_exc())
                    
    if not results:
        print(f"No sequences found specifically matching PDB ID prefix '{prefix_to_match}' in local FASTA DB.")
    return results

# Example Usage (for testing this module directly):
if __name__ == '__main__':
    print("--- Testing Local PDB FASTA Parser ---  (Ensure rcsb_pdb_sequence_db.txt is accessible)")

    # Test 1: Initialize DB (will print status)
    db = _get_and_initialize_fasta_db()
    if not db:
        print("DB initialization failed. Exiting tests.")
        exit()

    # Test 2: Fetch a specific sequence (use a known ID from your example if DB is small, else generic)
    # Assuming '101m_A' is in your rcsb_pdb_sequence_db.txt as per your example
    pdb_to_test = "101m" # PDB ID from user example
    chain_to_test = "A"  # Chain from user example
    print(f"\nFetching sequence for {pdb_to_test}_{chain_to_test}...")
    seq_details = get_sequence_by_pdb_chain_id(pdb_to_test, chain_to_test)
    if seq_details:
        print(f"Details for {pdb_to_test}_{chain_to_test}:")
        print(f"  ID: {seq_details['id']}")
        print(f"  PDB: {seq_details['pdb_id']}")
        print(f"  Chain: {seq_details['chain_id']}")
        print(f"  Length: {seq_details['length']}")
        print(f"  Description: {seq_details['description']}")
        print(f"  Sequence: {seq_details['sequence'][:30]}... (first 30 chars)")
    else:
        print(f"Could not fetch details for {pdb_to_test}_{chain_to_test}.")

    # Test 3: Fetch all sequences for a PDB ID
    print(f"\nFetching all sequences for PDB ID {pdb_to_test}...")
    all_seqs = get_all_sequences_by_pdb_id(pdb_to_test)
    if all_seqs:
        print(f"Found {len(all_seqs)} sequences for PDB ID {pdb_to_test}:")
        for item in all_seqs:
            print(f"  ID: {item['id']}, Length: {item['length']}, Desc: {item['description']}")
    else:
        print(f"No sequences found for PDB ID {pdb_to_test}.")

    # Test 4: Non-existent PDB ID
    print("\nFetching sequence for non-existent PDB ID '0xxx_Z'...")
    non_existent_seq = get_sequence_by_pdb_chain_id("0xxx", "Z")
    if not non_existent_seq:
        print("Correctly returned None for non-existent PDB ID '0xxx_Z'.")

    print("\nFetching all sequences for non-existent PDB ID '0xxx'...")
    non_existent_all_seqs = get_all_sequences_by_pdb_id("0xxx")
    if not non_existent_all_seqs: # Expecting empty list
        print("Correctly returned empty list for non-existent PDB ID '0xxx'.")
        
    # Test 5: DB path if file is missing (manual test by temporarily renaming the DB file)
    # print("\n--- End of Tests ---") 