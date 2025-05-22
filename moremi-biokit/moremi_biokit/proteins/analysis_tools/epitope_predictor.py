"""
Predicts B-cell epitopes using the IEDB API.

Provides functionality to query the IEDB B-cell epitope prediction tools.
"""

import iedb
import requests
import pandas as pd
from typing import List, Dict, Any, Optional, Union

# Try to import iedb, but handle if not available for basic type checking
# In a real scenario, iedb would be a dependency.
try:
    import iedb # type: ignore
except ImportError:
    iedb = None # Allows for type checking without runtime error if not installed

IEDB_DEFAULT_METHOD = "Bepipred-2.0" # Bepipred-2.0 is often a good default
VALID_IEDB_BCELL_METHODS: List[str] = [
    "Chou-Fasman", "Emini", "Karplus-Schulz", 
    "Kolaskar-Tongaonkar", "Parker", "Bepipred", "Bepipred-2.0"
]

def _extract_peptide_sequences(peptides: List[Dict[str, Any]]) -> Optional[List[str]]:
    """Extracts peptide sequences from a list of epitope prediction result dictionaries.

    Internal helper function.

    Args:
        peptides: A list of dictionaries, where each dictionary represents an 
                  epitope and is expected to contain a 'Peptide' key.

    Returns:
        A list of peptide sequence strings, or None if the input is None.
    """
    if peptides is None:
        return None
    # Ensure that we handle potential non-dict items gracefully if needed
    return [epitope.get('Peptide', '') for epitope in peptides if isinstance(epitope, dict)]

def _extract_epitope_peptides_from_iedb_df(
    iedb_results_df: pd.DataFrame
) -> List[Dict[str, Any]]:
    """Helper function to parse IEDB DataFrame and extract epitope peptide details.
    
    Args:
        iedb_results_df (pd.DataFrame): DataFrame from iedb.query_bcell_epitope,
            expected columns: 'Position' (0-indexed), 'Residue', 'Score', 'Assignment' ('E' for epitope).
    
    Returns:
        List[Dict[str, Any]]: A list of dictionaries, where each dictionary 
            represents a predicted epitope peptide with keys:
            - "start_position" (int): 1-based start position of the epitope.
            - "end_position" (int): 1-based end position of the epitope.
            - "peptide_sequence" (str): The amino acid sequence of the epitope.
            - "length" (int): Length of the epitope peptide.
            - "average_score" (float): Average prediction score for the residues in this peptide.
    """
    peptides: List[Dict[str, Any]] = []
    if iedb_results_df.empty:
        return peptides

    current_peptide_residues: List[str] = []
    current_peptide_scores: List[float] = []
    peptide_start_pos: Optional[int] = None # 0-indexed

    # Ensure DataFrame is sorted by position if not already guaranteed
    df = iedb_results_df.sort_values(by="Position").reset_index(drop=True)

    for idx, row in df.iterrows():
        is_epitope_residue = row["Assignment"] == "E"
        current_pos = int(row["Position"]) # 0-indexed from IEDB direct output usually

        if is_epitope_residue:
            if peptide_start_pos is None:
                peptide_start_pos = current_pos
            current_peptide_residues.append(str(row["Residue"]))
            current_peptide_scores.append(float(row["Score"]))
        
        # If not an epitope residue OR it's the last residue in the DataFrame
        if (not is_epitope_residue or idx == len(df) - 1) and current_peptide_residues:
            end_pos = current_pos -1 if not is_epitope_residue else current_pos # 0-indexed end
            if peptide_start_pos is not None:
                peptides.append({
                    "start_position": peptide_start_pos + 1, # Convert to 1-based for reporting
                    "end_position": end_pos + 1,       # Convert to 1-based for reporting
                    "peptide_sequence": "".join(current_peptide_residues),
                    "length": len(current_peptide_residues),
                    "average_score": round(sum(current_peptide_scores) / len(current_peptide_scores), 3)
                })
            # Reset for next peptide
            current_peptide_residues = []
            current_peptide_scores = []
            peptide_start_pos = None
            
    return peptides

def get_epitope_sequences_from_prediction(
    predicted_epitopes: List[Dict[str, Any]]
) -> List[str]:
    """Extract just the peptide sequences from a list of predicted epitope dicts.

    Args:
        predicted_epitopes (List[Dict[str, Any]]): A list of epitope dictionaries,
            as returned by _extract_epitope_peptides_from_iedb_df, where each dict 
            is expected to have a "peptide_sequence" key.

    Returns:
        List[str]: A list of epitope peptide sequences.
    """
    if not predicted_epitopes:
        return []
    return [epitope.get("peptide_sequence", "") for epitope in predicted_epitopes if epitope.get("peptide_sequence")]

def predict_bcell_epitopes(
    sequence: str,
    method: str = IEDB_DEFAULT_METHOD,
    window_size: Optional[int] = 9,
) -> Dict[str, Any]:
    """Predict B-cell epitopes in a protein sequence using IEDB analysis tools.

    Args:
        sequence (str): The protein sequence for epitope prediction.
        method (str, optional): The IEDB prediction method to use. 
            Defaults to IEDB_DEFAULT_METHOD ("Bepipred-2.0"). 
            See VALID_IEDB_BCELL_METHODS for available options.
        window_size (Optional[int], optional): Window size for the prediction method. 
            Not all methods use this. If None, IEDB default for the method is used. 
            Defaults to None.

    Returns:
        Dict[str, Any]: A dictionary containing the prediction results:
        - "iedb_raw_results_df" (pd.DataFrame): DataFrame with per-residue scores and assignments.
        - "predicted_epitopes" (List[Dict[str, Any]]): List of identified epitope peptides,
              each with start/end positions, sequence, length, and average score.
        - "epitope_sequences_list" (List[str]): A simple list of the predicted epitope peptide sequences.
        - "epitope_count" (int): Total number of identified epitope peptides.
        - "overall_average_score" (float): Average score of all residues marked as epitopes.
        - "parameters" (Dict[str, Any]): Input parameters used for the prediction.
        On error, returns a dictionary with an "error" key and a descriptive message.
    """
    if not iedb:
        return {"error": "IEDB library is not installed or available. Please install `iedb` to use this function."}

    if not isinstance(sequence, str) or not sequence.strip():
        return {"error": "Input sequence cannot be empty or invalid."}
    cleaned_sequence = sequence.strip().upper()

    if method not in VALID_IEDB_BCELL_METHODS:
        return {"error": f"Invalid IEDB method: '{method}'. Valid options are: {VALID_IEDB_BCELL_METHODS}"}

    parameters_used = {"method": method, "sequence_length": len(cleaned_sequence)}
    if window_size is not None:
        parameters_used["window_size"] = window_size

    try:
        # query_bcell_epitope can take sequence, method, and window_size
        # It typically returns a pandas DataFrame
        
        # TODO: DELETE THIS SECTIONS AFTER SOME REVIEWS
        # call_args = {"method": method, "sequence": cleaned_sequence}
        # if window_size is not None: # Only pass window_size if specified
        #     call_args["window_size"] = window_size
        # raw_iedb_df = iedb.query_bcell_epitope(**call_args) # type: ignore
        
        raw_iedb_df = iedb.query_bcell_epitope(method=method, sequence=sequence, window_size=window_size)

        if not isinstance(raw_iedb_df, pd.DataFrame) or raw_iedb_df.empty:
            # Some methods might return empty if no epitopes or on error, 
            # or sometimes a string message. iedb library behavior can vary.
            return {
                "error": "IEDB query returned no data or an unexpected format.",
                "iedb_raw_results_df": pd.DataFrame(), 
                "predicted_epitopes": [],
                "epitope_sequences_list": [],
                "epitope_count": 0,
                "overall_average_score": 0.0,
                "parameters": parameters_used
            }
        
        # Standardize DataFrame columns if necessary (IEDB usually consistent)
        # Expected columns: Position, Residue, Score, Assignment
        required_cols = ["Position", "Residue", "Score", "Assignment"]
        if not all(col in raw_iedb_df.columns for col in required_cols):
             return {"error": f"IEDB DataFrame missing required columns. Found: {raw_iedb_df.columns.tolist()}", 
                     "details": "Expected columns: 'Position', 'Residue', 'Score', 'Assignment'"}
        
        # Ensure correct types for key columns
        try:
            raw_iedb_df["Position"] = pd.to_numeric(raw_iedb_df["Position"])
            raw_iedb_df["Score"] = pd.to_numeric(raw_iedb_df["Score"])
            raw_iedb_df["Residue"] = raw_iedb_df["Residue"].astype(str)
            raw_iedb_df["Assignment"] = raw_iedb_df["Assignment"].astype(str).str.strip()
        except Exception as e:
            return {"error": f"Error standardizing IEDB DataFrame column types: {e}", "iedb_raw_results_df": raw_iedb_df.to_dict()}

        # Extract peptide segments marked as epitopes
        extracted_peptides = _extract_epitope_peptides_from_iedb_df(raw_iedb_df)
        epitope_sequences = get_epitope_sequences_from_prediction(extracted_peptides)
        
        # Calculate overall average score from epitope residues in the raw dataframe
        epitope_residue_scores = raw_iedb_df[raw_iedb_df["Assignment"] == "E"]["Score"]
        overall_avg_score = round(epitope_residue_scores.mean(), 3) if not epitope_residue_scores.empty else 0.0

        return {
            # "iedb_raw_results_df": raw_iedb_df.to_dict('records'), # Return as list of dicts for easier JSON serialization
            "predicted_epitopes": extracted_peptides,
            "epitope_sequences_list": epitope_sequences,
            "epitope_count": len(extracted_peptides),
            "overall_average_score": overall_avg_score,
            "parameters": parameters_used
        }

    except requests.exceptions.RequestException as req_err:
        return {"error": f"IEDB API request failed: {req_err}", "parameters": parameters_used}
    except ValueError as val_err: # Can be raised by IEDB or DataFrame processing
        return {"error": f"Data processing error: {val_err}", "parameters": parameters_used}
    except Exception as e:
        # Log full error: print(f"Unexpected error in predict_bcell_epitopes: {type(e).__name__} {e}")
        return {"error": f"An unexpected error occurred: {e}", "parameters": parameters_used}

if __name__ == '__main__':
    example_seq = "QGQLVKPAGVALGVTAVAEMLPTGSLVPKPVSVTTAVGSVLAK"
    
    print(f"Predicting epitopes for sequence: {example_seq[:20]}...")
    # Ensure IEDB is installed: pip install iedb
    if iedb is None:
        print("\n*** IEDB library not found. Please install it to run the example: pip install iedb ***")
    else:
        results_default = predict_bcell_epitopes(example_seq)
        print("\n--- Results (Default Method: Bepipred-2.0) ---")
        if "error" in results_default:
            print(f"Error: {results_default['error']}")
        else:
            print(f"Epitope Count: {results_default['epitope_count']}")
            print(f"Overall Average Score: {results_default['overall_average_score']}")
            print(f"Predicted Epitope Sequences: {results_default['epitope_sequences_list']}")
            # print("Detailed Predicted Epitopes:")
            # for ep in results_default['predicted_epitopes']:
            #     print(f"  - Seq: {ep['peptide_sequence']}, Start: {ep['start_position']}, End: {ep['end_position']}, Avg Score: {ep['average_score']}")

        results_emini = predict_bcell_epitopes(example_seq, method="Emini")
        print("\n--- Results (Method: Emini) ---")
        if "error" in results_emini:
            print(f"Error: {results_emini['error']}")
        else:
            print(f"Epitope Count: {results_emini['epitope_count']}")
            print(f"Overall Average Score: {results_emini['overall_average_score']}")
            print(f"Predicted Epitope Sequences: {results_emini['epitope_sequences_list']}") 