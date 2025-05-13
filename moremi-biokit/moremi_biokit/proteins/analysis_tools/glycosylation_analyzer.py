"""
Predicts potential N-linked and O-linked glycosylation sites in protein sequences.

Uses regular expressions for N-linked sites (N-X-S/T sequon, X != P) and heuristics 
based on surrounding residues for O-linked sites (on S or T).
"""

import re
from typing import Dict, Any, Tuple, List, Set

def _get_n_glyc_confidence(sequence: str, site_index: int, sequon_length: int = 3) -> str:
    """Estimate confidence of N-glycosylation prediction based on context.

    Args:
        sequence (str): The full protein sequence (uppercase).
        site_index (int): The 0-based starting index of the N-glycosylation sequon (N).
        sequon_length (int): The length of the sequon (e.g., 3 for N-X-S/T).

    Returns:
        str: Confidence level ('High', 'Medium', 'Low', 'Very Low').
    """
    # Check if the full sequon is within sequence bounds
    if site_index + sequon_length > len(sequence):
        return "Low" # Should not happen if called correctly after finding a match

    # Proline at X position (N-P-S/T) strongly inhibits N-glycosylation
    if sequence[site_index + 1] == 'P':
        return "Very Low"

    # Contextual analysis based on surrounding residues (simplified)
    # These sets can be refined based on literature for more accuracy.
    favorable_context_residues: Set[str] = set('DQSTVC') # Asp, Gln, Ser, Thr, Val, Cys
    unfavorable_context_residues: Set[str] = set('WFYM') # Trp, Phe, Tyr, Met (bulky, hydrophobic)

    # Define a window around the N-X-S/T sequon (e.g., +/- 2 residues from N)
    window_start = max(0, site_index - 2)
    window_end = min(len(sequence), site_index + sequon_length + 2)
    context_window = sequence[window_start:site_index] + sequence[site_index + sequon_length : window_end]

    favorable_count = sum(1 for aa in context_window if aa in favorable_context_residues)
    unfavorable_count = sum(1 for aa in context_window if aa in unfavorable_context_residues)

    if favorable_count > unfavorable_count + 1: # Needs more favorable than unfavorable
        return "High"
    elif favorable_count >= unfavorable_count:
        return "Medium"
    else:
        return "Low"

def predict_glycosylation(sequence: str, o_glyc_window_size: int = 5, o_glyc_favorable_threshold: int = 2) -> Dict[str, Any]:
    """Predict potential N-glycosylation and O-glycosylation sites.

    N-glycosylation sites are identified using the N-X-S/T sequon (X != P).
    O-glycosylation prediction is a simplified model based on S/T residues in a 
    favorable context (e.g., near Pro, Ser, Thr, Asp).

    Args:
        sequence (str): The protein sequence to analyze.
        o_glyc_window_size (int, optional): Window size for O-glycosylation context analysis. 
            Defaults to 5.
        o_glyc_favorable_threshold (int, optional): Minimum count of favorable residues 
            in the window for O-glycosylation. Defaults to 2.

    Returns:
        Dict[str, Any]: A dictionary with glycosylation prediction results:
        - "n_glycosylation" (Dict[str, Any]):
            - "sites" (List[Dict[str, Any]]): List of predicted N-glyc sites, each with:
                - "position" (int): 1-based position of N.
                - "sequon" (str): The N-X-S/T sequon found.
                - "confidence" (str): Estimated confidence ('High', 'Medium', 'Low', 'Very Low').
            - "count" (int): Total number of N-glyc sites.
        - "o_glycosylation" (Dict[str, Any]):
            - "sites" (List[Dict[str, Any]]): List of predicted O-glyc sites, each with:
                - "position" (int): 1-based position of S/T.
                - "residue" (str): The S or T residue.
                - "context_window" (str): Amino acid context around the site.
            - "count" (int): Total number of O-glyc sites.
        Returns an {"error": "message"} dict if the input sequence is invalid.
    """
    if not isinstance(sequence, str) or not sequence.strip():
        return {"error": "Input sequence cannot be empty or invalid."}
    
    cleaned_sequence = sequence.strip().upper()
    seq_len = len(cleaned_sequence)

    # N-glycosylation prediction
    n_glyc_pattern = r'N[^P][ST]' # N-X-S/T where X is not P
    predicted_n_sites: List[Dict[str, Any]] = []
    for match in re.finditer(n_glyc_pattern, cleaned_sequence):
        site_start_idx = match.start()
        sequon = match.group(0)
        confidence = _get_n_glyc_confidence(cleaned_sequence, site_start_idx, len(sequon))
        predicted_n_sites.append({
            "position": site_start_idx + 1, # 1-based
            "sequon": sequon,
            "confidence": confidence
        })

    # O-glycosylation prediction (simplified)
    # Amino acids often found near O-glycosylation sites (can be refined)
    o_glyc_favorable_context: Set[str] = set('PGSTDR') # Pro, Gly, Ser, Thr, Asp, Arg are common
    predicted_o_sites: List[Dict[str, Any]] = []

    for i in range(seq_len):
        current_residue = cleaned_sequence[i]
        if current_residue in ('S', 'T'):
            # Define context window around the S/T residue
            # Example: window of size 5 means 2 residues on each side + S/T itself
            half_window = (o_glyc_window_size -1) // 2
            context_start = max(0, i - half_window)
            context_end = min(seq_len, i + half_window + 1)
            context_seq_window = cleaned_sequence[context_start:context_end]

            favorable_residues_in_window = sum(
                1 for aa in context_seq_window if aa in o_glyc_favorable_context
            )

            if favorable_residues_in_window >= o_glyc_favorable_threshold:
                predicted_o_sites.append({
                    "position": i + 1, # 1-based
                    "residue": current_residue,
                    "context_window": context_seq_window
                })
                
    return {
        "n_glycosylation": {
            "sites": predicted_n_sites,
            "count": len(predicted_n_sites)
        },
        "o_glycosylation": {
            "sites": predicted_o_sites,
            "count": len(predicted_o_sites)
        }
    } 