"""
Predicts antibody stability based on amino acid composition.

This module provides a simplified model to estimate antibody stability,
represented by a stability score and an estimated melting temperature (Tm).
It considers the contribution of generally stabilizing and destabilizing residues.
"""

from typing import Dict, Any, Set

def predict_stability(sequence: str) -> Dict[str, Any]:
    """Predict protein stability based on amino acid composition.

    This function provides a simplified model to estimate protein stability,
    inspired by general principles of amino acid contributions to stability.
    It calculates a stability score and estimates a melting temperature.

    Args:
        sequence (str): The amino acid sequence of the protein.

    Returns:
        Dict[str, Any]: A dictionary containing the stability prediction results
            under the key "stability_result". This includes:
            - "melting_temperature_celsius" (float): Estimated melting temperature in °C.
            - "normalized_stability_score" (float): Calculated stability score (normalized).
            - "details" (Dict[str, int]): 
                - "stabilizing_residues_count" (int): Count of stabilizing residues.
                - "destabilizing_residues_count" (int): Count of destabilizing residues.
                - "sequence_length" (int): Length of the input sequence.
                - "unrecognized_residues_count" (int): Count of unrecognized residues.
            Returns an {"error": "message"} dictionary if the sequence is empty or invalid.
    """
    if not isinstance(sequence, str) or not sequence.strip():
        return {"error": "Input sequence cannot be empty or invalid."}
    
    cleaned_sequence = sequence.strip().upper() # Standardize to uppercase
    seq_len = len(cleaned_sequence)

    if seq_len == 0:
        return {"error": "Input sequence cannot be empty after stripping whitespace."}

    # Amino acids generally considered to contribute to stability
    stabilizing_residues: Set[str] = set('CFILLMPVWY') # Added P, V, W, Y for broader coverage
    # Amino acids generally considered to be destabilizing or neutral/context-dependent
    destabilizing_residues: Set[str] = set('DEGHKNQRST') # Added A, S, T (often neutral but can be destabilizing)
                                                    # G, K, N, Q, R more clearly destabilizing in many contexts

    stability_score: float = 0.0
    num_stabilizing: int = 0
    num_destabilizing: int = 0

    for aa in cleaned_sequence:
        if aa in stabilizing_residues:
            stability_score += 1.0
            num_stabilizing += 1
        elif aa in destabilizing_residues:
            stability_score -= 0.5 # Penalize destabilizing residues
            num_destabilizing += 1
        # Amino acids not in either set (e.g., non-standard) are ignored

    # Normalize score by sequence length
    # Avoid division by zero, though already checked by seq_len == 0
    normalized_stability_score = stability_score / seq_len if seq_len > 0 else 0.0

    # Estimate melting temperature (Tm)
    # This is a highly simplified model. Real Tm prediction is complex.
    # Base Tm (e.g., 60°C) adjusted by the normalized stability score.
    # The scaling factor (e.g., 20) determines sensitivity.
    estimated_melting_temp = 60.0 + (normalized_stability_score * 20.0)

    stability_result: Dict[str, Any] = {
        "melting_temperature_celsius": round(estimated_melting_temp, 2),
        "normalized_stability_score": round(normalized_stability_score, 4),
        "details": {
            "stabilizing_residues_count": num_stabilizing,
            "destabilizing_residues_count": num_destabilizing,
            "sequence_length": seq_len,
            "unrecognized_residues_count": seq_len - (num_stabilizing + num_destabilizing)
        }
    }
    
    return {"stability_result": stability_result} 