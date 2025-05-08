"""Module for predicting protein aggregation propensity."""

import numpy as np
from typing import Dict, List, Tuple, Any

# Amino acid scales (can be moved to a separate constants module if used elsewhere)
HYDROPHOBICITY_KYTE_DOOLITTLE: Dict[str, float] = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

BETA_SHEET_PROPENSITY_CHOU_FASMAN: Dict[str, float] = {
    'A': 0.83, 'R': 0.93, 'N': 0.65, 'D': 0.54, 'C': 1.19,
    'Q': 1.10, 'E': 0.37, 'G': 0.75, 'H': 0.87, 'I': 1.60,
    'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
    'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70
}

def predict_aggregation(
    sequence: str, 
    window_size: int = 7, 
    aggregation_threshold: float = 1.0,
    hydrophobicity_weight: float = 0.6,
    beta_propensity_weight: float = 0.4
) -> Dict[str, Any]:
    """Predict protein aggregation propensity using a sliding window approach.

    This method uses a simplified model inspired by algorithms like TANGO,
    combining hydrophobicity and beta-sheet propensity scores to identify
    potential aggregation-prone regions (APRs).

    Args:
        sequence (str): The amino acid sequence of the protein.
        window_size (int, optional): The size of the sliding window to analyze 
            regions. Defaults to 7.
        aggregation_threshold (float, optional): The score threshold above which a 
            window is considered aggregation-prone. Defaults to 1.0.
        hydrophobicity_weight (float, optional): The weight for the hydrophobicity 
            score component. Defaults to 0.6.
        beta_propensity_weight (float, optional): The weight for the beta-sheet 
            propensity score component. Defaults to 0.4.

    Returns:
        Dict[str, Any]: A dictionary containing aggregation prediction results:
            - "aggregation_result" (Dict[str, Any]):
                - "overall_aggregation_propensity" (str): "Low", "Moderate", or "High".
                - "aggregation_prone_regions" (List[Dict[str, Any]]): A list of 
                  identified APRs, each with "region_coordinates" (Tuple[int, int]) 
                  and "sequence_fragment" (str).
                - "average_aggregation_score" (float): The average score over all windows.
            - "identified_regions_count" (int): The number of merged APRs.
            Returns an {"error": "message"} dictionary for invalid input.
    """
    if not isinstance(sequence, str) or not sequence.strip():
        return {"error": "Input sequence cannot be empty or invalid."}
    cleaned_sequence = sequence.strip().upper()
    seq_len = len(cleaned_sequence)

    if seq_len < window_size:
        # Return a structure consistent with successful runs but indicating the issue.
        return {
            "error": f"Sequence length ({seq_len}) is less than window size ({window_size}).",
            "aggregation_result": {
                 "overall_aggregation_propensity": "Undetermined",
                 "aggregation_prone_regions": [],
                 "average_aggregation_score": 0.0
            },
            "identified_regions_count": 0
        }

    aggregation_scores: List[float] = []
    raw_aggregation_regions: List[Tuple[int, int]] = [] # (start_idx, end_idx_exclusive)

    for i in range(seq_len - window_size + 1):
        window_sequence = cleaned_sequence[i : i + window_size]
        
        current_hydrophobicity_score = sum(
            HYDROPHOBICITY_KYTE_DOOLITTLE.get(aa, 0.0) for aa in window_sequence
        ) / window_size
        
        current_beta_score = sum(
            BETA_SHEET_PROPENSITY_CHOU_FASMAN.get(aa, 0.0) for aa in window_sequence
        ) / window_size

        combined_score = (
            current_hydrophobicity_score * hydrophobicity_weight + 
            current_beta_score * beta_propensity_weight
        )
        aggregation_scores.append(combined_score)

        if combined_score > aggregation_threshold:
            raw_aggregation_regions.append((i, i + window_size))

    # Merge overlapping/adjacent regions
    merged_regions: List[Tuple[int, int]] = []
    if raw_aggregation_regions:
        # Sort by start position, then by end position
        sorted_regions = sorted(raw_aggregation_regions, key=lambda x: (x[0], x[1]))
        current_merged_start, current_merged_end = sorted_regions[0]

        for i in range(1, len(sorted_regions)):
            next_start, next_end = sorted_regions[i]
            if next_start <= current_merged_end: # Overlap or adjacent
                current_merged_end = max(current_merged_end, next_end)
            else: # Gap found, start a new merged region
                merged_regions.append((current_merged_start, current_merged_end))
                current_merged_start, current_merged_end = next_start, next_end
        merged_regions.append((current_merged_start, current_merged_end)) # Add the last region

    # Determine overall aggregation propensity based on average score
    average_score = float(np.mean(aggregation_scores)) if aggregation_scores else 0.0
    
    propensity_label: str
    if average_score > aggregation_threshold:
        propensity_label = "High"
    elif average_score > aggregation_threshold * 0.7: # Heuristic for moderate
        propensity_label = "Moderate"
    else:
        propensity_label = "Low"
    
    formatted_apr_list: List[Dict[str, Any]] = []
    if merged_regions:
        for start_idx, end_idx in merged_regions:
            formatted_apr_list.append({
                "region_coordinates": (start_idx + 1, end_idx), # 1-indexed for reporting
                "sequence_fragment": cleaned_sequence[start_idx:end_idx]
            })
    
    aggregation_output: Dict[str, Any] = {
        "overall_aggregation_propensity": propensity_label,
        "aggregation_prone_regions": formatted_apr_list,
        "average_aggregation_score": round(average_score, 4)
    }
    
    return {
        "aggregation_result": aggregation_output,
        "identified_regions_count": len(merged_regions)
    } 