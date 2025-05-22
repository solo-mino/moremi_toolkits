"""Module for analyzing protein sequences using ProtParam from Biopython."""

from typing import Dict, Tuple, Any, Union

from Bio.SeqUtils.ProtParam import ProteinAnalysis

def analyze_with_protparam(sequence: str) -> Dict[str, Union[Dict[str, Any], str]]:
    """Analyze a protein sequence using the Protparam tool.

    Calculates various physicochemical properties of a protein sequence.

    Args:
        sequence (str): The protein sequence to analyze.

    Returns:
        Dict[str, Union[Dict[str, Any], str]]: A dictionary containing the 
            results of the analysis under the key "protein_params". 
            If an error occurs, returns a dictionary with an "error" key 
            and the error message.
        The "protein_params" dictionary includes:
        - molecular_weight (str): Formatted molecular weight.
        - aromaticity (str): Formatted aromaticity value.
        - instability_index (Tuple[str, str]): Formatted instability index 
            and stability classification (e.g., ("10.0000", "stable")).
        - isoelectric_point (str): Formatted isoelectric point (pI).
        - gravy (str): Formatted GRAVY (Grand Average of Hydropathy) score.
        - hydrophobic_amino_acids (str): Percentage of hydrophobic amino acids.
        - hydrophilic_amino_acids (str): Percentage of hydrophilic amino acids.
        - predicted_solubility (str): "Soluble" or "Insoluble".
        - secondary_structure_fraction (Tuple[float, float, float]): Fraction 
            of helix, turn, and sheet.
    """
    try:
        # Ensure sequence is a string and remove leading/trailing whitespace
        cleaned_sequence = str(sequence).strip()
        if not cleaned_sequence:
            return {"error": "Input sequence cannot be empty."}

        analyzed_seq = ProteinAnalysis(cleaned_sequence)

        molecular_weight = analyzed_seq.molecular_weight()
        aromaticity = analyzed_seq.aromaticity()
        instability_index = analyzed_seq.instability_index()
        isoelectric_point = analyzed_seq.isoelectric_point()
        gravy = analyzed_seq.gravy()
        # Returns a tuple: (Helix, Turn, Sheet)
        secondary_structure_fraction = analyzed_seq.secondary_structure_fraction()
        aa_counts = analyzed_seq.count_amino_acids()
        # aliphatic_index = analyzed_seq.aliphatic_index()

        # Standard classifications for hydrophobic and hydrophilic amino acids
        hydrophobic_aas = ['A', 'I', 'L', 'M', 'F', 'W', 'V', 'P']
        hydrophilic_aas = ['R', 'K', 'D', 'E', 'N', 'Q']

        hydrophobic_count = sum(aa_counts.get(aa, 0) for aa in hydrophobic_aas)
        hydrophilic_count = sum(aa_counts.get(aa, 0) for aa in hydrophilic_aas)
        
        total_aas = len(cleaned_sequence) # Use length of cleaned sequence for total
        if total_aas == 0: # Should be caught by the empty check, but good for robustness
            return {"error": "Cannot analyze zero-length sequence after cleaning."}

        hydrophobic_percentage = (hydrophobic_count / total_aas) * 100 if total_aas > 0 else 0
        hydrophilic_percentage = (hydrophilic_count / total_aas) * 100 if total_aas > 0 else 0

        # GRAVY < 0 often suggests solubility
        solubility = "Soluble" if gravy < 0 else "Insoluble"
        
        stability_classification = "stable" if instability_index < 40 else "unstable"

        protparam_results_dict: Dict[str, Any] = {
            "molecular_weight": f"{molecular_weight:.2f}",
            "aromaticity": f"{aromaticity:.4f}",
            "instability_index": (f"{instability_index:.4f}", stability_classification),
            "isoelectric_point": f"{isoelectric_point:.4f}",
            "gravy": f"{gravy:.4f}",
            "hydrophobic_amino_acids_percentage": f"{hydrophobic_percentage:.2f}%",
            "hydrophilic_amino_acids_percentage": f"{hydrophilic_percentage:.2f}%",
            "predicted_solubility": solubility,
            "secondary_structure_fraction": secondary_structure_fraction # (helix, turn, sheet)
        }
        
        return protparam_results_dict

    except ValueError as ve: # Specific error for invalid amino acids in ProteinAnalysis
        return {"error": f"Invalid character in sequence: {ve}"}
    except Exception as e:
        # Log the exception for debugging: print(f"Unexpected error in analyze_with_protparam: {e}")
        return {"error": f"Error analyzing protein parameters: {e}"} 