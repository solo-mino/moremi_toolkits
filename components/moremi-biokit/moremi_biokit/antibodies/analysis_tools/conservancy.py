import pandas as pd
from typing import Union, List, Dict, Optional, Any
import os
from datetime import datetime

class EpitopeConservancyAnalyzer:
    def __init__(self, 
                 protein_sequences: Union[str, List[str]], 
                 epitopes: Union[str, List[str]], 
                 identity_threshold: float = 30.0):
        """
        Initialize the analyzer with protein sequences, epitopes, and identity threshold.

        Args:
            protein_sequences: Single protein sequence or list of protein sequences to analyze
            epitopes: Single epitope or list of epitope sequences to search for
            identity_threshold: Minimum sequence identity threshold (threshold >= 30)
        """
        # Normalize inputs to lists if they are strings
        self.protein_sequences = [protein_sequences] if isinstance(protein_sequences, str) else protein_sequences
        self.epitopes = [epitopes] if isinstance(epitopes, str) else epitopes

        # Ensure threshold is above 30
        if identity_threshold < 30:
            raise ValueError("Identity threshold must be above 30")

        self.identity_threshold = identity_threshold

        # Validate inputs
        self._validate_sequences()

    def _validate_sequences(self): 
        """
        Validate input sequences to ensure they are non-empty strings.
        """
        for seq_type, sequences in [
            ("Protein", self.protein_sequences), 
            ("Epitope", self.epitopes)
        ]:
            if not sequences:
                # Allow empty epitope list for the class, but the wrapper function should handle it
                if seq_type == "Epitope": 
                    continue 
                raise ValueError(f"No {seq_type} sequences provided")
            
            for seq in sequences:
                if not isinstance(seq, str):
                    raise TypeError(f"{seq_type} sequences must be strings")
                if not seq.strip():
                    raise ValueError(f"Empty {seq_type} sequence found")

    def calculate_identity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence identity between two sequences.

        Args:
            seq1: First sequence
            seq2: Second sequence

        Returns:
            Float representing percentage identity
        """
        # Ensure sequences are of equal length for pairwise comparison
        len1, len2 = len(seq1), len(seq2)
        if len1 == 0 or len2 == 0 or len1 != len2:
            return 0.0

        # Calculate matches
        matches = sum(c1 == c2 for c1, c2 in zip(seq1.upper(), seq2.upper()))
        return (matches / len1) * 100.0 # Ensure float division

    def find_subsequence_matches(self, epitope: str, protein: str) -> List[tuple]:
        """
        Find all subsequences in the protein that meet the identity threshold.

        Args:
            epitope: Epitope sequence to search for
            protein: Protein sequence to search in

        Returns:
            List of tuples containing (position, subsequence, identity)
        """
        epitope_len = len(epitope)
        protein_len = len(protein)
        matches = []

        if epitope_len == 0 or protein_len < epitope_len:
            return []

        # Iterate through all possible start positions for the subsequence in the protein
        for i in range(protein_len - epitope_len + 1):
            subsequence = protein[i:i + epitope_len]
            identity = self.calculate_identity(epitope, subsequence)

            if identity >= self.identity_threshold:
                matches.append((i, subsequence, identity)) # Store 0-based index

        return matches

    def analyze_epitope(self, epitope: str) -> Dict[str, str]:
        """
        Analyze a single epitope against all protein sequences.

        Args:
            epitope: Epitope sequence to analyze

        Returns:
            Dictionary containing analysis results
        """
        all_match_identities: List[float] = [] # Store all identities found across all proteins
        match_count_above_threshold = 0
        num_proteins = len(self.protein_sequences)

        if num_proteins == 0:
             return { # Handle case where no protein sequences are provided
                'sequence': epitope,
                'percent_matches': "0.00% (0/0)",
                'minimum_identity': "0.00%",
                'maximum_identity': "0.00%",
            }

        for protein in self.protein_sequences:
            matches = self.find_subsequence_matches(epitope.upper(), protein.upper())
            if matches:
                # Find the best identity match for this protein
                best_match_identity = max(m[2] for m in matches)
                all_match_identities.append(best_match_identity)
                if best_match_identity >= self.identity_threshold:
                    match_count_above_threshold += 1
            # If no match in a protein, it implicitly has 0% identity for stats
            # else: 
            #    all_match_identities.append(0.0) # Optional: include 0 if no match found

        percent_matches = (match_count_above_threshold / num_proteins) * 100.0

        # Calculate actual min and max identities from found matches
        min_identity = min(all_match_identities) if all_match_identities else 0.0
        max_identity = max(all_match_identities) if all_match_identities else 0.0

        return {
            'sequence': epitope,
            'percent_matches': f"{percent_matches:.2f}% ({match_count_above_threshold}/{num_proteins})",
            'minimum_identity': f"{min_identity:.2f}%",
            'maximum_identity': f"{max_identity:.2f}%",
        }

    def analyze_all_epitopes(self) -> pd.DataFrame:
        """
        Analyze all epitopes and return results as a DataFrame.

        Returns:
            pandas DataFrame containing analysis results
        """
        if not self.epitopes: # Handle case where no epitopes were provided
             return pd.DataFrame(columns=['Epitope #', 'Epitope Sequence', 
                                        f'Percent matches >= {self.identity_threshold}%', 
                                        'Minimum Identity', 'Maximum Identity'])
                                        
        results = []
        for i, epitope in enumerate(self.epitopes, 1):
            result = self.analyze_epitope(epitope)
            result['Epitope #'] = i # Use the key expected by the wrapper
            # Use keys consistent with the wrapper function's expectations for the final DataFrame
            results.append({
                'Epitope #': result['Epitope #'],
                'Epitope Sequence': result['sequence'],
                f'Percent of protein sequence matches at identity >= {int(self.identity_threshold)}%': result['percent_matches'],
                'Minimum Identity': result['minimum_identity'],
                'Maximum Identity': result['maximum_identity']
            })

        # Define columns explicitly for consistent output
        columns = [
            'Epitope #',
            'Epitope Sequence',
            f'Percent of protein sequence matches at identity >= {int(self.identity_threshold)}%',
            'Minimum Identity',
            'Maximum Identity'
        ]
        df = pd.DataFrame(results, columns=columns)
        return df

    @classmethod
    def from_fasta(cls, fasta_file: str, epitopes: Union[str, List[str]], identity_threshold: float = 30.0):
        """
        Alternative constructor to create analyzer from a FASTA file of protein sequences.

        Args:
            fasta_file: Path to the FASTA file containing protein sequences.
            epitopes: Single epitope or list of epitope sequences to search for.
            identity_threshold: Minimum sequence identity threshold.

        Returns:
            EpitopeConservancyAnalyzer instance.
        """
        try:
            import Bio.SeqIO
        except ImportError:
            raise ImportError("Biopython is required to read FASTA files. Please install it (`pip install biopython`).")

        # Read sequences from FASTA file
        sequences = [str(record.seq) for record in Bio.SeqIO.parse(fasta_file, "fasta")]
        
        if not sequences:
            raise ValueError("No sequences found in the FASTA file")

        return cls(sequences, epitopes, identity_threshold)

# --- Wrapper Function --- 

def predict_conservancy(
    protein_sequences: Union[str, List[str]],
    epitopes: Optional[Union[str, List[str]]],
    identity_threshold: float = 70.0,
    output_dir: Optional[str] = "conservancy_results",
    filename: str = 'epitope_conservancy_analysis',
    save_csv: bool = True,
    display_results: bool = True
) -> Dict[str, Any]:
    """Runs epitope conservancy analysis using EpitopeConservancyAnalyzer.

    This function acts as a convenient wrapper around the EpitopeConservancyAnalyzer class.
    It handles initialization, analysis, result formatting, and optional saving/display.

    Args:
        protein_sequences: Protein sequence(s) (single string or list).
        epitopes: Epitope sequence(s) to search (single string, list, or None).
                  If None or empty list, returns a default result indicating no epitopes.
        identity_threshold: Minimum sequence identity threshold (default: 70.0, must be >= 30).
        output_dir: Directory to save CSV results (default: "conservancy_results"). 
                    Set to None to disable automatic directory creation/saving.
        filename: Base name for the optional output CSV file.
        save_csv: Whether to save the results DataFrame to a CSV file.
        display_results: Whether to print the results DataFrame to the console.

    Returns:
        A dictionary containing:
            - 'results' (pd.DataFrame or str): DataFrame with conservancy analysis results, 
              or a string message if no epitopes were provided.
            - 'conservancy_score' (float): A score from 0 to 1 representing the overall 
              conservation across all analyzed epitopes.
        Returns a dictionary with an 'error' key if initialization or analysis fails.
    """
    # Handle case where no epitopes are provided
    if epitopes is None or (isinstance(epitopes, list) and not epitopes):
        print("Warning: No epitopes provided for conservancy analysis.")
        return {
            "results": "No epitopes provided for analysis.", 
            "conservancy_score": 0.0 
        }

    try:
        # Validate threshold before initializing analyzer
        if not isinstance(identity_threshold, (int, float)) or identity_threshold < 30:
            raise ValueError("Identity threshold must be a number >= 30")

        # Initialize and run analysis using the class
        analyzer = EpitopeConservancyAnalyzer(
            protein_sequences=protein_sequences,
            epitopes=epitopes,
            identity_threshold=identity_threshold
        )
        results_df = analyzer.analyze_all_epitopes()

        # --- Calculate Overall Conservancy Score --- 
        conservancy_score = 0.0
        if not results_df.empty:
            # Extract numeric match percentages (e.g., "75.00% (...) -> 75.00)
            match_col_name = f'Percent of protein sequence matches at identity >= {int(identity_threshold)}%'
            match_percentages = results_df[match_col_name].str.extract(r'(\d+\.\d+)').astype(float).iloc[:, 0] / 100.0
            
            # Get epitope lengths
            epitope_lengths = results_df['Epitope Sequence'].str.len()
            
            # Calculate weighted scores (match % * length)
            weighted_scores = match_percentages * epitope_lengths
            
            # Normalize scores (0 to 1) - relative to the max possible score for this set
            # Max possible score for an epitope = 1.0 (100% match) * length
            max_possible_scores = 1.0 * epitope_lengths
            total_max_possible = max_possible_scores.sum()
            
            if total_max_possible > 0:
                # Overall score is the sum of weighted scores divided by the sum of max possible scores
                conservancy_score = weighted_scores.sum() / total_max_possible
            else: # Handle case with zero-length epitopes or no matches
                conservancy_score = 0.0
        
        # Ensure score is within [0, 1]
        conservancy_score = max(0.0, min(1.0, conservancy_score))

        # --- Optional CSV Saving --- 
        if save_csv and output_dir:
            try:
                os.makedirs(output_dir, exist_ok=True)
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename_timestamp = f"{filename}_threshold_{identity_threshold}_{timestamp}.csv"
                filepath = os.path.join(output_dir, filename_timestamp)
                results_df.to_csv(filepath, index=False)
                print(f"Conservancy results saved to: {filepath}")
            except IOError as e:
                print(f"Warning: Could not save conservancy results to CSV: {e}")

        # --- Optional Console Display --- 
        if display_results:
            print("\nEpitope Conservancy Analysis Results:")
            try:
                # Attempt to print using to_string for better formatting in console
                print(results_df.to_string(index=False))
            except: 
                # Fallback if to_string fails (e.g., very wide columns)
                print(results_df)
            print(f"\nOverall Conservancy Score: {conservancy_score:.4f}")
            
        return {
                "results": results_df,
                "conservancy_score": round(conservancy_score, 4),
            }

    except (ValueError, TypeError) as e:
        # Catch initialization errors (invalid threshold, bad sequences)
        print(f"Error during conservancy analysis setup: {e}")
        return {"error": str(e)}
    except Exception as e:
        # Catch unexpected errors during analysis
        print(f"An unexpected error occurred during conservancy prediction: {e}")
        import traceback
        traceback.print_exc()
        return {"error": f"Unexpected analysis error: {e}"}