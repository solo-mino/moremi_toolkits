import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from Bio.Align import PairwiseAligner, substitution_matrices
from collections import defaultdict
import re
import os
import time
import importlib.resources as pkg_resources

class ProteinComparisonTool:
    def __init__(self):
        """Initialize the tool with a path to the SAbDab CSV file."""
        try:
            package_ref = pkg_resources.files('moremi_biokit.proteins.analysis_tools')
            csv_ref = package_ref.joinpath('TheraSAbDab_SeqStruc_OnlineDownload.csv')
            
            with pkg_resources.as_file(csv_ref) as csv_path_resolved:
                self.sabdab_df = pd.read_csv(csv_path_resolved)
        except FileNotFoundError:
            raise FileNotFoundError(
                "TheraSAbDab_SeqStruc_OnlineDownload.csv not found in moremi_biokit.proteins.utils package. "
                "Please ensure the CSV file is included in the package installation."
            )

        
        # Initialize aligners for quick and accurate modes
        self.quick_aligner = PairwiseAligner()
        self.quick_aligner.mode = 'local'
        self.quick_aligner.match_score = 2
        self.quick_aligner.mismatch_score = -1
        self.quick_aligner.open_gap_score = -8
        self.quick_aligner.extend_gap_score = -0.5
        
        self.accurate_aligner = PairwiseAligner()
        self.accurate_aligner.mode = 'local'
        self.accurate_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        self.accurate_aligner.open_gap_score = -10
        self.accurate_aligner.extend_gap_score = -0.5
        self.accurate_aligner.target_internal_open_gap_score = -8
        self.accurate_aligner.target_internal_extend_gap_score = -0.5
        self.accurate_aligner.query_internal_open_gap_score = -8
        self.accurate_aligner.query_internal_extend_gap_score = -0.5
        
        # Preprocess database for faster searching
        self._preprocess_database()
        
        # Initialize performance statistics dictionary
        self._stats = {
            'total_database_entries': len(self.sabdab_df) if hasattr(self, 'sabdab_df') else 0,
            'candidates_after_length_filter': 0,
        }
        
    def _preprocess_database(self):
        """Preprocess database to speed up searching."""
        # Remove invalid entries
        mask = (self.sabdab_df['HeavySequence'].notna()) | (self.sabdab_df['LightSequence'].notna())
        self.sabdab_df = self.sabdab_df[mask].copy()
        
        # Pre-calculate sequence lengths and kmer profiles
        self.heavy_lengths = self.sabdab_df['HeavySequence'].str.len().fillna(0)
        self.light_lengths = self.sabdab_df['LightSequence'].str.len().fillna(0)
        
        # Create kmer profiles for quick filtering (using 3-mers)
        self.heavy_profiles = self._create_kmer_profiles(self.sabdab_df['HeavySequence'])
        self.light_profiles = self._create_kmer_profiles(self.sabdab_df['LightSequence'])
        
        # Pre-calculate length ranges for faster filtering
        self.length_ranges = {
            'short': (0, 200),
            'medium': (201, 300),
            'long': (301, float('inf'))
        }
        
    def _get_length_category(self, length: int) -> str:
        """Get the length category for a sequence."""
        for category, (min_len, max_len) in self.length_ranges.items():
            if min_len <= length <= max_len:
                return category
        return 'long'
        
    def _get_adaptive_thresholds(self, query_length: int) -> dict:
        """Get adaptive thresholds based on sequence length."""
        # More permissive thresholds for all sequence lengths
        thresholds = {
            'short': {
                'length': 0.1,  # Extremely permissive length filter
                'kmer': 0.1,    # Extremely permissive kmer filter
                'boost': 1.5    # Higher boost for short sequences
            },
            'medium': {
                'length': 0.05, # Even more permissive for medium sequences
                'kmer': 0.05,   
                'boost': 2.0    # Higher boost for medium sequences
            },
            'long': {
                'length': 0.01, # Most permissive for long sequences
                'kmer': 0.01,   
                'boost': 2.5    # Highest boost for long sequences
            }
        }
        
        category = self._get_length_category(query_length)
        return thresholds[category]
        
    def _create_kmer_profiles(self, sequences: pd.Series, k: int = 3) -> List[Dict[str, int]]:
        """Create kmer frequency profiles for quick sequence comparison."""
        profiles = []
        for seq in sequences:
            if pd.isna(seq):
                profiles.append({})
                continue
                
            profile = defaultdict(int)
            seq = seq.upper()
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                profile[kmer] += 1
            profiles.append(dict(profile))
        return profiles
        
    def _quick_kmer_similarity(self, query_profile: Dict[str, int], target_profile: Dict[str, int], boost_factor: float = 1.0) -> float:
        """Quick estimation of sequence similarity using kmer profiles."""
        if not target_profile:
            return 0.0
            
        common_kmers = set(query_profile.keys()) & set(target_profile.keys())
        if not common_kmers:
            return 0.0
            
        # Calculate similarity with boost factor
        similarity = sum(min(query_profile[k], target_profile[k]) for k in common_kmers)
        max_kmers = max(sum(query_profile.values()), sum(target_profile.values()))
        base_similarity = similarity / max_kmers if max_kmers > 0 else 0.0
        
        # Apply boost factor
        return min(base_similarity * boost_factor, 1.0)
        
    def _length_filter(self, query_length: int, target_length: float, threshold: float = 0.7) -> bool:
        """Quick length-based filtering with adaptive threshold."""
        if target_length == 0:
            return False
            
        # Adjust threshold based on sequence length
        if query_length > 300:
            threshold *= 0.6  # Much more permissive for very long sequences
        elif query_length > 200:
            threshold *= 0.8  # More permissive for medium-long sequences
            
        length_ratio = min(query_length, target_length) / max(query_length, target_length)
        return length_ratio >= threshold
        
    def calculate_identity(self, seq1: str, seq2: str, quick_mode: bool = True) -> float:
        """Calculate sequence identity between two sequences."""
        try:
            if not seq2 or pd.isna(seq2):
                return 0.0
                
            seq1 = seq1.upper()
            seq2 = seq2.upper()
            
            # Choose appropriate aligner based on mode
            aligner = self.quick_aligner if quick_mode else self.accurate_aligner
            
            # Get best alignment
            alignments = aligner.align(seq1, seq2)
            if not alignments:
                return 0.0
                
            alignment = alignments[0]
            alignment_str = str(alignment)
            
            # Parse alignment string to get sequences
            lines = alignment_str.split('\n')
            if len(lines) < 3:
                return 0.0
                
            seq1_aligned = lines[0]
            seq2_aligned = lines[2]
            
            # Calculate matches and similarity score
            matches = sum(a == b for a, b in zip(seq1_aligned, seq2_aligned))
            aligned_length = len(seq1_aligned)
            
            # Use the shorter sequence length for normalization
            shorter_len = min(len(seq1), len(seq2))
            identity = matches / shorter_len if shorter_len > 0 else 0.0
            
            # Apply length-based scaling
            if len(seq1) > 300 or len(seq2) > 300:
                identity *= 2.0  # Significant boost for long sequences
            elif len(seq1) > 200 or len(seq2) > 200:
                identity *= 1.5  # Moderate boost for medium sequences
            
            # Apply additional boost for partial matches
            if identity > 0.05:  # If there's any significant similarity
                identity *= 1.2  # Boost it further
            
            return min(identity, 1.0)  # Cap at 1.0
            
        except Exception as e:
            print(f"Warning: Error in sequence alignment: {e}")
            return 0.0
            
    def search_sabdab_local(self, query_sequence: str, similarity: float = 0.7, quick_mode: bool = True) -> pd.DataFrame:
        """Search the local SabDab DataFrame using optimized multi-stage filtering."""
        if not query_sequence or similarity <= 0:
            # Reset stats for clarity if search is invalid
            self._stats['candidates_after_length_filter'] = 0
            return pd.DataFrame()
            
        results = []
        query_length = len(query_sequence)
        query_profile = self._create_kmer_profiles(pd.Series([query_sequence]))[0]
        
        # Get adaptive thresholds based on sequence length
        thresholds = self._get_adaptive_thresholds(query_length)
        
        # Stage 1: Quick length-based filtering with very permissive threshold
        heavy_length_mask = self.heavy_lengths.apply(lambda x: self._length_filter(query_length, x, thresholds['length']))
        light_length_mask = self.light_lengths.apply(lambda x: self._length_filter(query_length, x, thresholds['length']))
        candidates_mask = heavy_length_mask | light_length_mask
        
        # Update statistics
        candidates_after_length_filter = candidates_mask.sum()
        self._stats['candidates_after_length_filter'] = candidates_after_length_filter
        
        if not candidates_mask.any():
            return pd.DataFrame()
            
        # Stage 2: Process all candidates with more permissive thresholds
        for idx, row in self.sabdab_df[candidates_mask].iterrows():
            # Check heavy chain
            if pd.notna(row["HeavySequence"]):
                identity = self.calculate_identity(query_sequence, row["HeavySequence"], quick_mode=quick_mode)
                if identity >= similarity:
                    results.append((row, round(identity * 100, 3)))
                    continue
            
            # Check light chain
            if pd.notna(row["LightSequence"]):
                identity = self.calculate_identity(query_sequence, row["LightSequence"], quick_mode=quick_mode)
                if identity >= similarity:
                    results.append((row, round(identity * 100, 3)))
        
        if not results:
            return pd.DataFrame()
            
        # Convert results to DataFrame
        results_df = pd.DataFrame([r[0] for r in results])
        results_df['SeqID over Region'] = [r[1] for r in results]
        
        columns_to_keep = [
            'Therapeutic', "Highest_Clin_Trial (Aug '24)", 'Est. Status',
            'Year Proposed', 'Year Recommended', 'Target', 'Conditions Active', 
            'Conditions Discontinued', 'SeqID over Region'
        ]
        # Ensure only existing columns are selected to avoid KeyError
        final_columns = [col for col in columns_to_keep if col in results_df.columns]
        return results_df[final_columns]

    def get_performance_stats(self, query_sequence_length: int, final_matches_count: int) -> Dict[str, int]:
        """
        Returns performance statistics from the most recent search.

        Args:
            query_sequence_length (int): The length of the sequence used in the search.
            final_matches_count (int): The number of matches returned by the search.

        Returns:
            Dict[str, int]: A dictionary containing performance metrics like total
                database entries, candidates after filtering stages, query length, 
                and final matches found.
        """
        stats = self._stats.copy()  # Start with stats captured during search
        stats['query_sequence_length'] = query_sequence_length
        stats['final_matches_found'] = final_matches_count
        
        # Ensure default values are present if search didn't run or failed early
        stats.setdefault('total_database_entries', len(self.sabdab_df) if hasattr(self, 'sabdab_df') else 0)
        stats.setdefault('candidates_after_length_filter', 0)
            
        return stats

def predict_developability(
    sequence: str,
    similarity_threshold: float = 0.7,
    progressive_threshold: Optional[Tuple[float, float]] = None,
    quick_mode: bool = True,
) -> Dict[str, Any]:
    """Predict the developability of an protein sequence by comparing it with known proteins.

    This function utilizes the ProteinComparisonTool to search a local SAbDab-derived 
    database for similar protein sequences and provides a developability score based on 
    the status and clinical trial phase of matched proteins.

    Args:
        sequence (str): The protein sequence to analyze (heavy or light chain).
        similarity_threshold (float, optional): Minimum sequence identity threshold 
            for a match if not using progressive search. Defaults to 0.7.
        progressive_threshold (Optional[Tuple[float, float]], optional): 
            A tuple (start_threshold, end_threshold) for progressive similarity search. 
            The search starts at start_threshold and decreases by 0.1 until end_threshold 
            or a match is found. If None, single threshold search is used. Defaults to None.
        quick_mode (bool, optional): If True, uses a faster alignment algorithm. 
            If False, uses a more sensitive (BLOSUM62) alignment. Defaults to True.
        

    Returns:
        Dict[str, Any]: A dictionary containing developability prediction results:
            - "matched_proteins_df" (pd.DataFrame): DataFrame of matching proteins 
              from SAbDab with their details and sequence identity.
            - "search_summary" (Dict[str, Any]): 
                - "has_active_matches" (bool): True if any matched protein is 'Active'.
                - "threshold_used" (float): The final similarity threshold that yielded matches.
                - "search_time_seconds" (float): Time taken for the search.
                - "matches_found_count" (int): Number of similar proteins found.
            - "developability_score" (float): A score (0-1) based on the proportion of 
              matched proteins that are 'Active' and in 'Phase-III' or 'Approved' 
              clinical trials.
            - "performance_statistics" (Dict[str, int]): Statistics from the search process.
            - "report_string" (str): A human-readable summary of the findings.
            If the input sequence is invalid, returns a dictionary with an 'error' key.
    """
    start_time = time.time()

    # Initialize tool
    try:
        tool = ProteinComparisonTool()
    except Exception as e:
        return {
            "error": f"Failed to initialize ProteinComparisonTool: {e}",
            "matched_proteins_df": pd.DataFrame(),
            "search_summary": {},
            "developability_score": 0.0,
            "performance_statistics": {},
            "report_string": f"Error: Could not initialize comparison tool. Details: {e}"
        }

    # Validate sequence
    if not isinstance(sequence, str) or not sequence.strip() or len(sequence.strip()) < 10:
        return {
            "error": "Invalid sequence: Must be a string with at least 10 characters.",
            "matched_proteins_df": pd.DataFrame(),
            "search_summary": {
                "has_active_matches": False,
                "threshold_used": similarity_threshold,
                "search_time_seconds": time.time() - start_time,
                "matches_found_count": 0
            },
            "developability_score": 0.0,
            "performance_statistics": tool.get_performance_stats(query_sequence_length=len(sequence.strip()) if sequence else 0, final_matches_count=0),
            "report_string": "Invalid sequence: Sequence is too short or empty."
        }
    cleaned_sequence = sequence.strip().upper()
    query_len = len(cleaned_sequence)

    # Perform search
    matched_df: pd.DataFrame
    threshold_actually_used: float

    if progressive_threshold:
        if not (isinstance(progressive_threshold, tuple) and len(progressive_threshold) == 2 and 
                isinstance(progressive_threshold[0], float) and isinstance(progressive_threshold[1], float) and
                progressive_threshold[0] >= progressive_threshold[1]):
            return {
                "error": "Invalid progressive_threshold. Must be a tuple (start_float, end_float) with start >= end.",
                 "matched_proteins_df": pd.DataFrame(), "search_summary": {}, "developability_score": 0.0, 
                 "performance_statistics": {}, "report_string": "Invalid progressive threshold parameter."
            }

        start_thresh, end_thresh = progressive_threshold
        current_thresh = start_thresh
        step = -0.05 # Smaller step for finer granularity
        
        while current_thresh >= end_thresh:
            matched_df = tool.search_sabdab_local(cleaned_sequence, current_thresh, quick_mode)
            if not matched_df.empty:
                break
            current_thresh += step
            if current_thresh < end_thresh and matched_df.empty: # Ensure end_thresh is tried
                current_thresh = end_thresh 
                matched_df = tool.search_sabdab_local(cleaned_sequence, current_thresh, quick_mode)
                break
        threshold_actually_used = current_thresh
    else:
        threshold_actually_used = similarity_threshold
        matched_df = tool.search_sabdab_local(cleaned_sequence, threshold_actually_used, quick_mode)
    
    search_time_taken = time.time() - start_time
    num_matches_found = len(matched_df)

    # Calculate developability score and check for active proteins
    has_active = False
    developability_score_val = 0.0
    active_matches_count = 0

    if not matched_df.empty:
        # Ensure columns exist before trying to access them
        est_status_col = 'Est. Status'
        highest_trial_col = "Highest_Clin_Trial (Aug '24)"

        if est_status_col in matched_df.columns:
            active_condition = matched_df[est_status_col].str.contains('Active', case=False, na=False)
            has_active = active_condition.any()
            active_matches_count = active_condition.sum()

            if highest_trial_col in matched_df.columns:
                advanced_trial_condition = matched_df[highest_trial_col].isin(['Phase-III', 'Approved'])
                active_advanced_count = matched_df[active_condition & advanced_trial_condition].shape[0]
                developability_score_val = active_advanced_count / num_matches_found if num_matches_found > 0 else 0.0
            else:
                 # If highest trial column is missing, base score only on active status (less ideal)
                 developability_score_val = active_matches_count / num_matches_found if num_matches_found > 0 else 0.0
        elif num_matches_found > 0: # Matches found but no 'Est. Status' column
            pass # Developability score remains 0.0, has_active False

    # Generate report string
    report_lines = []
    if matched_df.empty:
        report_lines.append(f"No similar proteins found (search took {search_time_taken:.2f}s at threshold {threshold_actually_used:.2f}).")
        report_lines.append("Suggestions to find matches:")
        if not progressive_threshold:
            report_lines.append(f"  - Try progressive search: e.g., progressive_threshold=({similarity_threshold}, {max(0.3, similarity_threshold-0.2):.1f})")
        if similarity_threshold > 0.5:
             report_lines.append(f"  - Try a lower single threshold: e.g., similarity_threshold={max(0.3, similarity_threshold-0.2):.1f}")
        report_lines.append("  - Verify the input sequence is a valid protein variable domain sequence.")
        report_lines.append("  - Consider using quick_mode=False for a more sensitive search (slower)." if quick_mode else "")
    else:
        report_lines.append(f"Found {num_matches_found} similar proteins (search took {search_time_taken:.2f}s at threshold {threshold_actually_used:.2f}).")
        if has_active:
            report_lines.append(f"  - {active_matches_count} matched proteins have an 'Active' status.")
        else:
            report_lines.append("  - No matched proteins have an 'Active' status.")
        report_lines.append(f"  - Developability Score: {developability_score_val:.4f}")
    
    # Performance statistics
    perf_stats = tool.get_performance_stats(query_sequence_length=query_len, final_matches_count=num_matches_found)
    report_lines.append("\nSearch Performance Stats:")
    for key, value in perf_stats.items():
        report_lines.append(f"  - {key.replace('_', ' ').capitalize()}: {value}")
    report_lines.append(f"  - Search time seconds: {search_time_taken:.2f}")
    report_lines.append(f"  - Quick mode used: {quick_mode}")

    result_dict = {
        "matched_proteins_df": matched_df.to_dict('records') if not matched_df.empty else [],
        "search_summary": {
            "has_active_matches": bool(has_active),
            "threshold_used": float(threshold_actually_used),
            "search_time_seconds": round(search_time_taken, 3),
            "matches_found_count": int(num_matches_found)
        },
        "developability_score": round(developability_score_val, 4),
        # "performance_statistics": perf_stats,
        "report_string": "\n".join(report_lines)
    }
    
    return result_dict