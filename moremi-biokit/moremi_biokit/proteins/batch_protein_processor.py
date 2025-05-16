"""
Batch Protein Processor (Refactored for ProteinValidatorV2 and ProteinRankerV2)

This script processes multiple proteins from a sequence file, validates them
using ProteinValidatorV2, and then ranks them using the centralized 
`rank_proteins_from_metrics` function which leverages ProteinRankerV2 logic.
It supports dynamic weighting based on specified metrics and configurable output formats.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Union
import argparse

# Updated imports for V2 components
from .protein_validator_v2 import MetricCategory 
from .protein_ranker import rank_proteins_from_metrics, ScoringConfig # rank_proteins_from_metrics is the main entry point

# Configure logging - This will be the base logger.
# rank_proteins_from_metrics will set up its own run-specific file logger.
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - BatchProcessor - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
        # A general log file for the batch processor itself could be added here if desired,
        # distinct from the run-specific log created by rank_proteins_from_metrics.
    ]
)

class BatchProteinProcessor:
    def __init__(self, 
                 input_file: str, 
                 output_dir_base: Optional[str] = None, 
                 generate_pdf: bool = False, # Default changed
                 generate_csv: bool = True,
                 metrics_to_run_str: Optional[str] = None, # Comma-separated string
                 target_antigen_sequence: Optional[str] = None,
                 target_antigen_pdb_file_path: Optional[str] = None,
                 target_antigen_pdb_chain_id: Optional[str] = None,
                 antigen_pdb_download_path: Optional[str] = None
                 ):
        """Initialize the batch processor.
        
        Args:
            input_file (str): Path to file containing protein sequences (FASTA or one per line).
            output_dir_base (Optional[str]): Base directory to save outputs. 
                                         `rank_proteins_from_metrics` will create a timestamped sub-directory here.
                                         Defaults to "batch_protein_analysis_runs".
            generate_pdf (bool): Whether to generate PDF reports (passed to ranker). Defaults to False.
            generate_csv (bool): Whether to generate CSV reports (passed to ranker). Defaults to True.
            metrics_to_run_str (Optional[str]): Comma-separated string of MetricCategory names to run and score.
                                                If None, ranker uses default weights or all available metrics.
            target_antigen_sequence (Optional[str]): Antigen amino acid sequence.
            target_antigen_pdb_file_path (Optional[str]): Path to a local PDB file for the antigen.
            target_antigen_pdb_chain_id (Optional[str]): Antigen identifier 'PDBID_CHAIN'.
            antigen_pdb_download_path (Optional[str]): Directory for downloaded antigen PDBs.
                                                       If None, ranker uses a subdir in its run output.
        """
        self.input_file = Path(input_file)
        if not self.input_file.is_file():
            logging.error(f"Input file not found: {self.input_file}")
            raise FileNotFoundError(f"Input file not found: {self.input_file}")

        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.generate_pdf = generate_pdf
        self.generate_csv = generate_csv
        
        self.target_antigen_sequence = target_antigen_sequence
        self.target_antigen_pdb_file_path = target_antigen_pdb_file_path
        self.target_antigen_pdb_chain_id = target_antigen_pdb_chain_id
        self.antigen_pdb_download_path = antigen_pdb_download_path

        # Parse metrics_to_run_str into List[MetricCategory]
        self.parsed_metrics_to_run: Optional[List[MetricCategory]] = None
        if metrics_to_run_str:
            self.parsed_metrics_to_run = []
            metric_names = [name.strip() for name in metrics_to_run_str.split(',')]
            for name in metric_names:
                try:
                    # Find the MetricCategory enum member by its value (string name)
                    found_enum = False
                    for enum_member in MetricCategory:
                        if enum_member.value.lower() == name.lower():
                            self.parsed_metrics_to_run.append(enum_member)
                            found_enum = True
                            break
                    if not found_enum:
                        logging.warning(f"MetricCategory '{name}' not recognized. It will be ignored.")
                except ValueError: # Should not happen if MetricCategory values are simple strings
                    logging.warning(f"Invalid MetricCategory name: {name}. It will be ignored.")
            if not self.parsed_metrics_to_run: # If all names were invalid
                self.parsed_metrics_to_run = None 
        
        # Set up base output directory for all runs orchestrated by this batch processor instance
        if output_dir_base:
            self.output_dir_base = Path(output_dir_base)
        else:
            # Default base directory for batch runs if not specified by user
            self.output_dir_base = Path("batch_protein_analysis_runs") / f"batch_job_{self.timestamp}"
        self.output_dir_base.mkdir(parents=True, exist_ok=True)
        
        self.log_file_batch = self.output_dir_base / f"batch_processor_{self.timestamp}.log"
        batch_file_handler = logging.FileHandler(self.log_file_batch)
        batch_file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - BatchProcessor - %(message)s'))
        logging.getLogger().addHandler(batch_file_handler) # Add to root logger for this script

        logging.info(f"BatchProteinProcessor initialized. Base output directory: {self.output_dir_base}")
        if self.parsed_metrics_to_run:
            logging.info(f"Target metrics for ranking: {[m.value for m in self.parsed_metrics_to_run]}")


    def read_sequence_file(self) -> List[dict]:
        """Read protein sequences from input file.
        
        The file can be in FASTA format or a simple text file with one sequence per line.
        For FASTA format, the format should be:
        >protein_name [optional antigen_reference_for_this_protein]
        SEQUENCE
        
        For simple text format:
        SEQUENCE
        
        Returns:
            List of dictionaries containing 'sequence', 'name' (optional), 'antigen_ref' (optional).
        """
        sequences_data: List[dict] = []
        
        try:
            with open(self.input_file, 'r') as f_obj:
                # Read all lines, strip them, and filter out empty/comment lines upfront
                content_lines = [
                    line.strip() for line in f_obj 
                    if line.strip() and not line.strip().startswith("#") and not line.strip().startswith("//")
                ]

            if not content_lines:
                logging.warning(f"No non-empty, non-comment lines found in {self.input_file.name}.")
                return []

            # Determine if the file is likely FASTA by checking for '>'
            is_fasta = any(line.startswith(">") for line in content_lines)
            
            if is_fasta:
                logging.debug(f"FASTA format detected for {self.input_file.name} by BatchProcessor.")
                current_seq_lines: List[str] = []
                current_name: Optional[str] = None
                current_antigen_ref: Optional[str] = None
                
                header_count = 0

                for line_content in content_lines:
                    if line_content.startswith('>'):
                        header_count += 1
                        if current_seq_lines: # Save previous sequence
                            sequences_data.append({
                                'sequence': "".join(current_seq_lines),
                                'name': current_name,
                                'antigen_ref': current_antigen_ref
                            })
                        current_seq_lines = [] # Reset for next sequence
                        
                        header = line_content[1:].strip()
                        parts = header.split('[')
                        current_name = parts[0].strip() if parts[0].strip() else f"UnnamedSeq_Header{header_count}"
                        current_antigen_ref = parts[1].rstrip(']').strip() if len(parts) > 1 else None
                    else: # sequence line
                        if all(c.upper() in "ACDEFGHIKLMNPQRSTVWYX*" for c in line_content):
                            current_seq_lines.append(line_content)
                        else:
                            logging.warning(
                                f"Line in {self.input_file.name} (FASTA body for '{current_name}') "
                                f"contains non-standard characters and was skipped: {line_content[:30]}..."
                            )
                
                if current_seq_lines:
                    sequences_data.append({
                        'sequence': "".join(current_seq_lines),
                        'name': current_name or f"UnnamedSeq_FinalFASTA",
                        'antigen_ref': current_antigen_ref
                    })
            else: # Plain text, one sequence per line
                logging.debug(f"Plain text format detected for {self.input_file.name} by BatchProcessor.")
                for line_num, seq_line in enumerate(content_lines, 1):
                    if all(c.upper() in "ACDEFGHIKLMNPQRSTVWYX*" for c in seq_line):
                        sequences_data.append({
                            'sequence': seq_line,
                            'name': f"SeqFromFile_L{line_num}",
                            'antigen_ref': None
                        })
                    else:
                        logging.warning(
                            f"Line {line_num} in {self.input_file.name} (plain format) "
                            f"contains non-standard characters and was skipped: {seq_line[:30]}..."
                        )
        
        except Exception as e:
            logging.error(f"Error reading sequence file {self.input_file.name} in BatchProcessor: {e}", exc_info=True)
            # sequences_data will be whatever was collected before error, or empty if error was in open()
            # Consider whether to return [] or partial data depending on desired behavior on error
            return [] # Safest to return empty on any error during file read/parse
            
        logging.info(f"BatchProcessor.read_sequence_file read {len(sequences_data)} sequences from {self.input_file.name}")
        return sequences_data

    def process_batch(self):
        """Process all proteins using the rank_proteins_from_metrics function."""
        logging.info(f"🚀 Starting batch processing for {self.input_file.name}...")
        
        sequence_details_list = self.read_sequence_file()
        if not sequence_details_list:
            logging.error("No sequences read from input file. Batch processing aborted.")
            print("❌ No sequences found in input file. Aborting.")
            return

        # The rank_proteins_from_metrics function can take the input file path directly.
        # This is generally preferred as it allows the validator to handle various FASTA formats
        # or simple sequence lists internally.
        # The `sequence_details_list` (from `read_sequence_file`) can be used for enhanced logging here if needed.
        
        # `rank_proteins_from_metrics` will create its own timestamped sub-directory
        # within `self.output_dir_base`.
        
        actual_run_output_dir_name = f"ranking_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        # The path passed to rank_proteins_from_metrics is the PARENT directory for its run-specific folder.
        
        logging.info(f"Calling rank_proteins_from_metrics. Output will be in a subdirectory of: {self.output_dir_base}")
        
        try:
            final_ranked_df = rank_proteins_from_metrics(
                protein_sequences_input=self.input_file, # Pass file path directly
                output_dir=str(self.output_dir_base),    # This is the PARENT dir for the run
                config=None, # Allow user to pass ScoringConfig via CLI if feature is added
                generate_pdf=self.generate_pdf,
                generate_csv=self.generate_csv,
                metrics_to_run=self.parsed_metrics_to_run,
                target_antigen_sequence=self.target_antigen_sequence,
                target_antigen_pdb_file_path=self.target_antigen_pdb_file_path,
                target_antigen_pdb_chain_id=self.target_antigen_pdb_chain_id,
                antigen_pdb_download_path=self.antigen_pdb_download_path
            )

            run_specific_output_dir = self.output_dir_base / actual_run_output_dir_name 
            # Note: rank_proteins_from_metrics internally creates a slightly different timestamped name.
            # For precise path, one might need to list dirs or have rank_proteins_from_metrics return it.
            # For now, we'll just point to the base. The internal log of rank_proteins_from_metrics will have the exact path.

            if final_ranked_df is not None and not final_ranked_df.empty:
                logging.info(f"Successfully processed and ranked {len(final_ranked_df)} proteins.")
                print(f"📊 Successfully processed and ranked {len(final_ranked_df)} proteins.")
                print(f"🏆 Top 5 ranked proteins (from {self.input_file.name}):")
                for idx, row in final_ranked_df.head(5).iterrows():
                    # Attempt to match sequence back to a name from sequence_details_list for richer logging
                    protein_name = f"Protein (Seq: {row['sequence'][:20]}...)" # Default if no name match
                    for detail in sequence_details_list:
                        if detail['sequence'] == row['sequence']:
                            protein_name = detail['name'] or protein_name
                            break
                    score = row['total_score']
                    print(f"  {idx+1}. Name: {protein_name} - Score: {score:.4f}")
            elif final_ranked_df is not None: # Empty dataframe
                logging.info("Processing completed, but no proteins were successfully ranked (e.g., all failed validation or no valid metrics).")
                print("ℹ️ Processing completed, but no proteins were ranked. Check logs for details.")
            else: # final_ranked_df is None, indicating a more significant issue in rank_proteins_from_metrics
                logging.error("Processing failed. `rank_proteins_from_metrics` returned None.")
                print("❌ Processing failed critically. Check logs for details.")

        except Exception as e:
            logging.error(f"An error occurred during batch processing: {e}", exc_info=True)
            print(f"❌ An critical error occurred: {e}. Check logs in {self.output_dir_base}")

        finally:
            logging.info(f"✨ Batch processing finished for {self.input_file.name}.")
            logging.info(f"📁 Results and detailed logs are in a sub-directory within: {self.output_dir_base}")
            logging.info(f"📋 Batch processor log: {self.log_file_batch}")
            print(f"✨ Batch processing complete!")
            print(f"📁 Results and detailed logs from the run are in a subdirectory within: {self.output_dir_base}")
            print(f"📋 Batch processor log: {self.log_file_batch}")


def main():
    parser = argparse.ArgumentParser(
        description="Batch process protein sequences using ProteinValidatorV2 and ProteinRankerV2 logic.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_file", 
        help="Path to file containing protein sequences (FASTA or one sequence per line)."
    )
    parser.add_argument(
        "output_dir_base", 
        nargs='?', 
        default=None, 
        help="Base directory to save outputs. A run-specific timestamped subdirectory will be created here. Default: batch_protein_analysis_runs/batch_job_[timestamp]"
    )
    parser.add_argument(
        "--pdf", 
        action=argparse.BooleanOptionalAction, 
        default=False, # Default changed
        help="Generate PDF reports (--no-pdf to disable)."
    )
    parser.add_argument(
        "--csv", 
        action=argparse.BooleanOptionalAction, 
        default=True, 
        help="Generate CSV reports (--no-csv to disable)."
    )
    parser.add_argument(
        "--metrics-to-run", 
        type=str, 
        default=None, 
        help="Comma-separated string of metric category names to run and use for ranking (e.g., 'Structure,ProtParam,Binding Affinity'). If provided, these metrics get equal weight summing to 1.0. Otherwise, default ranker weights are used."
    )
    
    # Antigen parameters (passed to ProteinValidatorV2 via rank_proteins_from_metrics)
    antigen_group = parser.add_argument_group('Antigen Parameters (for metrics like Binding Affinity)')
    antigen_group.add_argument("--target-antigen-sequence", type=str, default=None, help="Directly provided antigen amino acid sequence.")
    antigen_group.add_argument("--target-antigen-pdb-file-path", type=str, default=None, help="Path to a local PDB file for the antigen.")
    antigen_group.add_argument("--target-antigen-pdb-chain-id", type=str, default=None, help="Antigen identifier as 'PDBID_CHAIN' (e.g., \"1XYZ_A\").")
    antigen_group.add_argument("--antigen-pdb-download-path", type=str, default=None, help="Specific directory to store downloaded/generated antigen PDBs. If not set, ranker uses a subdirectory in its run output.")

    args = parser.parse_args()
    
    try:
        processor = BatchProteinProcessor(
            input_file=args.input_file, 
            output_dir_base=args.output_dir_base,
            generate_pdf=args.pdf,
            generate_csv=args.csv,
            metrics_to_run_str=args.metrics_to_run,
            target_antigen_sequence=args.target_antigen_sequence,
            target_antigen_pdb_file_path=args.target_antigen_pdb_file_path,
            target_antigen_pdb_chain_id=args.target_antigen_pdb_chain_id,
            antigen_pdb_download_path=args.antigen_pdb_download_path
        )
        processor.process_batch()
    except FileNotFoundError:
        # Handled by __init__ logging, just exit gracefully for CLI
        sys.exit(1)
    except Exception as e:
        logging.getLogger().critical(f"Unhandled exception in BatchProteinProcessor execution: {e}", exc_info=True)
        print(f"A critical error occurred. Check logs. Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()