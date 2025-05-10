"""
Batch Antibody Processor for Production Use

This script processes multiple proteins from a sequence file, validates them,
generates comprehensive reports, and ranks them based on various properties.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Optional
import argparse

from .protein_validator import ProteinValidator
from .protein_ranker import ProteinRanker

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('protein_processing.log')
    ]
)

class BatchProteinProcessor:
    def __init__(self, input_file: str, 
                 output_dir: Optional[str] = None, 
                 generate_pdf: bool = True, 
                 generate_csv: bool = True,
                 target_antigen_sequence: Optional[str] = None,
                 target_antigen_pdb_file_path: Optional[str] = None,
                 target_antigen_pdb_chain_id: Optional[str] = None,
                 target_antigen_pdb_id_input: Optional[str] = None,
                 antigen_pdb_download_path: Optional[str] = None
                 ):
        """Initialize the batch processor.
        
        Args:
            input_file: Path to file containing protein sequences (one per line)
            output_dir: Directory to save outputs (default: creates timestamped dir)
            generate_pdf (bool): Whether to generate PDF reports.
            generate_csv (bool): Whether to generate CSV reports.
            target_antigen_sequence (Optional[str]): Directly provided antigen amino acid sequence.
            target_antigen_pdb_file_path (Optional[str]): Path to a local PDB file for the antigen.
            target_antigen_pdb_chain_id (Optional[str]): Antigen identifier as 'PDBID_CHAIN' (e.g., "1XYZ_A").
            target_antigen_pdb_id_input (Optional[str]): PDB ID of the target antigen for PDB file fetching.
            antigen_pdb_download_path (Optional[str]): Directory to store downloaded antigen PDBs.
                                                       If None, a subdir in output_dir will be used.
        """
        self.input_file = Path(input_file)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.generate_pdf = generate_pdf
        self.generate_csv = generate_csv
        
        # Store antigen parameters
        self.target_antigen_sequence = target_antigen_sequence
        self.target_antigen_pdb_file_path = target_antigen_pdb_file_path
        self.target_antigen_pdb_chain_id = target_antigen_pdb_chain_id
        self.target_antigen_pdb_id_input = target_antigen_pdb_id_input
        # self.antigen_pdb_download_path = antigen_pdb_download_path # Store if needed for other purposes, or pass directly

        # Set up output directory
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path("protein_analysis_results") / self.timestamp
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components, passing antigen parameters to ProteinValidator
        validator_pdb_path = self.output_dir / "pdbs" # Define explicitly for clarity
        validator_pdb_path.mkdir(parents=True, exist_ok=True)

        # Determine path for downloaded antigen PDBs
        if antigen_pdb_download_path:
            self.actual_antigen_pdb_download_path = Path(antigen_pdb_download_path)
        else:
            self.actual_antigen_pdb_download_path = self.output_dir / "downloaded_antigen_pdbs"
        self.actual_antigen_pdb_download_path.mkdir(parents=True, exist_ok=True)

        self.validator = ProteinValidator(
            pdb_files_path=str(validator_pdb_path), # For generated PDBs by predict_structure
            target_antigen_sequence=self.target_antigen_sequence,
            target_antigen_pdb_file_path=self.target_antigen_pdb_file_path,
            target_antigen_pdb_chain_id=self.target_antigen_pdb_chain_id,
            target_antigen_pdb_id_input=self.target_antigen_pdb_id_input,
            antigen_pdb_download_path=str(self.actual_antigen_pdb_download_path) # For downloaded antigen PDBs
        )
        self.ranker = ProteinRanker(generate_pdf=self.generate_pdf, generate_csv=self.generate_csv)
        self.ranker.set_output_directory(str(self.output_dir))
        
        # Set up logging file in output directory
        self.log_file = self.output_dir / "processing.log"
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logging.getLogger().addHandler(file_handler)

    def read_sequence_file(self) -> List[dict]:
        """Read protein sequences from input file.
        
        The file should be in FASTA format or a simple text file with one sequence per line.
        For FASTA format, the format should be:
        >protein_name [optional antigen]
        SEQUENCE
        
        For simple text format:
        SEQUENCE
        
        Returns:
            List of dictionaries containing sequence and optional antigen info
        """
        sequences = []
        
        with open(self.input_file, 'r') as f:
            lines = f.readlines()
            
        # Check if file is in FASTA format
        if lines and lines[0].startswith('>'):
            # Process FASTA format
            current_seq = None
            current_name = None
            current_antigen = None
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_seq:
                        sequences.append({
                            'sequence': current_seq,
                            'name': current_name,
                            'antigen': current_antigen
                        })
                    
                    # Parse header
                    header = line[1:].strip()
                    parts = header.split('[')
                    current_name = parts[0].strip()
                    current_antigen = parts[1].rstrip(']').strip() if len(parts) > 1 else None
                    current_seq = ""
                else:
                    # Append to sequence
                    current_seq = (current_seq or "") + line
            
            # Add the last sequence
            if current_seq:
                sequences.append({
                    'sequence': current_seq,
                    'name': current_name,
                    'antigen': current_antigen
                })
        else:
            # Simple one sequence per line format
            for line in lines:
                seq = line.strip()
                if seq:
                    sequences.append({
                        'sequence': seq,
                        'name': None,
                        'antigen': None
                    })
        
        return sequences

    def process_batch(self):
        """Process all proteins in the input file."""
        sequence_data = self.read_sequence_file()
        total_proteins = len(sequence_data)
        
        logging.info(f"🚀 Starting validation of {total_proteins} proteins...")
        print(f"\n🔍 Processing {total_proteins} proteins from {self.input_file.name}")
        
        # Create a temporary file with sequences for the validator
        temp_seq_file = self.output_dir / "sequences.txt"
        with open(temp_seq_file, 'w') as f:
            for data in sequence_data:
                f.write(f"{data['sequence']}\n")
        
        # Process proteins
        results = self.validator.process_proteins(str(temp_seq_file), str(self.output_dir))
        
        successful_metrics = []
        failed_proteins = []
        
        # Track results
        for i, result in enumerate(results):
            seq_data = sequence_data[i]
            seq = seq_data['sequence']
            name = seq_data['name'] or f"Antibody_{i+1}"
            
            if result.success:
                metrics = result.metrics
                # If antigen was provided in the file, use it
                if seq_data['antigen']:
                    metrics.antigen = seq_data['antigen']
                
                successful_metrics.append(metrics)
                formula = metrics.molecular_formula
                logging.info(f"✅ {name}: Validation successful - Formula: {formula}")
                print(f"✅ {name}: Validation successful - Formula: {formula}")
            else:
                error = result.error
                failed_proteins.append((seq, name, error))
                logging.error(f"❌ {name}: Validation failed - {error}")
                print(f"❌ {name}: Validation failed - {error}")
        
        if successful_metrics:
            print(f"\n📊 Ranking {len(successful_metrics)} proteins and generating reports...")
            rankings = self.ranker.rank_proteins(successful_metrics)
            
            # Save failed proteins report if any
            if failed_proteins:
                failed_file = self.output_dir / "failed_proteins.txt"
                with open(failed_file, 'w') as f:
                    f.write("Name\tSequence\tError\n")
                    for seq, name, error in failed_proteins:
                        f.write(f"{name}\t{seq}\t{error}\n")
                print(f"\n⚠️ {len(failed_proteins)} proteins failed processing")
                print(f"   See {failed_file} for details")
            
            print(f"\n✨ Processing complete!")
            print(f"📁 Results saved to: {self.output_dir}")
            print(f"📋 Log file: {self.log_file}")
            
            # Log completion message
            logging.info(f"✨ Processing complete!")
            logging.info(f"📁 Results saved to: {self.output_dir}")
            logging.info(f"📋 Log file: {self.log_file}")
            
            # Return top 5 proteins
            if len(rankings) > 0:
                print("\n🏆 Top 5 proteins:")
                top5 = rankings.head(5)
                for idx, row in top5.iterrows():
                    antigen = row['antigen'] or "Unknown"
                    score = row['total_score']
                    print(f"{idx+1}. Antigen: {antigen} - Score: {score:.4f}")
        else:
            print("\n❌ No valid proteins to process!")
            logging.error("❌ No valid proteins to process!")
            
            # Still show completion message
            print(f"\n✨ Processing complete!")
            print(f"📁 Results saved to: {self.output_dir}")
            print(f"📋 Log file: {self.log_file}")
            
            # Log completion message
            logging.info(f"✨ Processing complete!")
            logging.info(f"📁 Results saved to: {self.output_dir}")
            logging.info(f"📋 Log file: {self.log_file}")

def main():
    parser = argparse.ArgumentParser(description="Batch process protein sequences, validate, rank, and report.")
    parser.add_argument("input_file", help="Path to file containing protein sequences (FASTA or one sequence per line).")
    parser.add_argument("output_dir", nargs='?', default=None, help="Directory to save outputs (default: protein_analysis_results/[timestamp]).")
    parser.add_argument("--pdf", action=argparse.BooleanOptionalAction, default=True, help="Generate PDF reports (--no-pdf to disable).")
    parser.add_argument("--csv", action=argparse.BooleanOptionalAction, default=True, help="Generate CSV reports (--no-csv to disable).")
    
    # Add CLI arguments for antigen parameters
    parser.add_argument("--target-antigen-sequence", type=str, default=None, help="Directly provided antigen amino acid sequence.")
    parser.add_argument("--target-antigen-pdb-file-path", type=str, default=None, help="Path to a local PDB file for the antigen.")
    parser.add_argument("--target-antigen-pdb-chain-id", type=str, default=None, help="Antigen identifier as 'PDBID_CHAIN' (e.g., \"1XYZ_A\").")
    parser.add_argument("--target-antigen-pdb-id-input", type=str, default=None, help="PDB ID of the target antigen for PDB file fetching.")
    parser.add_argument("--antigen-pdb-download-path", type=str, default=None, help="Specific directory to store downloaded antigen PDBs. If not set, a subdirectory 'downloaded_antigen_pdbs' within the main output directory will be used.")

    args = parser.parse_args()
    
    processor = BatchProteinProcessor(
        args.input_file, 
        args.output_dir,
        generate_pdf=args.pdf,
        generate_csv=args.csv,
        target_antigen_sequence=args.target_antigen_sequence,
        target_antigen_pdb_file_path=args.target_antigen_pdb_file_path,
        target_antigen_pdb_chain_id=args.target_antigen_pdb_chain_id,
        target_antigen_pdb_id_input=args.target_antigen_pdb_id_input,
        antigen_pdb_download_path=args.antigen_pdb_download_path
    )
    processor.process_batch()

if __name__ == "__main__":
    main()