"""
Batch Molecule Processor for Production Use

This script processes multiple molecules from a SMILES file, validates them,
generates comprehensive reports, and ranks them based on various properties.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Dict, Any
from tqdm import tqdm
import argparse

from .small_molecule_validator_v3 import SmallMoleculeValidator
from .small_molecule_ranker_v4 import SmallMoleculeRankerV4

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        # TODO: Uncomment the next line to enable file logging
        # logging.FileHandler('molecule_processing.log')
    ]
)

class BatchMoleculeProcessor:
    def __init__(self, input_file: str, output_dir: Optional[str] = None, generate_pdf: Optional[bool] = True, generate_csv: Optional[bool] = True):
        """Initialize the batch processor.
        
        Args:
            input_file: Path to file containing SMILES strings (one per line)
            output_dir: Directory to save outputs (default: creates timestamped dir)
            generate_pdf: Whether to generate PDF reports for each molecule.
            generate_csv: Whether to generate CSV reports for each molecule.
        """
        self.input_file = Path(input_file)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Set up output directory
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path("molecule_analysis_results") / self.timestamp
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.validator = SmallMoleculeValidator()
        self.ranker = SmallMoleculeRankerV4(generate_pdf=generate_pdf, generate_csv=generate_csv)
        self.ranker.set_output_directory(str(self.output_dir))
        
        # Set up logging file in output directory
        self.log_file = self.output_dir / "processing.log"
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logging.getLogger().addHandler(file_handler)

    def read_smiles_file(self) -> List[str]:
        """Read SMILES strings from input file."""
        with open(self.input_file, 'r') as f:
            return [line.strip() for line in f if line.strip()]

    def process_batch(self):
        """Process all molecules in the input file."""
        smiles_list = self.read_smiles_file()
        total_molecules = len(smiles_list)
        
        logging.info(f"🚀 Starting validation of {total_molecules} molecules...")
        print(f"\n🔍 Processing {total_molecules} molecules from {self.input_file.name}")
        
        metrics_list = []
        failed_molecules = []
        
        for i, smiles in enumerate(tqdm(smiles_list, desc="Processing molecules", unit="mol")):
            mol_num = i + 1
            logging.info(f"\n⚡ Processing molecule {mol_num}/{total_molecules}")
            logging.info(f"├── 🧬 SMILES: {smiles}")
            
            try:
                print(f"\n📊 Validating molecule {mol_num}/{total_molecules}")
                print(f"├── 🧬 SMILES: {smiles}")
                print("├── 🔬 Calculating properties...")
                
                validation_result = self.validator.process_molecule(smiles)
                
                if validation_result.success:
                    metrics_list.append(validation_result.metrics)
                    mol_formula = validation_result.metrics.molecular_formula
                    logging.info(f"└── ✅ Validation successful")
                    logging.info(f"    └── 📝 Formula: {mol_formula}")
                    
                    print(f"└── ✅ Validation successful")
                    print(f"    └── 📝 Formula: {mol_formula}")
                else:
                    failed_molecules.append((smiles, validation_result.error))
                    logging.error(f"└── ❌ Validation failed: {validation_result.error}")
                    print(f"└── ❌ Validation failed: {validation_result.error}")
            
            except Exception as e:
                failed_molecules.append((smiles, str(e)))
                logging.error(f"└── ❌ Error processing molecule: {str(e)}")
                print(f"└── ❌ Error processing molecule: {str(e)}")
        
        if metrics_list:
            print("\n📊 Ranking molecules and generating reports...")
            rankings = self.ranker.rank_molecules(metrics_list)
            
            # Save failed molecules report if any
            if failed_molecules:
                failed_file = self.output_dir / "failed_molecules.txt"
                with open(failed_file, 'w') as f:
                    f.write("SMILES\tError\n")
                    for smiles, error in failed_molecules:
                        f.write(f"{smiles}\t{error}\n")
                print(f"\n⚠️ {len(failed_molecules)} molecules failed processing")
                print(f"   See {failed_file} for details")
            
            print(f"\n✨ Processing complete!")
            print(f"📁 Results saved to: {self.output_dir}")
            print(f"📋 Log file: {self.log_file}")
        else:
            print("\n❌ No valid molecules to process!")

    def get_ranked_results_as_dict(self) -> List[Dict[str, Any]]:
        """Returns the final ranked results as a list of dictionaries.

        This method retrieves the ranked data from the internal SmallMoleculeRankerV4 instance.
        It should be called after `process_batch()` has completed.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries representing the ranked molecules.
                                 Returns an empty list if no ranking data is available.
        
        Example:
            >>> # Assuming 'processor' is an instance of BatchMoleculeProcessor
            >>> # and process_batch() has been called.
            >>> results_dict_list = processor.get_ranked_results_as_dict()
            >>> if results_dict_list:
            ...     print(f"Batch processing found {len(results_dict_list)} ranked molecules.")
        """
        if self.ranker and self.ranker.df is not None and not self.ranker.df.empty:
            return self.ranker.get_ranked_data_as_dict()
        return []

def main():
    parser = argparse.ArgumentParser(
        description="Process a batch of SMILES strings, validate, rank, and generate reports."
    )
    parser.add_argument("input_file", help="Path to file containing SMILES strings (one per line)")
    parser.add_argument("-o", "--output_dir", help="Directory to save outputs (default: creates timestamped dir)", default=None)
    
    # Report generation flags
    parser.add_argument("--pdf", action="store_true", help="Generate PDF reports for each molecule. If no report type is specified, both PDF and CSV are generated by default.")
    parser.add_argument("--csv", action="store_true", help="Generate CSV reports for each molecule. If no report type is specified, both PDF and CSV are generated by default.")
    parser.add_argument("--no-pdf", action="store_false", dest="generate_pdf_explicit", help="Explicitly disable PDF report generation.")
    parser.add_argument("--no-csv", action="store_false", dest="generate_csv_explicit", help="Explicitly disable CSV report generation.")

    args = parser.parse_args()

    # Determine report generation based on flags
    # If neither --pdf nor --csv is specified, but --no-pdf or --no-csv IS, respect the --no-* flags.
    # If only specific formats are requested (e.g., --pdf), generate only those.
    # If no format-specific flags are given at all, default to generating both.

    generate_pdf = True
    generate_csv = True

    if hasattr(args, 'generate_pdf_explicit') and not args.generate_pdf_explicit:
        generate_pdf = False
    elif args.pdf: # if --pdf is present and --no-pdf is not
        generate_pdf = True
        if not args.csv and not hasattr(args, 'generate_csv_explicit'): # if only --pdf is specified
            generate_csv = False
            
    if hasattr(args, 'generate_csv_explicit') and not args.generate_csv_explicit:
        generate_csv = False
    elif args.csv: # if --csv is present and --no-csv is not
        generate_csv = True
        if not args.pdf and not hasattr(args, 'generate_pdf_explicit'): # if only --csv is specified
            generate_pdf = False

    processor = BatchMoleculeProcessor(
        args.input_file, 
        args.output_dir,
        generate_pdf=generate_pdf,
        generate_csv=generate_csv
    )
    processor.process_batch()

if __name__ == "__main__":
    main()
