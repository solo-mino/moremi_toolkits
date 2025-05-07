"""
Batch Molecule Processor for Production Use

This script processes multiple molecules from a SMILES file, validates them,
generates comprehensive reports, and ranks them based on various properties.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Optional
from tqdm import tqdm

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
    def __init__(self, input_file: str, output_dir: Optional[str] = None):
        """Initialize the batch processor.
        
        Args:
            input_file: Path to file containing SMILES strings (one per line)
            output_dir: Directory to save outputs (default: creates timestamped dir)
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
        self.ranker = SmallMoleculeRankerV4()
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

def main():
    if len(sys.argv) < 2:
        print("Usage: python batch_molecule_processor.py <smiles_file> [output_dir]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    processor = BatchMoleculeProcessor(input_file, output_dir)
    processor.process_batch()

if __name__ == "__main__":
    main()
