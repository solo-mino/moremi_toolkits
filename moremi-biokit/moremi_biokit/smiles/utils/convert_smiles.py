#!/usr/bin/env python3
import os
import logging
from pathlib import Path
from typing import List, Tuple
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import py3Dmol

class MoleculeConverter:
    def __init__(self, smiles_file: str, output_dir: str):
        self.smiles_file = Path(smiles_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.setup_logging()
        
    def setup_logging(self):
        """Configure logging"""
        log_file = self.output_dir / 'conversion.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def read_smiles(self) -> List[str]:
        """Read SMILES strings from input file"""
        try:
            with open(self.smiles_file, 'r') as f:
                return [line.strip() for line in f if line.strip()]
        except Exception as e:
            self.logger.error(f"Error reading SMILES file: {e}")
            raise
            
    def get_molecular_formula(self, mol: rdkit.Chem.rdchem.Mol) -> str:
        """Get molecular formula for the molecule"""
        try:
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            # Clean up the formula to make it suitable for a filename
            formula = formula.replace('(', '').replace(')', '')
            return formula
        except Exception as e:
            self.logger.error(f"Error getting molecular formula: {e}")
            return f"molecule_{hash(Chem.MolToSmiles(mol)) % 10000}"
            
    def convert_smile_to_3d(self, smile: str) -> Tuple[bool, rdkit.Chem.rdchem.Mol]:
        """Convert SMILE to 3D structure with enhanced embedding options"""
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smile)
            if mol is None:
                raise ValueError(f"Failed to parse SMILES: {smile}")
                
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Try different embedding parameters
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useSmallRingTorsions = True
            params.useBasicKnowledge = True
            params.enforceChirality = True
            params.numThreads = 4  # Use multiple threads if available
            
            # Try with more iterations
            params.maxIterations = 2000  # Default is 500
            
            # First attempt with ETKDG v3
            success = AllChem.EmbedMolecule(mol, params)
            
            # If that fails, try with distance geometry and increased attempts
            if success == -1:
                self.logger.warning(f"Initial embedding failed for SMILES: {smile}. Trying alternative approach...")
                
                # Try with UFF pre-optimization
                mp = AllChem.MMFFGetMoleculeProperties(mol)
                if mp is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
                    ff.Initialize()
                    ff.Minimize()
                
                # Try with increased max attempts
                params.maxIterations = 5000
                success = AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)
                
                if success and len(success) > 0:
                    # Use the first successful conformation
                    conf_id = success[0]
                    mol = Chem.Mol(mol)
                    success = 0  # Mark as successful
                else:
                    # Try one more approach with ETDKG basic settings
                    params = AllChem.ETKDG()
                    params.randomSeed = 42
                    success = AllChem.EmbedMolecule(mol, params)
            
            if success == -1:
                raise ValueError(f"Failed to generate 3D conformation for SMILES: {smile} after multiple attempts")
                
            # Optimize the structure with MMFF
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                # If MMFF fails, try UFF
                self.logger.warning("MMFF optimization failed, trying UFF instead")
                AllChem.UFFOptimizeMolecule(mol)
            
            return True, mol
            
        except Exception as e:
            self.logger.error(f"Error converting SMILE to 3D: {e}")
            return False, None
            
    def save_as_pdb(self, mol: rdkit.Chem.rdchem.Mol, output_name: str) -> bool:
        """Save molecule as PDB file"""
        try:
            output_path = self.output_dir / f"{output_name}.pdb"
            writer = Chem.PDBWriter(str(output_path))
            writer.write(mol)
            writer.close()
            self.logger.info(f"Saved PDB file: {output_path}")
            return True
        except Exception as e:
            self.logger.error(f"Error saving PDB file {output_name}: {e}")
            return False
            
    def visualize_molecule(self, mol: rdkit.Chem.rdchem.Mol, output_name: str) -> bool:
        """Create visualization using py3DMol"""
        try:
            # Create py3DMol view
            view = py3Dmol.view(width=800, height=600)
            
            # Convert molecule to PDB format in memory
            pdb = Chem.MolToPDBBlock(mol)
            
            # Add molecule to view
            view.addModel(pdb, "pdb")
            
            # Style the visualization
            view.setStyle({
                'stick': {},
                'sphere': {'radius': 0.5}
            })
            
            # Center and zoom
            view.zoomTo()
            
            # Save to HTML file
            html_path = self.output_dir / f"{output_name}.html"
            view.save(str(html_path))
            
            self.logger.info(f"Saved visualization: {html_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating visualization for {output_name}: {e}")
            return False
            
    def convert_all(self):
        """Convert all SMILES to 3D structures"""
        smiles = self.read_smiles()
        self.logger.info(f"🔍 Starting conversion process for {len(smiles)} molecules")
        
        for i, smile in enumerate(smiles):
            molecule_num = i + 1
            self.logger.info(f"\n└── Molecule {molecule_num}/{len(smiles)}")
            self.logger.info(f"    ├── SMILES: {smile}")
            
            # Convert to 3D
            success, mol = self.convert_smile_to_3d(smile)
            if not success or mol is None:
                self.logger.error("    └── Failed to convert structure")
                continue
            
            # Get molecular formula for naming
            formula = self.get_molecular_formula(mol)
            self.logger.info(f"    ├── Formula: {formula}")
            output_name = formula
            
            if Path(self.output_dir / f"{output_name}.pdb").exists():
                output_name = f"{formula}_{molecule_num}"
                self.logger.info(f"    ├── Using unique name: {output_name}")
            
            # Save as PDB
            if not self.save_as_pdb(mol, output_name):
                self.logger.error("    └── Failed to save PDB")
                continue
            
            # Create visualization
            if self.visualize_molecule(mol, output_name):
                self.logger.info("    └── Visualization complete")
            
        self.logger.info("\n✨ Conversion process completed")

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert SMILES to 3D PDB files with visualizations',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--smiles', required=True,
                      help='Input file containing SMILES strings')
    parser.add_argument('--output', required=True,
                      help='Output directory for PDB files and visualizations')
    
    args = parser.parse_args()
    
    try:
        converter = MoleculeConverter(args.smiles, args.output)
        converter.convert_all()
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
        raise

if __name__ == "__main__":
    main()