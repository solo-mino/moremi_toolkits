"""
Small Molecule Metrics Collector V3

This module implements a comprehensive metrics collection system based on 
SwissADME and ADMETlab metrics as defined in the metrics breakdown document.
This version focuses on collecting metrics without pass/fail judgments for ranking purposes.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Union
import pandas as pd
from enum import Enum
import os

from rdkit import Chem
from rdkit.Chem import Descriptors

# Corrected imports for property calculators
from .property_calculators import (
    ADMETPredictor,
    PhysiochemicalProperties,
    LipophilicityProperties,
    DruglikenessProperties,
    MedicinalChemistryProperties
)

class MetricCategory(Enum):
    """Categories of metrics based on the comprehensive breakdown"""

    PHYSICOCHEMICAL = "Physicochemical"
    MEDICINAL_CHEMISTRY = "Medicinal Chemistry"  
    LIPOPHILICITY = "Lipophilicity"  # Display only
    DRUGLIKENESS = "Druglikeness"
    ABSORPTION = "Absorption"
    DISTRIBUTION = "Distribution"
    METABOLISM = "Metabolism"
    EXCRETION = "Excretion"
    TOXICITY = "Toxicity"

@dataclass
class MetricRanges:
    """
    Reference ranges for molecular properties based on the comprehensive metrics breakdown.
    These ranges are used for normalization and reporting, not for pass/fail criteria.
    """
    
    # Physicochemical Properties
    physicochemical = {
        'tpsa': (0, 140),  # Å² (Ertl et al., 2000)
        'logp': {
            'optimal': (1.0, 3.0),       # Optimal range
            'ranges': [
                (-float('inf'), 0),      # Poor membrane permeability
                (0, 1),                  # Suboptimal
                (1, 3),                  # Optimal
                (3, 5),                  # Good permeability, decreasing solubility
                (5, float('inf'))        # Poor solubility, toxicity concerns
            ]
        }
    }
    
    # Lipophilicity - Display only, does not affect ranking
    lipophilicity = {
        'ilogp': None,
        'xlogp3': None,
        'wlogp': None,
        'mlogp': None,
        'silicos_it': None,
        'consensus_logp': None
    }
    
    # Drug-likeness
    druglikeness = {
        'lipinski': {
            'hbd': 5,                    # H-bond donors
            'hba': 10,                   # H-bond acceptors
            'mw': 500,                   # Molecular weight
            'logp': 5,                   # LogP
            'molar_refractivity': (40, 130)  # Molar refractivity
        },
        'bioavailability': {
            'high': 0.7                  # High bioavailability threshold
        }
    }
    
    # Medicinal Chemistry
    medicinal_chemistry = {
        'synthetic_accessibility': {
            'very_easy': (1, 2),
            'easy': (2, 3),
            'moderate': (3, 5),
            'difficult': (5, 7),
            'very_difficult': (7, 10)
        },
        'qed': {
            'high': 0.67,
            'moderate': (0.35, 0.67),
            'low': 0.35
        }
    }
    
    # Absorption
    absorption = {
        'caco2': {
            'high': -4.7,                # log cm/s
            'moderate': (-5.0, -4.7),
            'low': -5.0
        },
        'pampa': {
            'high': -4.7,                # log cm/s
            'moderate': (-5.0, -4.7),
            'low': -5.0
        },
        'mdck': {
            'high': -4.7,                # log cm/s
            'moderate': (-5.0, -4.7),
            'low': -5.0
        },
        'hia': {
            'high': 0.8,                 # >80%
            'moderate': (0.3, 0.8),
            'low': 0.3
        },
        'pgp_substrate': {
            'substrate': 0.5,
            'non_substrate': 0.3
        },
        'pgp_inhibitor': {
            'inhibitor': 0.5,
            'non_inhibitor': 0.3
        }
    }
    
    # Distribution
    distribution = {
        'vdss': 0.71,                    # L/kg
        'ppb': {
            'high': 90,                  # >90%
            'moderate': (70, 90),
            'low': 70
        },
        'bbb': {
            'high': 0.8,                 # >80%
            'moderate': (0.3, 0.8),
            'low': 0.3
        },
        'fu': {
            'high_binding': 0.1,
            'moderate': (0.1, 0.3),
            'low_binding': 0.3
        }
    }
    
    # Excretion
    excretion = {
        'half_life': {
            'long': 24,                  # hours
            'moderate': (6, 24),
            'short': 6
        },
        'clearance': {
            'low': 3,                    # mL/min/kg
            'intermediate': (3, 8),
            'high': 8
        }
    }
    
    # Metabolism
    metabolism = {
        'cyp_inhibition': {
            'strong': 1,                 # μM
            'moderate': (1, 10),
            'weak': 10
        }
    }
    
    # Toxicity
    toxicity = {
        'herg': {
            'safe': 10,                  # μM
            'intermediate': (1, 10),
            'high_risk': 1
        }
    }
    

@dataclass
class MoleculeMetrics:
    """Container for all calculated metrics for a molecule"""
    smiles: str
    molecular_formula: str
    molecular_weight: float
    
    # Combined Categories
    physicochemical: Dict[str, float]
    medicinal_chemistry: Dict[str, float]
    
    # Other Categories
    lipophilicity: Dict[str, float]
    druglikeness: Dict[str, float]
    absorption: Dict[str, float]
    distribution: Dict[str, float]
    metabolism: Dict[str, float]
    excretion: Dict[str, float]
    toxicity: Dict[str, float]
    
    warnings: List[str]
    
    def to_dict(self) -> Dict:
        """Convert metrics to dictionary format"""
        return {
            'smiles': self.smiles,
            'molecular_formula': self.molecular_formula,
            'molecular_weight': self.molecular_weight,
            'metrics': {
                'physicochemical': self.physicochemical,
                'medicinal_chemistry': self.medicinal_chemistry,
                'lipophilicity': self.lipophilicity,
                'druglikeness': self.druglikeness,
                'absorption': self.absorption,
                'distribution': self.distribution,
                'metabolism': self.metabolism,
                'excretion': self.excretion,
                'toxicity': self.toxicity
            },
            'warnings': self.warnings
        }

@dataclass
class ProcessingResult:
    """Container for molecule processing results"""
    smiles: str
    metrics: Optional[MoleculeMetrics]
    error: Optional[str] = None
    success: bool = True

    def __str__(self) -> str:
        if self.success:
            return f"Success: {self.smiles}"
        return f"Failed: {self.smiles} - Error: {self.error}"

def calculate_molecular_weight(mol):
    """Calculate molecular weight from molecular formula"""
    if not mol:
        return 0
    try:
        return Descriptors.ExactMolWt(mol)
    except:
        return 0

class SmallMoleculeValidator:
    """
    Enhanced small molecule metrics collector that gathers comprehensive metrics
    without making pass/fail judgments.
    """
    
    def __init__(self, ranges: Optional[MetricRanges] = None):
        """Initialize validator with metric ranges"""
        self.ranges = ranges or MetricRanges()
        
    def process_molecule(self, smiles: str) -> ProcessingResult:
        """Process a single molecule and calculate all metrics"""
        try:
            mol = Chem.MolFromSmiles(smiles.strip("."))
            if not mol:
                return ProcessingResult(
                    smiles=smiles,
                    metrics=None,
                    success=False,
                    error="Invalid SMILES"
                )

            
            # Calculate basic properties
            molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            molecular_weight = Descriptors.ExactMolWt(mol)  # Calculate molecular weight
            
            # Convert mol to SMILES for consistent format across all calculators
            
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            # Calculate ADMETlab metrics
            admet_results = ADMETPredictor(canonical_smiles).calculate_admet_properties()
            
            # Calculate Phyiscochemical metrics
            swissadme_tpsa = PhysiochemicalProperties(canonical_smiles).calculate_physiochemical_properties()
            lipophilicity = LipophilicityProperties(canonical_smiles).calculate_lipophilicity()

            physicochemical = {
                'tpsa': swissadme_tpsa.get('TPSA', 0.0),
                'logp': lipophilicity.get('Consensus', 0.0),  # Using consensus logP for physicochemical properties
                'qed': admet_results.get('qed', 0.0)
            }
            
            # Calculate lipophilicity metrics (display only)
            lipophilicity_metrics = {
                'ilogp': lipophilicity.get('iLOGP', 0.0),
                'xlogp3': lipophilicity.get('XLOGP3', 0.0),
                'wlogp': lipophilicity.get('WLOGP', 0.0),
                'mlogp': lipophilicity.get('MLOGP', 0.0),
                'silicos_it': lipophilicity.get('SILICOS-IT', 0.0),
                'consensus_logp': lipophilicity.get('Consensus', 0.0)
            }
            
            # Calculate druglikeness metrics
            druglikeness = DruglikenessProperties(canonical_smiles).calculate_druglikeness()
            druglikeness_metrics = {
                'lipinski_violations': druglikeness['lipinski'].get('violations', 0),
                'bioavailability_score': druglikeness.get('bioavailability_score', 0.0)
            }
            
            # Calculate medicinal chemistry metrics
            medicinal_chem = MedicinalChemistryProperties(canonical_smiles).calculate_medicinal_chemistry()
            medicinal_chemistry = {
                'synthetic_accessibility': medicinal_chem.get('synthetic_accessibility', 0.0),
                'qed': admet_results.get('qed', 0.0)
            }

            # Process ADMET metrics
            absorption = {
                'caco2': admet_results.get('caco2_permeability', 0.0),
                'pampa': admet_results.get('pampa_permeability', 0.0),
                'mdck': admet_results.get('mdck_permeability', 0.0),
                'hia': admet_results.get('hia', 0.0),
                'pgp_substrate': admet_results.get('pgp_substrate', 0.0),
                'pgp_inhibitor': admet_results.get('pgp_inhibitor', 0.0)
            }
            
            distribution = {
                'vdss': admet_results.get('vdss', 0.0),
                'ppb': admet_results.get('ppb', 0.0),
                'bbb': admet_results.get('bbb', 0.0),
                'fu': admet_results.get('fraction_unbound', 0.0)
            }
            
            metabolism = {
                'cyp2c9_inhibition': admet_results.get('cyp2c9_inhibition', 0.0),
                'cyp2d6_inhibition': admet_results.get('cyp2d6_inhibition', 0.0),
                'cyp3a4_inhibition': admet_results.get('cyp3a4_inhibition', 0.0)
            }
            
            excretion = {
                'half_life': admet_results.get('half_life', 0.0),
                'clearance': admet_results.get('plasma_clearance', 0.0)
            }
            
            toxicity = {
                'herg': admet_results.get('herg_inhibition', 0.0)
            }
            
            # Initialize metrics container with all required fields
            metrics = MoleculeMetrics(
                smiles=canonical_smiles,
                molecular_formula=molecular_formula,
                molecular_weight=molecular_weight,
                physicochemical=physicochemical,
                medicinal_chemistry=medicinal_chemistry,
                lipophilicity=lipophilicity_metrics,
                druglikeness=druglikeness_metrics,
                absorption=absorption,
                distribution=distribution,
                metabolism=metabolism,
                excretion=excretion,
                toxicity=toxicity,
                warnings=[]
            )
            
            return ProcessingResult(
                smiles=canonical_smiles,
                metrics=metrics,
                success=True
            )
            
        except Exception as e:
            return ProcessingResult(
                smiles=smiles,
                metrics=None,
                success=False,
                error=str(e)
            )

    def process_molecules(self, input_source: Union[str, List[str]], output_dir: str) -> List[ProcessingResult]:
        """
        Process multiple molecules from a file or a list of SMILES strings.
        
        Args:
            input_source: Either a path to file containing SMILES strings or a list of SMILES strings
            output_dir: Directory to save results
        
        Returns:
            List of ProcessingResult objects containing both successful and failed molecules
        """
        os.makedirs(output_dir, exist_ok=True)
        results = []

        # Handle input based on its type
        if isinstance(input_source, str):
            # Process from file
            with open(input_source, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    # Skip empty lines and comments
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Extract SMILES (assuming first column is SMILES)
                    smiles = line.split()[0]
                    
                    result = self.process_molecule(smiles)
                    results.append(result)
        
        elif isinstance(input_source, list):
            # Process from list of SMILES
            for smiles in input_source:
                if smiles and isinstance(smiles, str):
                    smiles = smiles.strip()
                    result = self.process_molecule(smiles)
                    results.append(result)
        
        else:
            raise TypeError("input_source must be either a file path (str) or a list of SMILES strings")
        
        # Print summary
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        
        print("\nProcessing Summary:")
        print(f"Total molecules: {len(results)}")
        print(f"Successfully processed: {len(successful)}")
        print(f"Failed: {len(failed)}")
        
        if failed:
            print("\nFailed Molecules:")
            for result in failed:
                print(f"SMILES: {result.smiles}")
                print(f"Error: {result.error}")
                print("-" * 80)
        
        return results

    def get_successful_metrics(self, results: List[ProcessingResult]) -> List[MoleculeMetrics]:
        """Extract only successful metrics from processing results"""
        return [r.metrics for r in results if r.success]

    def validate_molecules(self, smiles_list):
        """Validate a list of SMILES strings and return valid molecules with metrics."""
        print(f"\nStarting validation of {len(smiles_list)} molecules...")
        valid_metrics = []
        invalid_count = 0
        
        for idx, smiles in enumerate(smiles_list, 1):
            try:
                print(f"\nValidating molecule {idx}/{len(smiles_list)}")
                print(f"├── SMILES: {smiles}")
                
                # Create molecule from SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"└── Invalid SMILES string")
                    invalid_count += 1
                    continue
                
                # Calculate basic properties
                print(f"├── Calculating basic properties...")
                metrics = MoleculeMetrics(
                    smiles=smiles,
                    molecular_formula=Chem.rdMolDescriptors.CalcMolFormula(mol),
                    molecular_weight=Descriptors.ExactMolWt(mol),
                    mol=mol,
                    physicochemical={},
                    medicinal_chemistry={},
                    lipophilicity={},
                    druglikeness={},
                    absorption={},
                    distribution={},
                    metabolism={},
                    excretion={},
                    toxicity={},
                    warnings=[]
                )
                
                valid_metrics.append(metrics)
                print(f"└── Validation successful")
                print(f"    └── Formula: {metrics.molecular_formula}")
                
            except Exception as e:
                print(f"└── Error during validation: {str(e)}")
                invalid_count += 1
                continue
        
        print(f"\nValidation complete:")
        print(f"├── Total molecules: {len(smiles_list)}")
        print(f"├── Valid molecules: {len(valid_metrics)}")
        print(f"└── Invalid molecules: {invalid_count}\n")
        
        return valid_metrics
    
if __name__ == "__main__":
    # Example usage
    validator = SmallMoleculeValidator()
    results = validator.process_molecules("molecules.smi", "metrics_results") 
