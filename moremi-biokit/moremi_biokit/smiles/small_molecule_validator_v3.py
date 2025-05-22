"""
Small Molecule Metrics Collector V3

This module implements a comprehensive metrics collection system based on 
SwissADME and ADMETlab metrics as defined in the metrics breakdown document.
This version focuses on collecting metrics without pass/fail judgments for ranking purposes.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Union, Any
import pandas as pd
from enum import Enum
import os
import logging
import json # Added for JSON operations

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

# Initialize logger for this module
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'MoleculeMetrics':
        """Reconstructs a MoleculeMetrics object from a dictionary (e.g., from JSON)."""
        metrics_data = data.get('metrics', {})
        # Helper to safely get nested metric dicts, defaulting to empty dict if key missing
        def get_metric_category(cat_enum_val: str) -> Dict[str, float]:
            return metrics_data.get(cat_enum_val, metrics_data.get(cat_enum_val.lower().replace(" ", "_"), {}))

        return cls(
            smiles=data.get('smiles', ''),
            molecular_formula=data.get('molecular_formula', 'N/A'),
            molecular_weight=data.get('molecular_weight', 0.0),
            
            physicochemical=get_metric_category(MetricCategory.PHYSICOCHEMICAL.value),
            medicinal_chemistry=get_metric_category(MetricCategory.MEDICINAL_CHEMISTRY.value),
            lipophilicity=get_metric_category(MetricCategory.LIPOPHILICITY.value),
            druglikeness=get_metric_category(MetricCategory.DRUGLIKENESS.value),
            absorption=get_metric_category(MetricCategory.ABSORPTION.value),
            distribution=get_metric_category(MetricCategory.DISTRIBUTION.value),
            metabolism=get_metric_category(MetricCategory.METABOLISM.value),
            excretion=get_metric_category(MetricCategory.EXCRETION.value),
            toxicity=get_metric_category(MetricCategory.TOXICITY.value),
            
            warnings=data.get('warnings', [])
        )

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
    
    def __init__(self):
        """Initialize validator with metric ranges"""
        
        
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

    def _append_molecule_result_to_realtime_csv(self, result: ProcessingResult, csv_path: str):
        """Appends a single molecule ProcessingResult to a CSV file, creating/writing headers if needed."""
        record_to_append = {}
        if result.success and result.metrics:
            metric_data_flat = result.metrics.to_dict() # Get base dict
            record_to_append['smiles'] = metric_data_flat.get('smiles')
            record_to_append['molecular_formula'] = metric_data_flat.get('molecular_formula')
            record_to_append['molecular_weight'] = metric_data_flat.get('molecular_weight')
            record_to_append['success'] = True
            record_to_append['error'] = None
            record_to_append['warnings'] = "; ".join(metric_data_flat.get('warnings', []))

            # Flatten nested metrics from the 'metrics' key
            nested_metrics = metric_data_flat.get('metrics', {})
            for cat_key, cat_value in nested_metrics.items():
                if isinstance(cat_value, dict):
                    for sub_key, sub_value in cat_value.items():
                        record_to_append[f'{cat_key}_{sub_key}'] = sub_value
                else:
                    record_to_append[cat_key] = cat_value # Should not happen often with current structure
        else: # Failed or no metrics
            record_to_append = {
                'smiles': result.smiles,
                'success': False,
                'error': result.error or "No metrics generated",
                'molecular_formula': "N/A",
                'molecular_weight': None,
                'warnings': ""
            }

        # Convert all list/dict values in record to string (though most should be flat now)
        for key, value in record_to_append.items():
            if isinstance(value, (list, dict)):
                record_to_append[key] = str(value)

        df_record = pd.DataFrame([record_to_append])
        file_exists = os.path.isfile(csv_path)
        try:
            df_record.to_csv(csv_path, mode='a', header=not file_exists, index=False)
        except Exception as e_csv:
            logger.error(f"Error appending to molecule realtime CSV {csv_path}: {e_csv}")

    def process_molecules(self, 
                          input_source: Union[str, List[str]], 
                          output_dir: str, # This is for final outputs, not realtime backups
                          realtime_csv_backup_path: Optional[str] = None,
                          realtime_json_backup_path: Optional[str] = None
                          ) -> List[ProcessingResult]:
        """
        Process multiple molecules from a file or a list of SMILES strings.
        Incrementally saves all attempts to realtime_csv_backup_path and successful MoleculeMetrics to realtime_json_backup_path.
        
        Args:
            input_source: Either a path to file containing SMILES strings or a list of SMILES strings
            output_dir: Base directory to save final non-realtime results (e.g., final CSV from save_metrics_to_csv).
            realtime_csv_backup_path (Optional[str]): Path to save all validation attempts incrementally (CSV).
            realtime_json_backup_path (Optional[str]): Path to save successful MoleculeMetrics incrementally (JSON).
        
        Returns:
            List of ProcessingResult objects containing both successful and failed molecules
        """
        # Ensure output_dir for final results exists (though this method doesn't directly save final files itself anymore, it's good practice)
        os.makedirs(output_dir, exist_ok=True)
        
        all_processing_results: List[ProcessingResult] = []
        current_successful_metrics_list: List[MoleculeMetrics] = []
        smiles_list_to_process: List[str] = []

        if isinstance(input_source, str):
            try:
                with open(input_source, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        smiles_list_to_process.append(line.split()[0])
            except FileNotFoundError:
                logger.error(f"Input SMILES file not found: {input_source}")
                return [] # Early exit if file not found
            except Exception as e_file:
                logger.error(f"Error reading SMILES file {input_source}: {e_file}")
                return []
        elif isinstance(input_source, list):
            smiles_list_to_process = [s.strip() for s in input_source if isinstance(s, str) and s.strip()]
        else:
            logger.error("input_source must be either a file path (str) or a list of SMILES strings")
            raise TypeError("input_source must be either a file path (str) or a list of SMILES strings")

        if not smiles_list_to_process:
            logger.warning("No SMILES strings to process.")
            # Touch backup files if paths are provided, so they exist even if empty
            if realtime_csv_backup_path and not os.path.exists(realtime_csv_backup_path):
                 pd.DataFrame([]).to_csv(realtime_csv_backup_path, index=False)
            if realtime_json_backup_path and not os.path.exists(realtime_json_backup_path):
                with open(realtime_json_backup_path, 'w') as f_json_empty:
                    json.dump([], f_json_empty)
            return []

        # Ensure backup directories exist if paths are provided
        if realtime_csv_backup_path:
            # os.path.dirname can return empty string if path is just a filename
            csv_backup_dir = os.path.dirname(realtime_csv_backup_path)
            if csv_backup_dir: os.makedirs(csv_backup_dir, exist_ok=True)
        if realtime_json_backup_path:
            json_backup_dir = os.path.dirname(realtime_json_backup_path)
            if json_backup_dir: os.makedirs(json_backup_dir, exist_ok=True)

        for idx, smiles_str in enumerate(smiles_list_to_process, 1):
            logger.info(f"Processing molecule {idx}/{len(smiles_list_to_process)}: {smiles_str}")
            result = self.process_molecule(smiles_str)
            all_processing_results.append(result)

            if realtime_csv_backup_path:
                try:
                    self._append_molecule_result_to_realtime_csv(result, realtime_csv_backup_path)
                except Exception as e_csv_append:
                    logger.error(f"Failed to append to molecule realtime CSV backup {realtime_csv_backup_path} for {smiles_str}: {e_csv_append}")
            
            if result.success and result.metrics:
                current_successful_metrics_list.append(result.metrics)
                if realtime_json_backup_path:
                    try:
                        dict_list_to_save = [pm.to_dict() for pm in current_successful_metrics_list]
                        with open(realtime_json_backup_path, 'w') as f_json:
                            json.dump(dict_list_to_save, f_json, indent=2)
                    except Exception as e_json_write:
                        logger.error(f"Failed to update molecule realtime JSON backup {realtime_json_backup_path} for {smiles_str}: {e_json_write}")

        successful_count = len(current_successful_metrics_list)
        failed_count = len(all_processing_results) - successful_count
        
        logger.info("\nProcessing Summary:")
        logger.info(f"Total molecules attempted: {len(all_processing_results)}")
        logger.info(f"Successfully processed: {successful_count}")
        logger.info(f"Failed: {failed_count}")
        
        if failed_count > 0:
            logger.warning("\nFailed Molecules Details:")
            for res_fail in all_processing_results:
                if not res_fail.success:
                    logger.warning(f"  SMILES: {res_fail.smiles}, Error: {res_fail.error}")
        
        return all_processing_results

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
    
    def get_metrics_as_dict_list(self, results: List[ProcessingResult]) -> List[Dict]:
        """Converts a list of successful ProcessingResult objects to a list of dictionaries.

        Each dictionary in the list represents the metrics of a successfully processed molecule,
        with nested metric categories flattened for easier conversion to tabular formats.

        Args:
            results (List[ProcessingResult]): A list of ProcessingResult objects,
                typically obtained from `process_molecule` or `process_molecules`.

        Returns:
            List[Dict]: A list of dictionaries, where each dictionary contains
                the flattened metrics for a molecule. Returns an empty list if
                no molecules were processed successfully or if the input list is empty.
        """
        metrics_list = []
        for result in results:
            if result.success and result.metrics:
                metric_dict = result.metrics.to_dict()
                
                # Flatten the 'metrics' dictionary
                flat_dict = {
                    'smiles': metric_dict.get('smiles'),
                    'molecular_formula': metric_dict.get('molecular_formula'),
                    'molecular_weight': metric_dict.get('molecular_weight'),
                }
                
                # Add all individual metrics from nested dictionaries
                for category, category_metrics in metric_dict.get('metrics', {}).items():
                    if isinstance(category_metrics, dict):
                        for key, value in category_metrics.items():
                            flat_dict[f'{category}_{key}'] = value
                    else:
                        # Handle cases where a category might not be a dict (though unlikely based on current structure)
                        flat_dict[category] = category_metrics 
                
                flat_dict['warnings'] = ', '.join(metric_dict.get('warnings', [])) # Join warnings into a string
                
                metrics_list.append(flat_dict)
            elif not result.success:
                logger.warning(f"Skipping failed molecule: {result.smiles} - Error: {result.error}")
        return metrics_list

    def save_metrics_to_csv(self, results: List[ProcessingResult], output_dir: str, filename: str = "molecule_metrics.csv") -> Optional[str]:
        """Saves the metrics of successfully processed molecules to a CSV file.

        The method first converts the processing results to a list of dictionaries
        using `get_metrics_as_dict_list`, then creates a Pandas DataFrame,
        and finally saves it to the specified CSV file.

        Args:
            results (List[ProcessingResult]): A list of ProcessingResult objects.
            output_dir (str): The directory where the CSV file will be saved.
                The directory will be created if it doesn't exist.
            filename (str, optional): The name of the CSV file.
                Defaults to "molecule_metrics.csv".

        Returns:
            Optional[str]: The full path to the saved CSV file if successful,
                otherwise None.
        
        Raises:
            IOError: If there's an issue writing the file.
            Exception: For other potential errors during DataFrame creation or saving.
        """
        if not results:
            logger.warning("No results provided to save_metrics_to_csv. Skipping file generation.")
            return None

        metrics_dict_list = self.get_metrics_as_dict_list(results)

        if not metrics_dict_list:
            logger.info("No successful metrics to save to CSV.")
            return None

        try:
            df = pd.DataFrame(metrics_dict_list)
            
            # Ensure output directory exists
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                logger.info(f"Created output directory: {output_dir}")

            output_path = os.path.join(output_dir, filename)
            
            df.to_csv(output_path, index=False)
            logger.info(f"Successfully saved molecule metrics to {output_path}")
            return output_path

        except IOError as e:
            logger.error(f"IOError saving metrics to CSV: {e}")
            raise
        except Exception as e:
            logger.error(f"An unexpected error occurred while saving metrics to CSV: {e}")
            raise
        return None # Should not be reached if exceptions are raised

if __name__ == "__main__":
    # Example usage
    validator = SmallMoleculeValidator()
    results = validator.process_molecules("molecules.smi", "metrics_results") 
