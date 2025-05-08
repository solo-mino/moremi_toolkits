"""
Antibody Metrics Collector V3

This module implements a comprehensive metrics collection system for antibodies based on 
the metrics breakdown document. This version focuses on collecting all metrics for validation
and reporting purposes, while also calculating weighted scores for ranking.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Any
import pandas as pd
from enum import Enum
import os
import logging
from datetime import datetime
from .utils import pdb_fetcher
import random

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

from .analysis_tools import (
    perform_blast, analyze_with_protparam, predict_stability, predict_aggregation,
    predict_glycosylation, predict_structure, predict_bcell_epitopes,
    predict_developability, predict_conservancy, predict_binding_affinity, predict_immunogenicity
)



class MetricCategory(Enum):
    """Categories of metrics based on the comprehensive breakdown"""
    BLAST = "BLAST"
    PROTPARAM = "ProtParam"
    IMMUNOGENICITY = "Immunogenicity"
    STABILITY = "Stability"
    AGGREGATION = "Aggregation"
    GLYCOSYLATION = "Glycosylation"
    BINDING_AFFINITY = "Binding Affinity"
    STRUCTURE = "Structure"
    EPITOPE = "Epitope"
    DEVELOPABILITY = "Developability"
    CONSERVANCY = "Conservancy"
    

@dataclass
class MetricRanges:
    """
    Reference ranges for antibody properties based on the comprehensive metrics breakdown.
    These ranges are used for normalization and reporting, not for pass/fail criteria.
    """
    
    # Binding Affinity (Weight: 0.30)
    binding_affinity = {
        'kd_range': (1e-3, 1e-12),  # Dissociation constant range
        'optimal': 1e-12  # Optimal Kd value
    }
    
    # GMQE Score (Weight: 0.20)
    gmqe = {
        'range': (0.6, 1.0),
        'optimal': 1.0
    }
    
    # Glycosylation Sites (Weight: 0.10)
    glycosylation = {
        'n_sites_range': (1, 4),
        'optimal': 4
    }
    
    # Aggregation (Weight: 0.10)
    aggregation = {
        'propensity': ['Low', 'Medium', 'High'],
        'regions_range': (0, 8),
        'optimal_propensity': 'Low',
        'optimal_regions': 0
    }
    
    # ProtParam (Weight: 0.10)
    protparam = {
        'gravy_range': (-1.5, -0.5),
        'optimal_gravy': -0.5,
        'solubility': ['Soluble', 'Insoluble'],
        'optimal_solubility': 'Soluble'
    }
    
    # Immunogenicity (Weight: 0.05)
    immunogenicity = {
        'score_range': (0, 1),
        'optimal': 1.0
    }
    
    # Conservancy (Weight: 0.05)
    conservancy = {
        'score_range': (0, 1),
        'optimal': 1.0
    }
    
    # Stability (Weight: 0.05)
    stability = {
        'tm_range': (65, 90),  # Melting temperature range in °C
        'optimal': 90
    }
    
    # Epitope (Weight: 0.03)
    epitope = {
        'score_range': (0, 1),
        'optimal': 1.0
    }
    
    # Developability (Weight: 0.02)
    developability = {
        'score_range': (0, 1),
        'optimal': 1.0
    }

@dataclass
class AntibodyMetrics:
    """Container for all calculated metrics for an antibody"""
    sequence: str
    antigen: str
    molecular_weight: float
    molecular_formula: str
    
    # Combined Categories
    blast: Dict[str, Any]
    protparam: Dict[str, Any]
    immunogenicity: Dict[str, Any]
    stability: Dict[str, Any]
    aggregation: Dict[str, Any]
    glycosylation: Dict[str, Any]
    binding_affinity: Dict[str, Any]
    structure: Dict[str, Any]
    epitope: Dict[str, Any]
    developability: Dict[str, Any]
    conservancy: Dict[str, Any]
    
    # Calculated scores
    weighted_scores: Dict[str, float]
    total_score: float
    
    warnings: List[str]
    
    def to_dict(self) -> Dict:
        """Convert metrics to dictionary format"""
        return {
            'sequence': self.sequence,
            'antigen':self.antigen,
            'molecular_weight': self.molecular_weight,
            'molecular_formula': self.molecular_formula,
            'metrics': {
                'blast': self.blast,
                'protparam': self.protparam,
                'immunogenicity': self.immunogenicity,
                'stability': self.stability,
                'aggregation': self.aggregation,
                'glycosylation': self.glycosylation,
                'binding_affinity': self.binding_affinity,
                'structure': self.structure,
                'epitope': self.epitope,
                'developability': self.developability,
                'conservancy': self.conservancy
            },
            'weighted_scores': self.weighted_scores,
            'total_score': self.total_score,
            'warnings': self.warnings
        }

@dataclass
class ProcessingResult:
    """Container for antibody processing results"""
    sequence: str
    metrics: Optional[AntibodyMetrics]
    error: Optional[str] = None
    success: bool = True

    def __str__(self) -> str:
        if self.success:
            return f"Success: {self.sequence[:20]}..."
        return f"Failed: {self.sequence[:20]}... - Error: {self.error}"

class AntibodyValidator:
    """
    Enhanced antibody metrics collector that gathers comprehensive metrics
    and calculates weighted scores for ranking.
    """
    
    def __init__(self, ranges: Optional[MetricRanges] = None, pdb_files_path: Optional[str] = None, target_antigen_pdb_id: Optional[str] = None):
        """Initialize validator with metric ranges and optional target antigen.
        
        Args:
            ranges (Optional[MetricRanges]): Metric ranges for scoring.
            pdb_files_path (Optional[str]): Path to store/find PDB files.
            target_antigen_pdb_id (Optional[str]): Specific PDB ID of the target antigen. 
                                                   If None, a random internal antigen is chosen.
        """
        self.ranges = ranges or MetricRanges()
        self.pdb_files_path = pdb_files_path
        
        self.target_antigen_pdb_id: Optional[str] = None
        self.target_antigen_pdb_path: Optional[str] = None

        if target_antigen_pdb_id:
            self.target_antigen_pdb_id = target_antigen_pdb_id
        else:
            internal_pdb_ids = pdb_fetcher.list_internal_pdb_ids()
            if internal_pdb_ids:
                self.target_antigen_pdb_id = random.choice(internal_pdb_ids)
                logging.info(f"No target_antigen_pdb_id provided, randomly selected: {self.target_antigen_pdb_id}")
            else:
                logging.warning("No internal PDBs found for default target antigen selection.")

        if self.target_antigen_pdb_id:
            if self.pdb_files_path:
                try:
                    os.makedirs(self.pdb_files_path, exist_ok=True)
                    antigen_content = pdb_fetcher.read_internal_pdb(self.target_antigen_pdb_id)
                    if antigen_content:
                        destination_path = os.path.join(self.pdb_files_path, f"{self.target_antigen_pdb_id}.pdb")
                        with open(destination_path, 'w') as f_out:
                            f_out.write(antigen_content)
                        self.target_antigen_pdb_path = destination_path
                        logging.info(f"Using internal antigen '{self.target_antigen_pdb_id}', materialised at '{self.target_antigen_pdb_path}'")
                    else:
                        logging.warning(f"Could not read content for internal antigen PDB '{self.target_antigen_pdb_id}'. Binding affinity prediction might be affected.")
                except Exception as e:
                    logging.error(f"Error preparing antigen PDB '{self.target_antigen_pdb_id}' in '{self.pdb_files_path}': {e}")
            else:
                logging.warning("pdb_files_path not set in AntibodyValidator; cannot prepare target antigen PDB path locally. Will rely on direct path from pdb_fetcher if available for binding.")
                # Attempt to get a direct path if pdb_files_path is not for writing
                fetched_info = pdb_fetcher.fetch_pdb(identifier=self.target_antigen_pdb_id, use_internal=True) # output_dir not specified
                if fetched_info and fetched_info.get("status") == "success" and fetched_info.get("pdb_path"):
                    self.target_antigen_pdb_path = fetched_info.get("pdb_path")
                    logging.info(f"Using direct path for internal antigen '{self.target_antigen_pdb_id}': {self.target_antigen_pdb_path}")
                else:
                    logging.warning(f"Could not obtain a direct path for internal antigen PDB '{self.target_antigen_pdb_id}'.")

        elif not target_antigen_pdb_id: # Only log if it wasn't explicitly set to None/empty
            logging.info("No target_antigen_pdb_id selected or available.")
        
    def calculate_weighted_scores(self, metrics: Dict[str, Any]) -> Tuple[Dict[str, float], float]:
        """Calculate weighted scores for each metric category"""
        weighted_scores = {}
        
        # Binding Affinity (0.30)
        if 'binding_affinity' in metrics:
            kd = float(metrics['binding_affinity'].get('dissociation_constant', '1e-3').replace('e-', ''))
            if 3 <= kd <= 12:  # Within valid range
                weighted_scores['binding_affinity'] = (12 - kd) / 9 * 0.30
            else:
                weighted_scores['binding_affinity'] = 0.0
        
        # GMQE Score (0.20)
        if 'structure' in metrics and 'gmqe_score' in metrics['structure']:
            gmqe = metrics['structure']['gmqe_score']
            if 0.6 <= gmqe <= 1.0:
                weighted_scores['gmqe'] = (gmqe - 0.6) / 0.4 * 0.20
            else:
                weighted_scores['gmqe'] = 0.0
        
        # Glycosylation Sites (0.10)
        if 'glycosylation' in metrics:
            n_sites = metrics['glycosylation'].get('n_glyc_sites_count', 0)
            if 1 <= n_sites <= 4:
                weighted_scores['glycosylation'] = (n_sites - 1) / 3 * 0.10
            else:
                weighted_scores['glycosylation'] = 0.0
        
        # Aggregation (0.10)
        if 'aggregation' in metrics:
            agg_result = metrics['aggregation']
            propensity = agg_result.get('aggregation_propensity', 'High')
            regions = len(agg_result.get('aggregation_prone_regions', []))
            
            # Propensity score (0.05)
            if propensity == 'Low':
                weighted_scores['aggregation_propensity'] = 0.05
            elif propensity == 'Medium':
                weighted_scores['aggregation_propensity'] = 0.025
            else:
                weighted_scores['aggregation_propensity'] = 0.0
            
            # Regions score (0.05)
            if 0 <= regions <= 8:
                weighted_scores['aggregation_regions'] = (8 - regions) / 8 * 0.05
            else:
                weighted_scores['aggregation_regions'] = 0.0
        
        # ProtParam (0.10)
        if 'protparam' in metrics:
            protparam = metrics['protparam']
            gravy = float(protparam.get('gravy', 0))
            solubility = protparam.get('predicted_solubility', 'Insoluble')
            
            # GRAVY score (0.05)
            if -1.5 <= gravy <= -0.5:
                weighted_scores['gravy'] = (gravy + 1.5) / 1.0 * 0.05
            else:
                weighted_scores['gravy'] = 0.0
            
            # Solubility score (0.05)
            weighted_scores['solubility'] = 0.05 if solubility == 'Soluble' else 0.0
        
        # Immunogenicity (0.05)
        if 'immunogenicity' in metrics:
            score = metrics['immunogenicity'].get('score', 0)
            if 0 <= score <= 1:
                weighted_scores['immunogenicity'] = score * 0.05
            else:
                weighted_scores['immunogenicity'] = 0.0
        
        # Conservancy (0.05)
        if 'conservancy' in metrics:
            score = metrics['conservancy'].get('conservancy_score', 0)
            if 0 <= score <= 1:
                weighted_scores['conservancy'] = score * 0.05
            else:
                weighted_scores['conservancy'] = 0.0
        
        # Stability (0.05)
        if 'stability' in metrics:
            tm = metrics['stability'].get('melting_temperature', 0)
            if 65 <= tm <= 90:
                weighted_scores['stability'] = (tm - 65) / 25 * 0.05
            else:
                weighted_scores['stability'] = 0.0
        
        # Epitope (0.03)
        if 'epitope' in metrics:
            score = metrics['epitope'].get('score', 0)
            if 0 <= score <= 1:
                weighted_scores['epitope'] = score * 0.03
            else:
                weighted_scores['epitope'] = 0.0
        
        # Developability (0.02)
        if 'developability' in metrics:
            score = metrics['developability'].get('developability_score', 0)
            if 0 <= score <= 1:
                weighted_scores['developability'] = score * 0.02
            else:
                weighted_scores['developability'] = 0.0
        
        total_score = sum(weighted_scores.values())
        return weighted_scores, total_score
    
    def process_antibody(self, sequence: str) -> ProcessingResult:
        """Process a single antibody sequence and calculate all metrics"""
        try:
            # Basic sequence validation
            if not sequence or len(sequence) < 10:
                logging.warning(f"Invalid sequence: sequence is too short or empty")
                return ProcessingResult(
                    sequence=sequence,
                    metrics=None,
                    success=False,
                    error="Invalid sequence: sequence is too short or empty"
                )
            
            # Calculate basic properties
            analyzed_seq = ProteinAnalysis(sequence)
            molecular_weight = analyzed_seq.molecular_weight()
            molecular_formula = f"C{len(sequence)}H{len(sequence)*2}N{len(sequence)}O{len(sequence)}"
            
            # Collect all metrics
            metrics = {}
            
            # Print tree-style validation process header
            logging.info(f"\n🧪 Processing antibody sequence: {sequence[:20]}...")
            logging.info(f"├── Calculating basic properties...")
            print(f"\n🧪 Processing antibody sequence: {sequence[:20]}...")
            print(f"├── Calculating basic properties...")
            
            # BLAST analysis
            logging.info(f"├── 🧬 Running BLAST analysis....")
            print("├── 🧬 Running BLAST analysis...")
            blast_result = perform_blast(sequence)
            metrics['blast'] = blast_result.get('blast_result', {})
            logging.info(f"│   └── ✓ BLAST complet")
            print("│   └── ✓ BLAST complete")
            
            # ProtParam analysis
            logging.info(f"├── 🔍 Calculating ProtParam properties...")
            print("├── 🔍 Analyzing ProtParam properties...")
            protparam_result = analyze_with_protparam(sequence)
            metrics['protparam'] = protparam_result.get('protein_params', {})
            logging.info(f"│   └── ✓ ProtParam complete")
            print("│   └── ✓ ProtParam complete")
            
            # Immunogenicity prediction
            logging.info(f"├── 🦠 Assessing immunogenicity...")
            print("├── 🦠 Assessing immunogenicity...")
            immunogenicity_result = predict_immunogenicity(sequence)
            metrics['immunogenicity'] = immunogenicity_result
            logging.info(f"│   └── ✓ Immunogenicity assessment complete")
            print("│   └── ✓ Immunogenicity assessment complete")
            
            # Stability prediction
            logging.info(f"├── 🔥 Evaluating stability...")
            print("├── 🔥 Evaluating stability...")
            stability_result = predict_stability(sequence)
            metrics['stability'] = stability_result.get('stability_result', {})
            logging.info(f"│   └── ✓ Stability evaluation complete")
            print("│   └── ✓ Stability evaluation complete")
            
            # Aggregation prediction
            logging.info(f"├── 🧱 Predicting aggregation propensity...")
            print("├── 🧱 Predicting aggregation propensity...")
            aggregation_result = predict_aggregation(sequence)
            metrics['aggregation'] = aggregation_result.get('aggregation_result', {})
            logging.info(f"│   └── ✓ Aggregation prediction complete")
            print("│   └── ✓ Aggregation prediction complete")
            
            # Glycosylation prediction
            logging.info(f"├── 🍭  Identifying glycosylation sites...")
            print("├── 🍭 Identifying glycosylation sites...")
            glycosylation_result = predict_glycosylation(sequence)
            metrics['glycosylation'] = glycosylation_result
            logging.info(f"│   └── ✓ Glycosylation sites identified")
            print("│   └── ✓ Glycosylation sites identified")
            
            # Structure prediction
            logging.info(f"├── 🧩 Generating structural model...")
            print("├── 🧩 Generating structural model...")
            structure_result = predict_structure(sequence, filename=molecular_formula, directory=self.pdb_files_path)
            metrics['structure'] = structure_result
            logging.info(f"│   └── ✓ Structural model generated.")
            print("│   └── ✓ Structural model generated")
            
            # Binding affinity prediction
            logging.info(f"├── 🔗 Calculating binding affinity...")
            print("├── 🔗 Calculating binding affinity...")
            if self.target_antigen_pdb_path and structure_result.get("pdb_file_path"):
                binding_result = predict_binding_affinity( self.target_antigen_pdb_path, structure_result["pdb_file_path"])
                metrics['binding_affinity'] = binding_result
                logging.info(f"│   └── ✓ Binding affinity calculated against {self.target_antigen_pdb_id}")
                print(f"│   └── ✓ Binding affinity calculated against {self.target_antigen_pdb_id}")
            else:
                metrics['binding_affinity'] = {"error": "Receptor(Antigen) PDB or antibody PDB path not available"}
                logging.warning(f"│   └── ⚠️ Skipping binding affinity: Ligand(Antibody) PDB path: {structure_result.get('pdb_file_path')}, Receptor(Antigen) PDB path: {self.target_antigen_pdb_path}")
                print(f"│   └── ⚠️ Skipping binding affinity: Receptor(Antigen) PDB path ({self.target_antigen_pdb_id}) not prepared or Ligand(Antibody) PDB not available.")
            
            # Epitope prediction
            logging.info(f"├── 🎯 Predicting epitopes...")
            print("├── 🎯 Predicting epitope regions...")
            epitope_result = predict_bcell_epitopes(sequence)
            metrics['epitope'] = epitope_result
            logging.info(f"│   └── ✓ Epitope prediction complete")
            print("│   └── ✓ Epitope prediction complete")
            
            # Conservancy prediction
            logging.info(f"├── 🌐 Analyzing sequence conservancy...")
            print("├── 🌐 Analyzing sequence conservancy...")
            conservancy_result = predict_conservancy(sequence, epitope_result.get('list_of_epitopes', []), save_csv = False, display_results = False)
            metrics['conservancy'] = conservancy_result
            logging.info(f"│   └── ✓ Conservancy analysis complete")
            print("│   └── ✓ Conservancy analysis complete")
            
            # Developability prediction
            logging.info(f"├── 🔧 Assessing developability...")
            print("├── 🔧 Assessing developability...")
            developability_result = predict_developability(sequence)
            metrics['developability'] = developability_result
            print("│   └── ✓ Developability assessment complete")
            
            # Calculate weighted scores
            logging.info(f"└── 📊 Calculating weighted scores...")
            print("└── 📊 Calculating final scores")
            weighted_scores, total_score = self.calculate_weighted_scores(metrics)
            
            # Initialize metrics container
            antibody_metrics = AntibodyMetrics(
                sequence=sequence,
                antigen=self.target_antigen_pdb_id or "Unknown", # Use PDB ID as antigen identifier
                molecular_weight=molecular_weight,
                molecular_formula=molecular_formula,
                blast=metrics['blast'],
                protparam=metrics['protparam'],
                immunogenicity=metrics['immunogenicity'],
                stability=metrics['stability'],
                aggregation=metrics['aggregation'],
                glycosylation=metrics['glycosylation'],
                binding_affinity=metrics['binding_affinity'],
                structure=metrics['structure'],
                epitope=metrics['epitope'],
                developability=metrics['developability'],
                conservancy=metrics['conservancy'],
                weighted_scores=weighted_scores,
                total_score=total_score,
                warnings=[]
            )
            
            # logging.info(f"✅ Successfully processed antibody with score: {total_score:.3f}")
            logging.info(f"✅ Successfully processed antibody with score")
            # print(f"✅ Successfully processed antibody with score: {total_score:.3f}")
            print(f"✅ Successfully processed antibody with score")
            
            return ProcessingResult(
                sequence=sequence,
                metrics=antibody_metrics,
                success=True
            )
            
        except Exception as e:
            logging.error(f"❌ Error during processing: {str(e)}")
            print(f"❌ Error during processing: {str(e)}")
            return ProcessingResult(
                sequence=sequence,
                metrics=None,
                success=False,
                error=str(e)
            )

    def process_antibodies(self, input_file: str, output_dir: str = None) -> List[ProcessingResult]:
        """Process multiple antibody sequences from a file and save results to output directory"""
        try:
            # Create results directory
            results_dir = output_dir or os.path.join(os.path.dirname(input_file), "results")
            os.makedirs(results_dir, exist_ok=True)
            
            # Read all sequences
            with open(input_file, 'r') as f:
                lines = [line.strip() for line in f.readlines()]
            
            # Filter out empty lines and comments
            sequences = [line for line in lines if line and not line.startswith("#") and not line.startswith("//")]
            
            if not sequences:
                logging.warning(f"No valid sequences found in {input_file}")
                return []
                
            logging.info(f"Found {len(sequences)} antibody sequences in {input_file}")
            
            results = []
            success_count = 0
            error_count = 0
            
            # Process each antibody
            logging.info(f"🌲 Starting batch processing of {len(sequences)} antibodies...")
            print(f"\n🌲 Processing {len(sequences)} antibodies from {input_file}")
            
            for i, sequence in enumerate(sequences):
                try:
                    logging.info(f"Processing antibody {i+1}/{len(sequences)}:")
                    result = self.process_antibody(sequence)
                    results.append(result)
                    
                    if result.success:
                        success_count += 1
                        logging.info(f"✅ Successfully processed antibody {i+1} with score: {result.metrics.total_score:.3f}")
                        print(f"│  ├── {i+1}/{len(sequences)}: ✅ Success - Score: {result.metrics.total_score:.3f}")
                    else:
                        error_count += 1
                        logging.warning(f"❌ Failed to process antibody {i+1}: {result.error}")
                        print(f"│  ├── {i+1}/{len(sequences)}: ❌ Error - {result.error}")
                
                except Exception as e:
                    error_count += 1
                    error_msg = f"Unexpected error: {str(e)}"
                    logging.error(f"❌ Exception processing antibody {i+1}: {error_msg}")
                    print(f"│  ├── {i+1}/{len(sequences)}: ❌ Exception - {error_msg}")
                    # Create a failed result
                    results.append(ProcessingResult(
                        sequence=sequence,
                        metrics=None,
                        success=False,
                        error=error_msg
                    ))
            
            # Print summary tree
            print("│")
            print(f"└── 📊 Summary: {success_count} successful, {error_count} failed")
            
            for i, result in enumerate(results):
                status = "✅" if result.success else "❌"
                detail = f"Score: {result.metrics.total_score:.3f}" if result.success else f"Error: {result.error}"
                print(f"    ├── {i+1}: {status} {detail}")
            
            logging.info(f"Completed processing {len(sequences)} antibodies")
            logging.info(f"✅ Success: {success_count}, ❌ Failed: {error_count}")
            
            # Save results
            self.save_results(results, results_dir)
            logging.info(f"💾 Results saved to {results_dir}")
            
            return results
            
        except Exception as e:
            logging.error(f"❌ Error in batch processing: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
            
            # Return whatever results we've got so far, if any
            return results if 'results' in locals() else []
    
    def save_result(self, metrics: AntibodyMetrics, output_dir: str):
        """Save individual antibody metrics to a file"""
        
        # Convert metrics to dictionary and save
        metrics_dict = metrics.to_dict()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Correctly construct the directory path
        json_dir = os.path.join(output_dir, "json")

        # Create the directory if it doesn't exist
        os.makedirs(json_dir, exist_ok=True)

        output_file = os.path.join(json_dir, f"{metrics_dict['molecular_formula']}_{timestamp}.json")
        pd.DataFrame([metrics_dict]).to_json(output_file, orient='records', indent=2)
    
    def save_summary(self, results: List[ProcessingResult], output_dir: str):
        """Save summary of all antibody processing results"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        summary_file = os.path.join(output_dir, f"antibody_validation_summary_{timestamp}.txt")
        
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        
        with open(summary_file, 'w') as f:
            f.write("Antibody Validation Summary\n")
            f.write("=========================\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("Overall Statistics:\n")
            f.write("-----------------\n")
            f.write(f"Total antibodies processed: {len(results)}\n")
            f.write(f"Successfully processed: {len(successful)}\n")
            f.write(f"Failed: {len(failed)}\n\n")
            
            if successful:
                f.write("Score Statistics:\n")
                f.write("----------------\n")
                scores = [r.metrics.total_score for r in successful]
                f.write(f"Average score: {np.mean(scores):.4f}\n")
                f.write(f"Maximum score: {max(scores):.4f}\n")
                f.write(f"Minimum score: {min(scores):.4f}\n\n")
                
                f.write("Top 10 Highest Scoring Antibodies:\n")
                f.write("--------------------------------\n")
                top_10 = sorted(successful, key=lambda x: x.metrics.total_score, reverse=True)[:10]
                for result in top_10:
                    f.write(f"Sequence: {result.sequence[:20]}...\n")
                    f.write(f"Score: {result.metrics.total_score:.4f}\n")
                    f.write("-" * 50 + "\n")
        
        if failed:
            f.write("\nFailed Antibodies:\n")
            f.write("-----------------\n")
            for result in failed:
                    f.write(f"Sequence: {result.sequence[:20]}...\n")
                    f.write(f"Error: {result.error}\n")
                    f.write("-" * 50 + "\n")

    def get_successful_metrics(self, results: List[ProcessingResult]) -> List[AntibodyMetrics]:
        """Extract only successful metrics from processing results"""
        successful = []
        failed = []
        
        for i, result in enumerate(results):
            if result.success and result.metrics is not None:
                successful.append(result.metrics)
            else:
                error_msg = result.error if result.error else "Unknown error"
                seq_preview = result.sequence[:30] + "..." if result.sequence and len(result.sequence) > 30 else result.sequence
                failed.append((i, seq_preview, error_msg))
        
        # Log information about failed antibodies
        if failed:
            logging.warning(f"❌ {len(failed)} antibodies failed processing:")
            for i, seq, error in failed:
                logging.warning(f"  - Antibody #{i+1} (seq: {seq}): {error}")
        
        logging.info(f"✅ Successfully processed {len(successful)} out of {len(results)} antibodies")
        return successful

    def validate_antibodies(self, sequence_list: List[str]) -> List[AntibodyMetrics]:
        """
        Validate a list of antibody sequences and return valid antibodies with metrics.
        
        Args:
            sequence_list: List of antibody sequences to validate
            
        Returns:
            List of AntibodyMetrics objects for valid antibodies
        """
        print(f"\nStarting validation of {len(sequence_list)} antibodies...")
        valid_metrics = []
        invalid_count = 0
        
        for idx, sequence in enumerate(sequence_list, 1):
            try:
                print(f"\nValidating antibody {idx}/{len(sequence_list)}")
                print(f"├── Sequence: {sequence[:20]}...")
                
                # Basic sequence validation
                if not sequence or len(sequence) < 10:
                    print(f"└── Invalid sequence: sequence is too short or empty")
                    invalid_count += 1
                    continue
                
                # Calculate basic properties
                print(f"├── Calculating basic properties...")
                analyzed_seq = ProteinAnalysis(sequence)
                molecular_weight = analyzed_seq.molecular_weight()
                molecular_formula = f"C{len(sequence)}H{len(sequence)*2}N{len(sequence)}O{len(sequence)}"
                
                # Initialize metrics container with basic properties
                current_antigen_id = self.target_antigen_pdb_id or "Unknown"
                metrics_obj = AntibodyMetrics(
                    sequence=sequence,
                    antigen=current_antigen_id, # Use PDB ID as antigen identifier
                    molecular_weight=molecular_weight,
                    molecular_formula=molecular_formula,
                    blast={},
                    protparam={},
                    immunogenicity={},
                    stability={},
                    aggregation={},
                    glycosylation={},
                    binding_affinity={},
                    structure={},
                    epitope={},
                    developability={},
                    conservancy={},
                    weighted_scores={},
                    total_score=0.0,
                    warnings=[]
                )
                
                # Calculate all metrics
                print(f"├── Running BLAST analysis...")
                blast_result = perform_blast(sequence)
                metrics_obj.blast = blast_result.get('blast_result', {})
                
                print(f"├── Running ProtParam analysis...")
                protparam_result = analyze_with_protparam(sequence)
                metrics_obj.protparam = protparam_result.get('protein_params', {})
                
                print(f"├── Running immunogenicity prediction...")
                metrics_obj.immunogenicity = predict_immunogenicity(sequence)
                
                print(f"├── Running stability prediction...")
                stability_result = predict_stability(sequence)
                metrics_obj.stability = stability_result.get('stability_result', {})
                
                print(f"├── Running aggregation prediction...")
                aggregation_result = predict_aggregation(sequence)
                metrics_obj.aggregation = aggregation_result.get('aggregation_result', {})
                
                print(f"├── Running glycosylation prediction...")
                metrics_obj.glycosylation = predict_glycosylation(sequence)
                
                print(f"├── Running structure prediction...")
                structure_prediction_result = predict_structure(sequence, filename=molecular_formula, directory=self.pdb_files_path)
                metrics_obj.structure = structure_prediction_result
                
                print(f"├── Running binding affinity prediction...")
                if self.target_antigen_pdb_path and structure_prediction_result.get("pdb_file_path"):
                    binding_result = predict_binding_affinity(structure_prediction_result['pdb_file_path'], self.target_antigen_pdb_path)
                    metrics_obj.binding_affinity = binding_result
                else:
                    metrics_obj.binding_affinity = {"error": "Antigen PDB or antibody PDB path not available for validation method"}
                    logging.warning(f"Skipping binding affinity in validate_antibodies for {sequence[:20]}")

                print(f"├── Running epitope prediction...")
                metrics_obj.epitope = predict_bcell_epitopes(sequence)
                
                print(f"├── Running conservancy prediction...")
                metrics_obj.conservancy = predict_conservancy(sequence, metrics_obj.epitope.get('list_of_epitopes', []))
                
                print(f"├── Running developability prediction...")
                metrics_obj.developability = predict_developability(sequence)
                
                # Calculate weighted scores
                print(f"├── Calculating weighted scores...")
                # Need to pass the metrics from metrics_obj to calculate_weighted_scores
                temp_metrics_for_scoring = {
                    'binding_affinity': metrics_obj.binding_affinity,
                    'structure': metrics_obj.structure,
                    'glycosylation': metrics_obj.glycosylation,
                    'aggregation': metrics_obj.aggregation,
                    'protparam': metrics_obj.protparam,
                    'immunogenicity': metrics_obj.immunogenicity,
                    'conservancy': metrics_obj.conservancy,
                    'stability': metrics_obj.stability,
                    'epitope': metrics_obj.epitope,
                    'developability': metrics_obj.developability
                }
                weighted_scores, total_score = self.calculate_weighted_scores(temp_metrics_for_scoring)
                metrics_obj.weighted_scores = weighted_scores
                metrics_obj.total_score = total_score
                
                valid_metrics.append(metrics_obj)
                print(f"└── Validation successful")
                print(f"    ├── Formula: {metrics_obj.molecular_formula}")
                print(f"    ├── Weight: {metrics_obj.molecular_weight:.2f} Da")
                print(f"    └── Total Score: {metrics_obj.total_score:.4f}")
                
            except Exception as e:
                print(f"└── Error during validation: {str(e)}")
                invalid_count += 1
                continue
        
        print(f"\nValidation complete:")
        print(f"├── Total antibodies: {len(sequence_list)}")
        print(f"├── Valid antibodies: {len(valid_metrics)}")
        print(f"└── Invalid antibodies: {invalid_count}\n")
        
        return valid_metrics

    def save_results(self, results: List[ProcessingResult], output_dir: str):
        """
        Save all results to the output directory, handling any errors gracefully.
        
        Args:
            results: List of processing results
            output_dir: Directory to save results
        """
        try:
            # Create output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)
            
            # Save summary first
            self.save_summary(results, output_dir)
            
            # Save individual successful results
            successful_count = 0
            for result in results:
                if result.success and result.metrics:
                    try:
                        self.save_result(result.metrics, output_dir)
                        successful_count += 1
                    except Exception as e:
                        logging.error(f"❌ Error saving result for {result.sequence[:20]}...: {str(e)}")
            
            logging.info(f"💾 Saved {successful_count} individual results to {output_dir}")
            
            # Save failed antibodies to a separate file
            failed_results = [r for r in results if not r.success]
            if failed_results:
                failed_file = os.path.join(output_dir, "failed_antibodies_details.json")
                with open(failed_file, 'w') as f:
                    import json
                    # Convert to serializable format
                    failed_data = [{
                        "sequence": r.sequence,
                        "error": r.error
                    } for r in failed_results]
                    json.dump(failed_data, f, indent=2)
                logging.info(f"❌ Saved {len(failed_results)} failed antibody details to {failed_file}")
                
        except Exception as e:
            logging.error(f"❌ Error saving results: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())

if __name__ == "__main__":
    # Example usage
    validator = AntibodyValidator()
    results = validator.process_antibodies("antibodies.txt", "validation_results")
