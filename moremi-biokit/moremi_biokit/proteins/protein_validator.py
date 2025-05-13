"""
Antibody Metrics Collector V3

This module implements a comprehensive metrics collection system for proteins based on 
the metrics breakdown document. This version focuses on collecting all metrics for validation
and reporting purposes, while also calculating weighted scores for ranking.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Any, Union
import pandas as pd
from enum import Enum
import os
import logging
from datetime import datetime
from moremi_biokit import pdb_fetcher
from moremi_biokit.connectors import rcsb
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
    Reference ranges for protein properties based on the comprehensive metrics breakdown.
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
class ProteinMetrics:
    """Container for all calculated metrics for an protein"""
    sequence: str
    antigen: str
    antigen_id: Optional[str]
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
            'antigen_id': self.antigen_id,
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
    """Container for protein processing results"""
    sequence: str
    metrics: Optional[ProteinMetrics]
    error: Optional[str] = None
    success: bool = True

    def __str__(self) -> str:
        if self.success:
            return f"Success: {self.sequence[:20]}..."
        return f"Failed: {self.sequence[:20]}... - Error: {self.error}"

class ProteinValidator:
    """
    Enhanced protein metrics collector that gathers comprehensive metrics
    and calculates weighted scores for ranking.
    """
    
    def __init__(self, 
                 ranges: Optional[MetricRanges] = None, 
                 pdb_files_path: Optional[str] = None, # Directory to store PDBs generated by predict_structure
                 target_antigen_sequence: Optional[str] = None,
                 target_antigen_pdb_file_path: Optional[str] = None, # Direct path to a local antigen PDB file
                 target_antigen_pdb_chain_id: Optional[str] = None, # e.g., "1XYZ_A" for download and seq fetch
                 target_antigen_pdb_id_input: Optional[str] = None, # Legacy param for PDB ID only (for download)
                 antigen_pdb_download_path: Optional[str] = None # Directory to store downloaded antigen PDBs
                 ):
        """Initialize validator with metric ranges and flexible antigen specification.
        
        Args:
            ranges (Optional[MetricRanges]): Metric ranges for scoring.
            pdb_files_path (Optional[str]): Directory to store PDB files generated by
                the structure prediction tools (e.g., ESMFold).
            target_antigen_sequence (Optional[str]): Directly provided antigen amino acid sequence.
            target_antigen_pdb_file_path (Optional[str]): Path to a local PDB file for the antigen.
            target_antigen_pdb_chain_id (Optional[str]): Antigen identifier as 'PDBID_CHAIN' (e.g., "1XYZ_A").
                                                        The PDB file will be fetched, and sequence for the chain extracted.
            target_antigen_pdb_id_input (Optional[str]): PDB ID of the target antigen for PDB file fetching.
                                                        Considered if other specific antigen params are not set.
            antigen_pdb_download_path (Optional[str]): Directory to store downloaded antigen PDB files.
                                                       If not provided, downloads might be attempted to a default
                                                       location or fail if `pdb_fetcher.fetch_pdb` requires an output_dir.
        """
        self.ranges = ranges or MetricRanges()
        self.pdb_files_path = pdb_files_path # For ESMFold generated PDBs
        self.antigen_pdb_download_path = antigen_pdb_download_path # For downloaded antigen PDBs

        # Initialize core antigen properties
        self.target_antigen_sequence: Optional[str] = target_antigen_sequence
        self.target_antigen_pdb_path: Optional[str] = None
        self.target_antigen_pdb_id: Optional[str] = None
        self.target_antigen_chain_id: Optional[str] = None

        pdb_id_to_materialize = None

        # 1. Prioritize direct PDB file path for antigen
        if target_antigen_pdb_file_path:
            if os.path.isfile(target_antigen_pdb_file_path):
                self.target_antigen_pdb_path = target_antigen_pdb_file_path
                pdb_file_basename = os.path.basename(target_antigen_pdb_file_path)
                self.target_antigen_pdb_id = pdb_file_basename.split('.')[0].lower()
                logging.info(f"Using provided local antigen PDB file: {self.target_antigen_pdb_path} (Infered ID: {self.target_antigen_pdb_id})")
                if self.target_antigen_sequence:
                    logging.info(f"Using provided target_antigen_sequence for PDB file {self.target_antigen_pdb_id}.")
                else:
                    logging.warning(f"Local antigen PDB file {self.target_antigen_pdb_path} provided, but no corresponding target_antigen_sequence. Sequence-dependent metrics might be affected if sequence is not derivable from chain ID.")
                
            else:
                logging.warning(f"Provided target_antigen_pdb_file_path '{target_antigen_pdb_file_path}' is not a valid file. Ignoring.")

        # 2. If no local PDB path, check for PDBID_CHAIN for download/seq fetch
        elif target_antigen_pdb_chain_id:
            try:
                pdb_id_val, chain_id_val = target_antigen_pdb_chain_id.split('_', 1)
                self.target_antigen_pdb_id = pdb_id_val.lower()
                self.target_antigen_chain_id = chain_id_val.upper()
                pdb_id_to_materialize = self.target_antigen_pdb_id
                logging.info(f"Antigen specified by PDB_CHAIN_ID: {self.target_antigen_pdb_id}_{self.target_antigen_chain_id}")

                if self.target_antigen_sequence is None:
                    logging.info(f"Fetching antigen sequence for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id}...")
                    seq_details = rcsb.get_pdb_chain_sequence_details(self.target_antigen_pdb_id, self.target_antigen_chain_id)
                    if seq_details and seq_details.get('sequence'):
                        self.target_antigen_sequence = seq_details['sequence']
                        logging.info(f"Successfully fetched antigen sequence (len: {len(self.target_antigen_sequence)}) for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id}.")
                    else:
                        logging.warning(f"Could not fetch antigen sequence for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id}.")
                else:
                    logging.info(f"Using provided target_antigen_sequence for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id}.")
            except ValueError:
                logging.error(f"Invalid format for target_antigen_pdb_chain_id: '{target_antigen_pdb_chain_id}'. Expected 'PDBID_CHAIN'. Parameter ignored.")
        
        # 3. Else if legacy PDB ID input is given (and no local file or PDB_CHAIN_ID was processed)
        elif target_antigen_pdb_id_input:
            self.target_antigen_pdb_id = target_antigen_pdb_id_input.lower()
            pdb_id_to_materialize = self.target_antigen_pdb_id
            logging.info(f"Antigen PDB ID for file materialization set from target_antigen_pdb_id_input: {self.target_antigen_pdb_id}")
            if self.target_antigen_sequence:
                 logging.info(f"Using provided target_antigen_sequence for PDB ID {self.target_antigen_pdb_id}.")
            else:
                logging.warning(f"Antigen PDB ID {self.target_antigen_pdb_id} provided via input param, but no chain specified for sequence fetching, and no sequence provided directly. Antigen sequence is unknown.")

        # 4. Materialize PDB file if an ID for materialization is set and no local path was directly given
        if pdb_id_to_materialize and not self.target_antigen_pdb_path:
            self._materialize_antigen_pdb_file(pdb_id_to_materialize)

        # 5. Fallback to random internal PDB if no antigen PDB ID/path determined yet
        if not self.target_antigen_pdb_id and not self.target_antigen_pdb_path:
            logging.info("No specific target antigen PDB provided by user, attempting to use a random internal PDB ID.")
            internal_pdb_ids = pdb_fetcher.list_internal_pdb_ids()
            if internal_pdb_ids:
                random_pdb_id = random.choice(internal_pdb_ids)
                self.target_antigen_pdb_id = random_pdb_id # Set the ID
                logging.info(f"Randomly selected internal PDB ID for antigen: {self.target_antigen_pdb_id}")
                self._materialize_antigen_pdb_file(self.target_antigen_pdb_id)
                # For random PDB, sequence is unknown unless fetched, but we don't have a chain.
                if not self.target_antigen_sequence:
                    logging.warning(f"Using random PDB {self.target_antigen_pdb_id} as antigen. Antigen sequence is unknown as no chain was specified.")
            else:
                logging.warning("No internal PDBs found for default target antigen selection, and no specific antigen provided by user.")
        
        # Final logging of antigen state
        log_msg_pdb_id = self.target_antigen_pdb_id if self.target_antigen_pdb_id else "Not set/inferred"
        log_msg_pdb_path = self.target_antigen_pdb_path if self.target_antigen_pdb_path else "Not set"
        log_msg_chain_id = self.target_antigen_chain_id if self.target_antigen_chain_id else "Not set"
        log_msg_seq_status = 'Provided/Fetched' if self.target_antigen_sequence else 'Unknown'
        
        logging.info(
            f"ProteinValidator antigen setup: PDB_ID='{log_msg_pdb_id}', PDB_Path='{log_msg_pdb_path}', "
            f"Chain_ID='{log_msg_chain_id}', Sequence_Status='{log_msg_seq_status}'"
        )

    def _materialize_antigen_pdb_file(self, pdb_id_to_materialize: str):
        """
        Helper function to fetch/locate and set the PDB file path for a given antigen PDB ID.
        Uses pdb_fetcher.fetch_pdb for external fetching.
        The download location is determined by self.antigen_pdb_download_path.
        Sets self.target_antigen_pdb_path.
        """
        if not pdb_id_to_materialize:
            logging.warning("_materialize_antigen_pdb_file called with no PDB ID.")
            return

        logging.info(f"Attempting to materialize PDB file for antigen ID: {pdb_id_to_materialize} using pdb_fetcher.fetch_pdb")

        output_directory_for_antigen_pdb = self.antigen_pdb_download_path

        if not output_directory_for_antigen_pdb:
            logging.warning(
                f"antigen_pdb_download_path is not set in ProteinValidator. "
                f"PDB file for antigen ID '{pdb_id_to_materialize}' might not be downloaded if fetch_pdb requires an output_dir. "
                f"Attempting fetch without explicit output directory (may use internal cache or fail)."
            )
            # Attempt fetch without output_dir, relying on pdb_fetcher's internal logic or cache
            fetched_info = pdb_fetcher.fetch_pdb(
                target=pdb_id_to_materialize,
                input_type='identifier',
                use_internal_db=True, # Allow checking internal if no download path specified
                output_dir=None 
            )
            if fetched_info and fetched_info.get("status") == "success" and fetched_info.get("pdb_path"):
                self.target_antigen_pdb_path = fetched_info.get("pdb_path")
                logging.info(f"Using PDB path (potentially internal or pre-existing) for antigen '{pdb_id_to_materialize}': {self.target_antigen_pdb_path}")
            else:
                logging.warning(f"Could not obtain a PDB path for antigen ID '{pdb_id_to_materialize}' via fetch_pdb without an explicit antigen_pdb_download_path.")
            return

        try:
            os.makedirs(output_directory_for_antigen_pdb, exist_ok=True)
            logging.info(f"Ensured antigen PDB download path exists: {output_directory_for_antigen_pdb}")

            fetched_info = pdb_fetcher.fetch_pdb(
                target=pdb_id_to_materialize,
                input_type='identifier',
                source_db='rcsb', 
                use_internal_db=False, # Prioritize external fetch when download path is given
                output_dir=output_directory_for_antigen_pdb,
                db_file_format='pdb'
            )

            if fetched_info and fetched_info.get("status") == "success" and fetched_info.get("pdb_path"):
                self.target_antigen_pdb_path = fetched_info.get("pdb_path")
                logging.info(f"Successfully fetched/materialized PDB for antigen '{pdb_id_to_materialize}' to '{self.target_antigen_pdb_path}' using pdb_fetcher.fetch_pdb.")
            else:
                status_msg = fetched_info.get('status', 'unknown') if fetched_info else 'fetch_pdb call failed'
                error_msg = fetched_info.get('message', 'No specific error message.') if fetched_info else ''
                logging.warning(
                    f"Could not fetch/materialize PDB for antigen ID '{pdb_id_to_materialize}' using pdb_fetcher.fetch_pdb "
                    f"to '{output_directory_for_antigen_pdb}'. Status: {status_msg}. Message: {error_msg}"
                )
        except Exception as e:
            logging.error(f"Exception during PDB materialization for '{pdb_id_to_materialize}' to '{output_directory_for_antigen_pdb}': {e}")

    # TODO: will be deleted after some confirming some review 
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
    
    def process_protein(self, sequence: str) -> ProcessingResult:
        """Process a single protein sequence and calculate all metrics"""
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
            logging.info(f"\n🧪 Processing protein sequence: {sequence[:20]}...")
            logging.info(f"├── Calculating basic properties...")
            print(f"\n🧪 Processing protein sequence: {sequence[:20]}...")
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
            structure_result = predict_structure(sequence, output_pdb_filename_prefix=molecular_formula, output_directory=self.pdb_files_path)
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
                metrics['binding_affinity'] = {"error": "Receptor(Antigen) PDB or protein PDB path not available"}
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
            protein_metrics = ProteinMetrics(
                sequence=sequence,
                antigen=self.target_antigen_sequence or "", # Use fetched/provided antigen sequence
                antigen_id=self.target_antigen_pdb_id or "Unknown", # Use PDB ID of the antigen
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
            
            # logging.info(f"✅ Successfully processed protein with score: {total_score:.3f}")
            logging.info(f"✅ Successfully processed protein with score")
            # print(f"✅ Successfully processed protein with score: {total_score:.3f}")
            print(f"✅ Successfully processed protein with score")
            
            return ProcessingResult(
                sequence=sequence,
                metrics=protein_metrics,
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

    def process_proteins(self, input_file: str, output_dir: str = None) -> List[ProcessingResult]:
        """Process multiple protein sequences from a file (text or FASTA) and save results to output directory."""
        try:
            # Create results directory
            results_dir = output_dir or os.path.join(os.path.dirname(input_file), "results")
            os.makedirs(results_dir, exist_ok=True)
            
            sequences = []
            current_sequence_lines = []
            # Read all sequences from the input file
            with open(input_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("//"):
                        continue # Skip empty lines and comments
                    
                    if line.startswith(">"): # FASTA header
                        if current_sequence_lines: # Save previous sequence
                            sequences.append("".join(current_sequence_lines))
                            current_sequence_lines = []
                        # We can store or log the header if needed, for now, we just use it as a separator
                    else: # Sequence line
                        current_sequence_lines.append(line)
            
            if current_sequence_lines: # Add the last sequence in the file
                sequences.append("".join(current_sequence_lines))

            if not sequences:
                logging.warning(f"No valid sequences found in {input_file} (checked for plain and FASTA format)")
                return []
                
            logging.info(f"Found {len(sequences)} protein sequences in {input_file} (parsed as plain/FASTA format)")
            
            results = []
            success_count = 0
            error_count = 0
            
            # Process each protein
            logging.info(f"🌲 Starting batch processing of {len(sequences)} proteins...")
            print(f"\n🌲 Processing {len(sequences)} proteins from {input_file}")
            
            for i, sequence in enumerate(sequences):
                try:
                    logging.info(f"Processing protein {i+1}/{len(sequences)}:")
                    result = self.process_protein(sequence)
                    results.append(result)
                    
                    if result.success:
                        success_count += 1
                        logging.info(f"✅ Successfully processed protein {i+1} with score: {result.metrics.total_score:.3f}")
                        print(f"│  ├── {i+1}/{len(sequences)}: ✅ Success - Score: {result.metrics.total_score:.3f}")
                    else:
                        error_count += 1
                        logging.warning(f"❌ Failed to process protein {i+1}: {result.error}")
                        print(f"│  ├── {i+1}/{len(sequences)}: ❌ Error - {result.error}")
                
                except Exception as e:
                    error_count += 1
                    error_msg = f"Unexpected error: {str(e)}"
                    logging.error(f"❌ Exception processing protein {i+1}: {error_msg}")
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
            
            logging.info(f"Completed processing {len(sequences)} proteins")
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
    
    # TODO: this so will be deleted after some confirmations
    def save_result(self, metrics: ProteinMetrics, output_dir: str):
        """Save individual protein metrics to a file"""
        
        # Convert metrics to dictionary and save
        metrics_dict = metrics.to_dict()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Correctly construct the directory path
        json_dir = os.path.join(output_dir, "json")

        # Create the directory if it doesn't exist
        os.makedirs(json_dir, exist_ok=True)

        output_file = os.path.join(json_dir, f"{metrics_dict['molecular_formula']}_{timestamp}.json")
        pd.DataFrame([metrics_dict]).to_json(output_file, orient='records', indent=2)
    
    # TODO: DELETE DO SOMETHING TO THIS GUY
    def save_summary(self, results: List[ProcessingResult], output_dir: str):
        """Save summary of all protein processing results"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        summary_file = os.path.join(output_dir, f"protein_validation_summary_{timestamp}.txt")
        
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        
        with open(summary_file, 'w') as f:
            f.write("Antibody Validation Summary\n")
            f.write("=========================\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("Overall Statistics:\n")
            f.write("-----------------\n")
            f.write(f"Total proteins processed: {len(results)}\n")
            f.write(f"Successfully processed: {len(successful)}\n")
            f.write(f"Failed: {len(failed)}\n\n")
            
            if successful:
                f.write("Score Statistics:\n")
                f.write("----------------\n")
                scores = [r.metrics.total_score for r in successful]
                f.write(f"Average score: {np.mean(scores):.4f}\n")
                f.write(f"Maximum score: {max(scores):.4f}\n")
                f.write(f"Minimum score: {min(scores):.4f}\n\n")
                
                f.write("Top 10 Highest Scoring proteins:\n")
                f.write("--------------------------------\n")
                top_10 = sorted(successful, key=lambda x: x.metrics.total_score, reverse=True)[:10]
                for result in top_10:
                    f.write(f"Sequence: {result.sequence[:20]}...\n")
                    f.write(f"Score: {result.metrics.total_score:.4f}\n")
                    f.write("-" * 50 + "\n")
        
        if failed:
            f.write("\nFailed proteins:\n")
            f.write("-----------------\n")
            for result in failed:
                    f.write(f"Sequence: {result.sequence[:20]}...\n")
                    f.write(f"Error: {result.error}\n")
                    f.write("-" * 50 + "\n")

    def get_successful_metrics(self, results: List[ProcessingResult]) -> List[ProteinMetrics]:
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
        
        # Log information about failed proteins
        if failed:
            logging.warning(f"❌ {len(failed)} proteins failed processing:")
            for i, seq, error in failed:
                logging.warning(f"  - Antibody #{i+1} (seq: {seq}): {error}")
        
        logging.info(f"✅ Successfully processed {len(successful)} out of {len(results)} proteins")
        return successful
    
    def validate_proteins(self, sequence_list: Union[str, List[str]]) -> List[ProteinMetrics]:
        """
        Validate a list of protein sequences or a single protein sequence string 
        by calling process_protein for each and returns a list of ProteinMetrics for successful validations.
        
        Args:
            sequence_list: A single protein sequence string or a list of protein sequences to validate.
            
        Returns:
            List of ProteinMetrics objects for valid proteins.
        """
        
        sequences_to_process: List[str]
        if isinstance(sequence_list, str):
            sequences_to_process = [sequence_list]
        elif isinstance(sequence_list, list):
            sequences_to_process = sequence_list
        else:
            logging.error("Invalid input type for sequence_list. Expected str or List[str].")
            print("Error: Invalid input for sequence_list. Must be a string or a list of strings.")
            return []

        print(f"\nStarting validation of {len(sequences_to_process)} proteins by leveraging process_protein...")
        valid_metrics_list = []
        invalid_count = 0
        
        for idx, sequence in enumerate(sequences_to_process, 1):
            print(f"\nValidating protein {idx}/{len(sequences_to_process)}: {sequence[:20]}...")
            result = self.process_protein(sequence) # Calls the main processing logic
            
            if result.success and result.metrics:
                valid_metrics_list.append(result.metrics)
                # process_protein already logs success, so we can keep this minimal
                print(f"└── Validation successful for protein {idx} (via process_protein). Score: {result.metrics.total_score:.3f}")
            else:
                invalid_count += 1
                # process_protein already logs failure
                print(f"└── Validation failed for protein {idx} (via process_protein). Error: {result.error}")
        
        print(f"\nValidation complete:")
        print(f"├── Total proteins attempted: {len(sequences_to_process)}")
        print(f"├── Successfully validated (metrics collected): {len(valid_metrics_list)}")
        print(f"└── Failed validation: {invalid_count}\n")
        
        return valid_metrics_list

    def get_metrics_as_dict_list(self, results: List[ProcessingResult]) -> List[Dict[str, Any]]:
        """
        Converts ProteinMetrics from a list of successful ProcessingResult objects into a list of dictionaries.

        Args:
            results: A list of ProcessingResult objects.

        Returns:
            A list of dictionaries, where each dictionary represents the metrics of a successfully processed protein.
        """
        dict_list = []
        for result in results:
            if result.success and result.metrics:
                try:
                    dict_list.append(result.metrics.to_dict())
                except Exception as e:
                    logging.error(f"Error converting metrics to dict for sequence {result.sequence[:20]}...: {e}")
        return dict_list

    def save_metrics_to_csv(self, results: List[ProcessingResult], output_csv_path: str):
        """
        Saves the ProteinMetrics from successful ProcessingResult objects to a CSV file.

        Args:
            results: A list of ProcessingResult objects.
            output_csv_path: The path (including filename) to save the CSV file.
        """
        if not results:
            logging.warning("No results provided to save_metrics_to_csv.")
            print("Warning: No results to save to CSV.")
            return

        metrics_dict_list = self.get_metrics_as_dict_list(results)
        
        if not metrics_dict_list:
            logging.warning("No successful metrics found to save to CSV.")
            print("Warning: No successful metrics to save to CSV.")
            return

        try:
            df = pd.DataFrame(metrics_dict_list)
            
            # Flatten nested dictionaries for better CSV representation if necessary
            # For example, if 'metrics' is a nested dict:
            if not df.empty and 'metrics' in df.columns:
                # This is a basic way to flatten; more sophisticated flattening might be needed
                # depending on the depth and structure of nested dicts.
                try:
                    metrics_expanded = pd.json_normalize(df['metrics'], sep='_')
                    df = df.drop(columns=['metrics'])
                    df = pd.concat([df, metrics_expanded], axis=1)
                except Exception as e:
                    logging.error(f"Could not fully flatten 'metrics' column for CSV: {e}. Proceeding with unflattened parts.")

            # Ensure output directory exists
            output_dir = os.path.dirname(output_csv_path)
            if output_dir: # If output_csv_path includes a directory
                os.makedirs(output_dir, exist_ok=True)
            
            df.to_csv(output_csv_path, index=False)
            logging.info(f"Successfully saved metrics for {len(metrics_dict_list)} proteins to CSV: {output_csv_path}")
            print(f"INFO: Metrics for {len(metrics_dict_list)} proteins saved to {output_csv_path}")
        except Exception as e:
            logging.error(f"Error saving metrics to CSV at {output_csv_path}: {e}")
            print(f"ERROR: Could not save metrics to CSV: {e}")

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
            
            # Save failed proteins to a separate file
            failed_results = [r for r in results if not r.success]
            if failed_results:
                failed_file = os.path.join(output_dir, "failed_proteins_details.json")
                with open(failed_file, 'w') as f:
                    import json
                    # Convert to serializable format
                    failed_data = [{
                        "sequence": r.sequence,
                        "error": r.error
                    } for r in failed_results]
                    json.dump(failed_data, f, indent=2)
                logging.info(f"❌ Saved {len(failed_results)} failed protein details to {failed_file}")
                
        except Exception as e:
            logging.error(f"❌ Error saving results: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())

if __name__ == "__main__":
    # Example usage
    # validator = ProteinValidator() # Original
    # Example with new parameters:
    # validator = ProteinValidator(pdb_files_path="validation_results/generated_pdbs", antigen_pdb_download_path="validation_results/downloaded_antigen_pdbs", target_antigen_pdb_chain_id="1XYZ_A")
    # validator = ProteinValidator(target_antigen_pdb_file_path="path/to/local/antigen.pdb", target_antigen_sequence="ANTIGENSEQUENCE")
    # validator = ProteinValidator(target_antigen_sequence="ANTIGENSEQUENCEONLY") # Might affect binding affinity if PDB path not derivable
    
    # Setup basic logging for testing
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Test case 1: Using PDB_CHAIN_ID
    print("\n--- Test Case 1: PDB_CHAIN_ID ---")
    validator1 = ProteinValidator(
        pdb_files_path="validation_results/generated_pdbs_test1", 
        antigen_pdb_download_path="validation_results/downloaded_antigen_pdbs_test1",
        target_antigen_pdb_chain_id="4R19_A" # Example PDB and chain
    )
    if validator1.target_antigen_sequence:
        print(f"Antigen sequence fetched for 4R19_A: {validator1.target_antigen_sequence[:30]}...")
    if validator1.target_antigen_pdb_path:
        print(f"Antigen PDB path for 4R19_A: {validator1.target_antigen_pdb_path}")
    # results1 = validator1.process_proteins("proteins.txt", "validation_results/test1_output") # Assuming proteins.txt exists

    # Test case 2: Using local PDB file and sequence
    print("\n--- Test Case 2: Local PDB and Sequence ---")
    # Create a dummy antigen.pdb for testing this case
    dummy_pdb_dir = "validation_results/dummy_antigen_local" # This path is for the *local* antigen, not a download path
    os.makedirs(dummy_pdb_dir, exist_ok=True)
    dummy_pdb_path = os.path.join(dummy_pdb_dir, "6M0J.pdb") 
    if not os.path.exists(dummy_pdb_path):
         with open(dummy_pdb_path, "w") as f:
             f.write("HEADER    DUMMY PDB FILE FOR TESTING 6M0J\n")
             f.write("ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\n")
    
    validator2 = ProteinValidator(
        pdb_files_path="validation_results/generated_pdbs_test2", # Separate path for generated pdbs
        # antigen_pdb_download_path is not needed here as we provide a direct file path
        target_antigen_pdb_file_path=dummy_pdb_path,
        target_antigen_sequence="MTQVPSNPPPVVGARHNFSLKECGF" 
    )
    print(f"Antigen sequence (provided) for local PDB: {validator2.target_antigen_sequence[:30]}...")
    print(f"Antigen PDB path (local): {validator2.target_antigen_pdb_path}")
    print(f"Antigen PDB ID (inferred from local): {validator2.target_antigen_pdb_id}")
    # results2 = validator2.process_proteins("proteins.txt", "validation_results/test2_output")

    # Test case 3: Only antigen sequence
    print("\n--- Test Case 3: Antigen Sequence Only ---")
    validator3 = ProteinValidator(target_antigen_sequence="JUSTASEQUENCEFORTESTING")
    print(f"Antigen sequence (provided): {validator3.target_antigen_sequence}")
    print(f"Antigen PDB path: {validator3.target_antigen_pdb_path}") # Expected to be None
    # results3 = validator3.process_proteins("proteins.txt", "validation_results/test3_output")

    # Test case 4: Legacy PDB ID input
    print("\n--- Test Case 4: Legacy PDB ID Input ---")
    validator4 = ProteinValidator(
        pdb_files_path="validation_results/generated_pdbs_test4",
        antigen_pdb_download_path="validation_results/downloaded_antigen_pdbs_test4",
        target_antigen_pdb_id_input="1A2Y" # Example PDB ID
    )
    print(f"Antigen PDB ID (from legacy input): {validator4.target_antigen_pdb_id}")
    print(f"Antigen PDB path: {validator4.target_antigen_pdb_path}")
    print(f"Antigen sequence: {validator4.target_antigen_sequence}") # Expected to be None or user-provided
    # results4 = validator4.process_proteins("proteins.txt", "validation_results/test4_output")

    # Test case 5: Default (random internal PDB)
    print("\n--- Test Case 5: Default (Random Internal PDB) ---")
    # Assuming internal PDBs are set up for pdb_fetcher.list_internal_pdb_ids() to work
    # validator5 = ProteinValidator(
    #     pdb_files_path="validation_results/generated_pdbs_test5",
    #     antigen_pdb_download_path="validation_results/downloaded_antigen_pdbs_test5" # Path if it attempts download for random
    # )
    # print(f"Antigen PDB ID (random): {validator5.target_antigen_pdb_id}")
    pass # End of __main__ example updates
