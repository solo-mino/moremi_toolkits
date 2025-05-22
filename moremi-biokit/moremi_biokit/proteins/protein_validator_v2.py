"""
Antibody Metrics Collector V3 - Refactored (V2)

This module implements a comprehensive metrics collection system for proteins based on 
the metrics breakdown document. This version focuses on collecting all metrics for validation
and reporting purposes, while also calculating weighted scores for ranking.
This V2 aims to improve modularity and efficiency, particularly in antigen handling.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any, Union
from pathlib import Path
import pandas as pd
from enum import Enum
import os
import logging
from datetime import datetime
from moremi_biokit.connectors import rcsb
import json # Added import for json

from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
    STRUCTURE = "Structure"
    BINDING_AFFINITY = "Binding Affinity"
    EPITOPE = "Epitope"
    CONSERVANCY = "Conservancy"
    DEVELOPABILITY = "Developability"
    
    
@dataclass
class ProteinMetrics:
    """Container for all calculated metrics for an protein"""
    sequence: str
    antigen: str # This might become Optional or be set later
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
    structure: Dict[str, Any]
    binding_affinity: Dict[str, Any]
    epitope: Dict[str, Any]
    conservancy: Dict[str, Any]
    developability: Dict[str, Any]
    
    # These scores will be calculated by the ranker, not the validator
    weighted_scores: Dict[str, float] = field(default_factory=dict)
    total_score: float = 0.0
    
    # Additional antigen information
    antigen_pdb_chain_id: Optional[str] = None
    
    warnings: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        """Convert metrics to dictionary format"""
        return {
            'sequence': self.sequence,
            'antigen': self.antigen,
            'antigen_id': self.antigen_id,
            'antigen_pdb_chain_id': self.antigen_pdb_chain_id,
            'molecular_weight': self.molecular_weight,
            'molecular_formula': self.molecular_formula,
            'metrics': {
                'blast': self.blast,
                'protparam': self.protparam,
                'immunogenicity': self.immunogenicity,
                'stability': self.stability,
                'aggregation': self.aggregation,
                'glycosylation': self.glycosylation,
                'structure': self.structure,
                'binding_affinity': self.binding_affinity,
                'epitope': self.epitope,
                'conservancy': self.conservancy,
                'developability': self.developability
            },
            'weighted_scores': self.weighted_scores,
            'total_score': self.total_score,
            'warnings': self.warnings
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ProteinMetrics':
        """Reconstructs a ProteinMetrics object from a dictionary (e.g., from JSON)."""
        metrics_data = data.get('metrics', {})
        return cls(
            sequence=data.get('sequence', ''),
            antigen=data.get('antigen', 'N/A'),
            antigen_id=data.get('antigen_id'),
            molecular_weight=data.get('molecular_weight', 0.0),
            molecular_formula=data.get('molecular_formula', 'N/A'),
            antigen_pdb_chain_id=data.get('antigen_pdb_chain_id'),
            
            blast=metrics_data.get(MetricCategory.BLAST.value.lower(), metrics_data.get('blast', {})), # Handle old and new key naming
            protparam=metrics_data.get(MetricCategory.PROTPARAM.value.lower(), metrics_data.get('protparam', {})),
            immunogenicity=metrics_data.get(MetricCategory.IMMUNOGENICITY.value.lower(), metrics_data.get('immunogenicity', {})),
            stability=metrics_data.get(MetricCategory.STABILITY.value.lower(), metrics_data.get('stability', {})),
            aggregation=metrics_data.get(MetricCategory.AGGREGATION.value.lower(), metrics_data.get('aggregation', {})),
            glycosylation=metrics_data.get(MetricCategory.GLYCOSYLATION.value.lower(), metrics_data.get('glycosylation', {})),
            structure=metrics_data.get(MetricCategory.STRUCTURE.value.lower(), metrics_data.get('structure', {})),
            binding_affinity=metrics_data.get(MetricCategory.BINDING_AFFINITY.value.lower().replace(" ", "_"), metrics_data.get('binding_affinity', {})),
            epitope=metrics_data.get(MetricCategory.EPITOPE.value.lower(), metrics_data.get('epitope', {})),
            conservancy=metrics_data.get(MetricCategory.CONSERVANCY.value.lower(), metrics_data.get('conservancy', {})),
            developability=metrics_data.get(MetricCategory.DEVELOPABILITY.value.lower(), metrics_data.get('developability', {})),
            
            weighted_scores=data.get('weighted_scores', {}),
            total_score=data.get('total_score', 0.0),
            warnings=data.get('warnings', [])
        )

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

class ProteinValidatorV2:
    """
    Enhanced protein metrics collector (V2) that gathers comprehensive metrics.
    This version defers antigen-specific setup until necessary.
    """
    
    def __init__(self, 
                 pdb_files_path: Optional[str] = None, # Directory to store PDBs generated by predict_structure for the *antibody*
                 metrics_to_run: Optional[List[MetricCategory]] = None # Specify which metric categories to run
                 ):
        """Initialize validator with metric ranges. Antigen setup is deferred.
        
        Args:
            pdb_files_path (Optional[str]): Directory to store PDB files generated by
                the structure prediction tools (e.g., ESMFold) for the antibody.
            metrics_to_run (Optional[List[MetricCategory]]): A list of MetricCategory enums to specify which
                                                             metric calculations should be performed. If None, all metrics are run.
        """
        self.pdb_files_path = pdb_files_path # For ESMFold generated PDBs (for the antibody)
        self.metrics_to_run = metrics_to_run

        # Antigen related properties - to be set by a dedicated method if antigen-dependent metrics are run.
        self.target_antigen_sequence: Optional[str] = None
        self.target_antigen_pdb_path: Optional[str] = None # Path to the antigen's PDB file (local, downloaded, or predicted)
        self.target_antigen_pdb_id: Optional[str] = None   # PDB ID of the antigen
        self.target_antigen_chain_id: Optional[str] = None # Chain ID for the antigen PDB
        self.antigen_pdb_download_path: Optional[str] = None # Path where antigen PDBs might be downloaded by set_antigen_context

        logging.info(
            f"ProteinValidatorV2 initialized. Antibody PDB output path: {self.pdb_files_path}. "
            f"Metrics to run: {'All' if not metrics_to_run else [m.value for m in metrics_to_run]}. "
            f"Antigen context will be set separately if needed."
        )

    def _predict_antigen_structure(self, antigen_sequence: str, output_directory: str):
        """
        Predicts the structure of the antigen sequence using structure_predictor.
        Sets self.target_antigen_pdb_path if prediction is successful.

        Args:
            antigen_sequence (str): The antigen sequence to predict.
            output_directory (str): The directory to save the predicted PDB file.

        Returns:
            Optional[str]: Path to the predicted PDB file, or None if failed.
        """
        if not antigen_sequence:
            logging.warning("❌ Cannot predict antigen structure: No antigen sequence provided")
            return None
            
        logging.info(f"🦠🔬Predicting structure for antigen sequence (len: {len(antigen_sequence)})")
        
        os.makedirs(output_directory, exist_ok=True)
        
        try:
            from .analysis_tools import structure_predictor # Keep import local to method if specific
            
            result = structure_predictor.predict_structure(
                sequence=antigen_sequence,
                output_directory=output_directory,
                output_pdb_filename_prefix="antigen_structure_pred" # More specific prefix
            )
            
            if result and result.get("pdb_file_path"):
                logging.info(f"✅ Successfully predicted antigen structure. PDB file: {result['pdb_file_path']}")
                return result["pdb_file_path"]
            else:
                logging.warning(f"❌ Failed to predict antigen structure for sequence: {antigen_sequence[:30]}...")
                return None
        except Exception as e:
            logging.error(f"❌ Error predicting antigen structure: {str(e)}")
            return None

    def _materialize_antigen_pdb_file(self, pdb_id_to_materialize: str, download_path_for_antigen_pdb: str):
        """
        Helper function to fetch/locate and return the PDB file path for a given antigen PDB ID.
        The download location is determined by download_path_for_antigen_pdb.
        
        Args:
            pdb_id_to_materialize (str): The PDB ID to fetch.
            download_path_for_antigen_pdb (str): Directory to download/store the PDB.

        Returns:
            Optional[str]: Path to the PDB file, or None if failed.
        """
        if not pdb_id_to_materialize:
            logging.warning("⚠️ _materialize_antigen_pdb_file called with no PDB ID.")
            return None

        logging.info(f"🔍 Attempting to materialize PDB file for antigen ID: {pdb_id_to_materialize} into {download_path_for_antigen_pdb} using rcsb.download_pdb_from_rcsb")

        try:
            os.makedirs(download_path_for_antigen_pdb, exist_ok=True)
            #logging.info(f"Ensured antigen PDB download path exists: {download_path_for_antigen_pdb}") # rcsb.download_pdb_from_rcsb will also ensure path

            # Directly use rcsb.download_pdb_from_rcsb
            # It handles its own logging for success/failure and returns the path or None.
            downloaded_path = rcsb.download_pdb_from_rcsb(
                pdb_id=pdb_id_to_materialize,
                output_dir=download_path_for_antigen_pdb,
                file_format='pdb' # Assuming 'pdb' format is desired for antigen context
            )

            if downloaded_path:
                logging.info(f"✅ Successfully downloaded PDB for antigen '{pdb_id_to_materialize}' to '{downloaded_path}' via rcsb.download_pdb_from_rcsb.")
                return downloaded_path
            else:
                # download_pdb_from_rcsb would have printed its own error (404, etc.)
                logging.warning(
                    f"⚠️ rcsb.download_pdb_from_rcsb failed to download PDB for antigen ID '{pdb_id_to_materialize}'. "
                    f"Check previous logs from the download function for details."
                )
                return None
        except Exception as e:
            # Catch any unexpected errors from os.makedirs or if download_pdb_from_rcsb itself raises an unexpected exception
            logging.error(f"❌ Exception during PDB materialization for '{pdb_id_to_materialize}' to '{download_path_for_antigen_pdb}': {e}")
            return None

    # Placeholder for the new method to set antigen context
    def set_antigen_context(self,
                            target_antigen_sequence: Optional[str] = None,
                            target_antigen_pdb_file_path: Optional[str] = None,
                            target_pdb_id: Optional[str] = None, # New parameter for PDB ID
                            target_antigen_pdb_chain_id: Optional[str] = None, # e.g., "4R19_A"
                            antigen_pdb_download_dir: Optional[str] = None) -> bool:
        """
        Sets the antigen PDB structure context, optionally with its sequence, for metrics like binding affinity.
        The primary goal is to obtain a PDB file path for the antigen.

        Precedence for determining antigen PDB and sequence:
        1.  **Local PDB File (`target_antigen_pdb_file_path`)**:
            - Uses the provided local PDB file. PDB path is set.
            - Sequence is obtained if `target_antigen_sequence` is also given, or attempted via RCSB 
              if `target_antigen_pdb_chain_id` is provided (using inferred PDB ID from filename).
        2.  **PDB ID (`target_pdb_id`)**:
            - Attempts to download the PDB file from RCSB using `target_pdb_id`.
            - If download succeeds: PDB path is set.
            - Sequence is obtained if `target_antigen_sequence` is provided. If `target_antigen_pdb_chain_id`
              is also provided and its PDB ID part matches `target_pdb_id`, its chain part is used
              with the downloaded PDB ID to fetch sequence from RCSB.
        3.  **PDB ID and Chain (`target_antigen_pdb_chain_id`)**:
            - Parses PDB ID and chain from `target_antigen_pdb_chain_id`.
            - Attempts to fetch sequence using `rcsb.fetch_sequence_details_by_pdb_chain_id` (MoreMi microservice).
            - If sequence is obtained, its 3D structure is predicted. If prediction succeeds, PDB path is set.
        4.  **Antigen Sequence Only (`target_antigen_sequence`):**
            - If only a sequence is provided (and no PDB path was determined above, or sequence was derived
              from microservice but structure prediction failed), its 3D structure is predicted.
            - If prediction succeeds, PDB path is set.

        Args:
            target_antigen_sequence (Optional[str]): Antigen amino acid sequence.
            target_antigen_pdb_file_path (Optional[str]): Path to a local PDB file for the antigen.
            target_pdb_id (Optional[str]): A PDB ID (e.g., "4R19") to download from RCSB.
            target_antigen_pdb_chain_id (Optional[str]): 'PDBID_CHAIN' (e.g., "4R19_A") for MoreMi microservice interaction.
            antigen_pdb_download_dir (Optional[str]): Directory for downloaded/predicted antigen PDBs.

        Returns:
            bool: True if an antigen PDB file path (`self.target_antigen_pdb_path`) was successfully set, False otherwise.
                  `self.target_antigen_sequence` will be populated if available but is not strictly required for success.
        """
        logging.info("🔍 Attempting to set antigen context. Primary goal: obtain PDB path.")

        if antigen_pdb_download_dir:
            self.antigen_pdb_download_path = antigen_pdb_download_dir
        elif self.pdb_files_path:
            self.antigen_pdb_download_path = os.path.join(self.pdb_files_path, "antigen_pdbs")
        else:
            self.antigen_pdb_download_path = os.path.join(os.getcwd(), "antigen_pdbs_v2_default")
        os.makedirs(self.antigen_pdb_download_path, exist_ok=True)
        logging.info(f"🔍 Antigen PDBs will be managed in: {self.antigen_pdb_download_path}")

        # Initialize/reset antigen context attributes
        self.target_antigen_sequence = None
        self.target_antigen_pdb_path = None
        self.target_antigen_pdb_id = None # This will store the PDB ID of the final chosen antigen structure
        self.target_antigen_chain_id = None # This will store the chain if applicable

        # --- Priority 1: Local PDB File Path --- 
        if target_antigen_pdb_file_path:
            logging.info(f"Priority 1: Checking local antigen PDB path: {target_antigen_pdb_file_path}")
            if os.path.isfile(target_antigen_pdb_file_path):
                self.target_antigen_pdb_path = target_antigen_pdb_file_path
                pdb_file_basename = os.path.basename(target_antigen_pdb_file_path)
                # Infer PDB ID from filename (e.g., "1abc.pdb" -> "1abc")
                self.target_antigen_pdb_id = pdb_file_basename.split('.')[0].lower()
                logging.info(f"✅ Successfully set antigen PDB path from local file: {self.target_antigen_pdb_path} (Inferred ID: {self.target_antigen_pdb_id})")

                if target_antigen_sequence:
                    self.target_antigen_sequence = target_antigen_sequence
                    logging.info(f"Using user-provided sequence for local PDB {self.target_antigen_pdb_id}.")
                elif target_antigen_pdb_chain_id: 
                    try:
                        # If chain_id is given, assume it corresponds to this local PDB
                        pdb_id_from_chain_arg, chain_val = target_antigen_pdb_chain_id.split('_', 1)
                        if pdb_id_from_chain_arg.lower() == self.target_antigen_pdb_id: # Ensure consistency
                            self.target_antigen_chain_id = chain_val.upper()
                            logging.info(f"Attempting to fetch sequence for local PDB ID '{self.target_antigen_pdb_id}' and chain '{self.target_antigen_chain_id}' (RCSB)...")
                            seq_details = rcsb.get_pdb_chain_sequence_details(self.target_antigen_pdb_id, self.target_antigen_chain_id)
                            if seq_details and seq_details.get('sequence'):
                                self.target_antigen_sequence = seq_details['sequence']
                                logging.info(f"Fetched sequence for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id} (len: {len(self.target_antigen_sequence)}).")
                            else:
                                logging.warning(f"⚠️Could not fetch sequence for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id} from RCSB for local PDB.")
                        else:
                            logging.warning(f"⚠️ PDB ID from target_antigen_pdb_chain_id ('{pdb_id_from_chain_arg}') does not match inferred PDB ID from local file ('{self.target_antigen_pdb_id}'). Chain not used for sequence fetching.")
                    except Exception as e:
                        logging.warning(f"⚠️ Error attempting to fetch sequence for local PDB via RCSB using provided chain ID: {e}")
                else:
                    logging.info("Local PDB provided. Sequence not explicitly given or fetched.")
            else:
                logging.warning(f"⚠️ Provided target_antigen_pdb_file_path '{target_antigen_pdb_file_path}' is not a valid file. Moving to next priority.")

        # --- Priority 2: PDB ID for Download ---
        if not self.target_antigen_pdb_path and target_pdb_id:
            logging.info(f"Priority 2: 🔍 Using PDB ID '{target_pdb_id}' for download from RCSB.")
            downloaded_pdb_path = self._materialize_antigen_pdb_file(target_pdb_id.lower(), self.antigen_pdb_download_path)
            if downloaded_pdb_path:
                self.target_antigen_pdb_path = downloaded_pdb_path
                self.target_antigen_pdb_id = target_pdb_id.lower() # Store the PDB ID used for download
                logging.info(f"✅ Successfully downloaded antigen PDB: {self.target_antigen_pdb_path} for ID {self.target_antigen_pdb_id}")

                if target_antigen_sequence:
                    self.target_antigen_sequence = target_antigen_sequence
                    logging.info(f"Using user-provided sequence for downloaded PDB {self.target_antigen_pdb_id}.")
                elif target_antigen_pdb_chain_id:
                    try:
                        pdb_id_from_chain_arg, chain_val = target_antigen_pdb_chain_id.split('_', 1)
                        if pdb_id_from_chain_arg.lower() == self.target_antigen_pdb_id: # Chain must match the downloaded PDB ID
                            self.target_antigen_chain_id = chain_val.upper()
                            logging.info(f"Attempting to fetch sequence for downloaded PDB ID '{self.target_antigen_pdb_id}' and chain '{self.target_antigen_chain_id}' (RCSB)...")
                            seq_details = rcsb.get_pdb_chain_sequence_details(self.target_antigen_pdb_id, self.target_antigen_chain_id)
                            if seq_details and seq_details.get('sequence'):
                                self.target_antigen_sequence = seq_details['sequence']
                                logging.info(f"Fetched sequence (len: {len(self.target_antigen_sequence)}) for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id} via RCSB.")
                            else:
                                logging.warning(f"⚠️ Could not fetch sequence for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id} from RCSB after PDB download.")
                        else:
                            logging.warning(f"⚠️ PDB ID from target_antigen_pdb_chain_id ('{pdb_id_from_chain_arg}') does not match downloaded PDB ID ('{self.target_antigen_pdb_id}'). Chain not used for sequence fetching.")
                    except Exception as e:
                        logging.warning(f"⚠️ Error processing target_antigen_pdb_chain_id ('{target_antigen_pdb_chain_id}') with downloaded PDB: {e}")
                else:
                    logging.info(f"PDB {self.target_antigen_pdb_id} downloaded. Sequence not explicitly provided or fetched via chain ID.")
            else:
                logging.warning(f"⚠️ Failed to download PDB for ID '{target_pdb_id}'. Moving to next priority.")

        # --- Priority 3: PDB ID and Chain (Microservice Sequence -> Predict Structure) ---
        if not self.target_antigen_pdb_path and target_antigen_pdb_chain_id:
            logging.info(f"Priority 3: 🔍 Using PDB_CHAIN_ID '{target_antigen_pdb_chain_id}' for MoreMi microservice sequence retrieval and prediction.")
            try:
                # Ensure format is PDBID_CHAIN for the microservice
                pdb_id_part_ms, chain_id_part_ms = target_antigen_pdb_chain_id.split('_', 1)
                if not pdb_id_part_ms or not chain_id_part_ms:
                    raise ValueError("PDB ID or Chain part empty in target_antigen_pdb_chain_id")

                logging.info(f"Attempting to fetch sequence via MoreMi microservice for '{target_antigen_pdb_chain_id}'.")
                seq_details_microservice = rcsb.fetch_sequence_details_by_pdb_chain_id(target_antigen_pdb_chain_id)

                if seq_details_microservice and seq_details_microservice.get('sequence'):
                    self.target_antigen_sequence = seq_details_microservice['sequence']
                    # Use PDB ID and Chain from microservice response for consistency
                    self.target_antigen_pdb_id = seq_details_microservice.get('pdb_id', pdb_id_part_ms.lower())
                    self.target_antigen_chain_id = seq_details_microservice.get('chain_id', chain_id_part_ms.upper())
                    logging.info(f"🔍 Fetched sequence (len: {len(self.target_antigen_sequence)}) for {self.target_antigen_pdb_id}_{self.target_antigen_chain_id} via microservice. Now predicting structure.")
                        
                    predicted_pdb_path = self._predict_antigen_structure(self.target_antigen_sequence, self.antigen_pdb_download_path)
                    if predicted_pdb_path:
                        self.target_antigen_pdb_path = predicted_pdb_path
                        # Update PDB ID for this predicted structure to reflect its origin
                        self.target_antigen_pdb_id = f"predicted_msvc_{self.target_antigen_pdb_id}_{self.target_antigen_chain_id}"
                        # Chain ID is already set from microservice if applicable
                        logging.info(f"✅ Successfully predicted structure from microservice sequence: {self.target_antigen_pdb_path}")
                    else:
                        logging.warning(f"⚠️ Failed to predict structure for sequence from microservice ({target_antigen_pdb_chain_id}). PDB path remains unset. Sequence (len: {len(self.target_antigen_sequence)}) is available if needed for sequence-only prediction.")
                        # Keep self.target_antigen_sequence, it might be used in Priority 4
                else:
                    logging.warning(f"⚠️ Could not fetch sequence from MoreMi microservice for '{target_antigen_pdb_chain_id}'. Cannot predict structure this way. Moving to next priority.")
            except ValueError as e:
                logging.error(f"❌ Invalid format for target_antigen_pdb_chain_id '{target_antigen_pdb_chain_id}' for microservice: {e}. Moving to next priority.")
            except Exception as e:
                logging.error(f"❌ Error processing PDB_CHAIN_ID '{target_antigen_pdb_chain_id}' with microservice: {e}. Moving to next priority.")

        # --- Priority 4: Antigen Sequence Only (Predict Structure) ---
        # This will also catch sequences derived from Priority 3 if structure prediction failed there but sequence was obtained
        if not self.target_antigen_pdb_path and (target_antigen_sequence or self.target_antigen_sequence):
            # Prioritize explicitly passed target_antigen_sequence if both exist
            current_sequence_to_predict = target_antigen_sequence if target_antigen_sequence else self.target_antigen_sequence

            logging.info(f"Priority 4: 🔍 Using provided/derived antigen sequence (len: {len(current_sequence_to_predict)}) for structure prediction as no PDB path resolved yet.")
            # Ensure self.target_antigen_sequence is set with the sequence we are about to use
            if not self.target_antigen_sequence: self.target_antigen_sequence = current_sequence_to_predict
            
            predicted_pdb_path = self._predict_antigen_structure(current_sequence_to_predict, self.antigen_pdb_download_path)
            if predicted_pdb_path:
                self.target_antigen_pdb_path = predicted_pdb_path
                self.target_antigen_pdb_id = "predicted_user_sequence" # Or a more specific ID if derived from microservice earlier
                if hasattr(self, 'target_antigen_chain_id') and self.target_antigen_chain_id and "predicted_msvc" in (self.target_antigen_pdb_id or ""):
                    # If chain and msvc origin were set, keep them, but PDB ID reflects it's now a prediction
                    pass # Chain ID is already set
                else: # Pure sequence input
                    self.target_antigen_chain_id = None 
                logging.info(f"✅ Successfully predicted structure from user-provided/derived antigen sequence: {self.target_antigen_pdb_path}")
            else:
                logging.warning(f"⚠️ Failed to predict structure for the user-provided/derived antigen sequence. PDB path remains unset.")
                self.target_antigen_sequence = None # If prediction fails, sequence alone without structure might not be useful for binding.
        
        # --- Final Outcome --- 
        if self.target_antigen_pdb_path:
            log_msg = (
                f"✅ Antigen context successfully set with PDB_Path='{self.target_antigen_pdb_path}'. "
                f"PDB_ID='{self.target_antigen_pdb_id or 'N/A'}', "
                f"Chain_ID='{self.target_antigen_chain_id or 'N/A'}'. "
            )
            if self.target_antigen_sequence:
                log_msg += f"Sequence (len {len(self.target_antigen_sequence)}) also obtained."
            else:
                log_msg += "Sequence not obtained or not applicable for this PDB source."
            logging.info(log_msg)
            return True
        else:
            logging.warning(
                f"Antigen context setup FAILED to obtain a PDB path. "
                "Binding affinity and other antigen-dependent metrics may not run or may be inaccurate."
            )
            # Ensure clean state if completely failed to get PDB path
            self.target_antigen_sequence = None
            self.target_antigen_pdb_path = None
            self.target_antigen_pdb_id = None
            self.target_antigen_chain_id = None
            return False

    def process_protein(self, sequence: str) -> ProcessingResult:
        """Process a single protein sequence and calculate all metrics"""
        protein_metrics = None
        try:
            # Basic sequence validation
            if not sequence or len(sequence) < 10: # Basic check
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
            # Simplified molecular formula - real one is complex. This is a placeholder.
            molecular_formula = f"C{analyzed_seq.count_amino_acids().get('C',0)}H{analyzed_seq.count_amino_acids().get('H',0)}N{analyzed_seq.count_amino_acids().get('N',0)}O{analyzed_seq.count_amino_acids().get('O',0)}S{analyzed_seq.count_amino_acids().get('S',0)}"
            
            # Collect all metrics
            metrics_data = {} # Using a temporary dict to store results before filling ProteinMetrics
            
            # Initialize metrics container early so we can add warnings
            # Antigen details will be filled based on self.target_antigen_sequence etc. if context was set
            protein_metrics = ProteinMetrics(
                sequence=sequence,
                antigen=self.target_antigen_sequence or "Antigen context not set", # Default if not set
                antigen_id=self.target_antigen_pdb_id or "Unknown",
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
                warnings=[]
            )
            
            # Set antigen_pdb_chain_id in ProteinMetrics if available from context
            if self.target_antigen_pdb_id and self.target_antigen_chain_id:
                protein_metrics.antigen_pdb_chain_id = f"{self.target_antigen_pdb_id}_{self.target_antigen_chain_id}"
            elif self.target_antigen_pdb_id == "predicted_antigen":
                protein_metrics.antigen_pdb_chain_id = "predicted_structure"
            
            # Print tree-style validation process header
            print(f"🧪 Processing protein sequence: {sequence[:20]}...")
            print(f"├── Calculating basic properties...")
            
            # BLAST analysis
            if self.metrics_to_run is None or MetricCategory.BLAST in self.metrics_to_run:
                logging.info(f"├── 🧬 Running BLAST analysis....")
                print("├── 🧬 Running BLAST analysis...")
                try:
                    blast_result = perform_blast(sequence)
                    protein_metrics.blast = blast_result.get('blast_result', {"error": "No BLAST result"})
                    logging.info(f"│   └── ✓ BLAST complete")
                    print("│   └── ✓ BLAST complete")
                except Exception as e:
                    error_msg = f"BLAST analysis failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.BLAST.value, error_msg)
                    protein_metrics.blast = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🧬 BLAST analysis {skip_msg.lower()}")
                print(f"├── 🧬 BLAST analysis {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.BLAST.value, skip_msg)
                protein_metrics.blast = {"status": skip_msg}
            
            # ProtParam analysis
            if self.metrics_to_run is None or MetricCategory.PROTPARAM in self.metrics_to_run:
                logging.info(f"├── 🔍 Calculating ProtParam properties...")
                print("├── 🔍 Analyzing ProtParam properties...")
                try:
                    protparam_result = analyze_with_protparam(sequence)
                    protein_metrics.protparam = protparam_result if protparam_result else {"error": "No ProtParam result"}
                    logging.info(f"│   └── ✓ ProtParam complete")
                    print("│   └── ✓ ProtParam complete")
                except Exception as e:
                    error_msg = f"ProtParam analysis failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.PROTPARAM.value, error_msg)
                    protein_metrics.protparam = {"error": error_msg}
            # This was an error in original logic, corrected the else placement
            else: # This else belongs to the "if self.metrics_to_run is None or MetricCategory.PROTPARAM in self.metrics_to_run"
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🔍 ProtParam analysis {skip_msg.lower()}")
                print(f"├── 🔍 ProtParam analysis {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.PROTPARAM.value, skip_msg)
                protein_metrics.protparam = {"status": skip_msg}
            
            # Immunogenicity prediction
            if self.metrics_to_run is None or MetricCategory.IMMUNOGENICITY in self.metrics_to_run:
                logging.info(f"├── 🦠 Assessing immunogenicity...")
                print("├── 🦠 Assessing immunogenicity...")
                try:
                    immunogenicity_result = predict_immunogenicity(sequence)
                    protein_metrics.immunogenicity = immunogenicity_result if immunogenicity_result else {"error": "No immunogenicity result"}
                    logging.info(f"│   └── ✓ Immunogenicity assessment complete")
                    print("│   └── ✓ Immunogenicity assessment complete")
                except Exception as e:
                    error_msg = f"Immunogenicity prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.IMMUNOGENICITY.value, error_msg)
                    protein_metrics.immunogenicity = {"error": error_msg}
            # This was an error in original logic, corrected the else placement
            else: # This else belongs to the "if self.metrics_to_run is None or MetricCategory.IMMUNOGENICITY in self.metrics_to_run"
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🦠 Immunogenicity assessment {skip_msg.lower()}")
                print(f"├── 🦠 Immunogenicity assessment {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.IMMUNOGENICITY.value, skip_msg)
                protein_metrics.immunogenicity = {"status": skip_msg}
            
            # Stability prediction
            if self.metrics_to_run is None or MetricCategory.STABILITY in self.metrics_to_run:
                logging.info(f"├── 🔥 Evaluating stability...")
                print("├── 🔥 Evaluating stability...")
                try:
                    stability_result = predict_stability(sequence)
                    protein_metrics.stability = stability_result.get('stability_result', {"error": "No stability result"})
                    logging.info(f"│   └── ✓ Stability evaluation complete")
                    print("│   └── ✓ Stability evaluation complete")
                except Exception as e:
                    error_msg = f"Stability prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.STABILITY.value, error_msg)
                    protein_metrics.stability = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🔥 Stability evaluation {skip_msg.lower()}")
                print(f"├── 🔥 Stability evaluation {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.STABILITY.value, skip_msg)
                protein_metrics.stability = {"status": skip_msg}
            
            # Aggregation prediction
            if self.metrics_to_run is None or MetricCategory.AGGREGATION in self.metrics_to_run:
                logging.info(f"├── 🧱 Predicting aggregation propensity...")
                print("├── 🧱 Predicting aggregation propensity...")
                try:
                    aggregation_result = predict_aggregation(sequence)
                    protein_metrics.aggregation = aggregation_result.get('aggregation_result', {"error": "No aggregation result"})
                    logging.info(f"│   └── ✓ Aggregation prediction complete")
                    print("│   └── ✓ Aggregation prediction complete")
                except Exception as e:
                    error_msg = f"Aggregation prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.AGGREGATION.value, error_msg)
                    protein_metrics.aggregation = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🧱 Aggregation prediction {skip_msg.lower()}")
                print(f"├── 🧱 Aggregation prediction {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.AGGREGATION.value, skip_msg)
                protein_metrics.aggregation = {"status": skip_msg}
            
            # Glycosylation prediction
            if self.metrics_to_run is None or MetricCategory.GLYCOSYLATION in self.metrics_to_run:
                logging.info(f"├── 🍭  Identifying glycosylation sites...")
                print("├── 🍭 Identifying glycosylation sites...")
                try:
                    glycosylation_result = predict_glycosylation(sequence)
                    protein_metrics.glycosylation = glycosylation_result if glycosylation_result else {"error": "No glycosylation result"}
                    logging.info(f"│   └── ✓ Glycosylation sites identified")
                    print("│   └── ✓ Glycosylation sites identified")
                except Exception as e:
                    error_msg = f"Glycosylation prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.GLYCOSYLATION.value, error_msg)
                    protein_metrics.glycosylation = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🍭 Glycosylation sites identification {skip_msg.lower()}")
                print(f"├── 🍭 Glycosylation sites identification {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.GLYCOSYLATION.value, skip_msg)
                protein_metrics.glycosylation = {"status": skip_msg}
            
            # Structure prediction (for the antibody itself)
            antibody_pdb_path = None # Will store path to the antibody's PDB
            if self.metrics_to_run is None or MetricCategory.STRUCTURE in self.metrics_to_run:
                logging.info(f"├── 🧩 Generating structural model (antibody)...")
                print("├── 🧩 Generating structural model (antibody)...")
                try:
                    # Use self.pdb_files_path for antibody structures
                    structure_result = predict_structure(sequence, output_pdb_filename_prefix=f"antibody_{molecular_formula}", output_directory=self.pdb_files_path)
                    protein_metrics.structure = structure_result if structure_result else {"error": "No structure result"}
                    if structure_result and structure_result.get("pdb_file_path"):
                        antibody_pdb_path = structure_result["pdb_file_path"]
                        logging.info(f"│   └── ✓ Antibody structural model generated: {antibody_pdb_path}")
                        print("│   └── ✓ Antibody structural model generated")
                    else:
                        logging.warning(f"│   └── ⚠️ Antibody structural model generation failed or no PDB path returned.")
                        print("│   └── ⚠️ Antibody structural model generation failed")
                except Exception as e:
                    error_msg = f"Antibody structure prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.STRUCTURE.value, error_msg)
                    protein_metrics.structure = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🧩 Antibody structural model generation {skip_msg.lower()}")
                print(f"├── 🧩 Antibody structural model generation {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.STRUCTURE.value, skip_msg)
                protein_metrics.structure = {"status": skip_msg}
            
            # Binding affinity prediction
            if self.metrics_to_run is None or MetricCategory.BINDING_AFFINITY in self.metrics_to_run:
                logging.info(f"├── 🔗 Calculating binding affinity...")
                print("├── 🔗 Calculating binding affinity...")
                # Check if antigen context and antibody structure are available
                if (self.target_antigen_pdb_path and antibody_pdb_path) or (self.target_antigen_sequence and antibody_pdb_path):
                    try:
                        # Note: predict_binding_affinity expects PDB paths.
                        # self.target_antigen_pdb_path is set by set_antigen_context
                        # antibody_pdb_path is from the structure prediction of the antibody
                        binding_result = predict_binding_affinity(
                            receptor_pdb_path=self.target_antigen_pdb_path, # Antigen is receptor
                            ligand_pdb_path=antibody_pdb_path      # Antibody is ligand
                        )
                        protein_metrics.binding_affinity = binding_result if binding_result else {"error": "No binding affinity result"}
                        antigen_source_log = f"antigen PDB {self.target_antigen_pdb_id or os.path.basename(self.target_antigen_pdb_path)}"
                        logging.info(f"│   └── ✓ Binding affinity calculated against {antigen_source_log}")
                        print(f"│   └── ✓ Binding affinity calculated against {antigen_source_log}")
                    except ValueError as e:
                        error_msg = f"Binding affinity calculation failed: {str(e)}"
                        logging.warning(f"│   └── ⚠️ {error_msg}")
                        print(f"│   └── ⚠️ {error_msg}")
                        self._add_warning(protein_metrics, MetricCategory.BINDING_AFFINITY.value, error_msg)
                        protein_metrics.binding_affinity = {"error": error_msg}
                    except Exception as e:
                        error_msg = f"Binding affinity calculation failed: {str(e)}"
                        logging.warning(f"│   └── ⚠️ {error_msg}")
                        print(f"│   └── ⚠️ {error_msg}")
                        self._add_warning(protein_metrics, MetricCategory.BINDING_AFFINITY.value, error_msg)
                        protein_metrics.binding_affinity = {"error": error_msg}
                else:
                    missing_items = []
                    skip_msg = "Skipped by user configuration."
                    if not self.target_antigen_pdb_path: missing_items.append("antigen PDB structure (context not set or failed)")
                    if not self.target_antigen_sequence: 
                        missing_items.append("antigen sequence (context not set or failed)")
                    if not antibody_pdb_path:
                        if protein_metrics.structure.get("status") == skip_msg:
                            missing_items.append("antibody PDB structure (skipped by user)")
                        else:
                             missing_items.append("antibody PDB structure (prediction failed or not available)")
                    
                    missing_str = " and ".join(missing_items) if missing_items else "required data (antigen PDB/sequence or antibody PDB)"
                    error_msg = f"Binding affinity calculation skipped: Missing {missing_str}"
                    protein_metrics.binding_affinity = {"error": error_msg}
                    self._add_warning(protein_metrics, MetricCategory.BINDING_AFFINITY.value, error_msg)
                    logging.warning(f"│   └── ⚠️ Skipping binding affinity. Antibody PDB: {antibody_pdb_path}, Antigen PDB: {self.target_antigen_pdb_path}, Antigen Seq: {bool(self.target_antigen_sequence)}")
                    print(f"│   └── ⚠️ Skipping binding affinity: missing {missing_str}")
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🔗 Binding affinity calculation {skip_msg.lower()}")
                print(f"├── 🔗 Binding affinity calculation {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.BINDING_AFFINITY.value, skip_msg)
                protein_metrics.binding_affinity = {"status": skip_msg}

            # Epitope prediction
            if self.metrics_to_run is None or MetricCategory.EPITOPE in self.metrics_to_run:
                logging.info(f"├── 🎯 Predicting epitopes...")
                print("├── 🎯 Predicting epitope regions...")
                try:
                    epitope_result = predict_bcell_epitopes(sequence)
                    protein_metrics.epitope = epitope_result if epitope_result else {"error": "No epitope result"}
                    logging.info(f"│   └── ✓ Epitope prediction complete")
                    print("│   └── ✓ Epitope prediction complete")
                except Exception as e:
                    error_msg = f"Epitope prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.EPITOPE.value, error_msg)
                    protein_metrics.epitope = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🎯 Epitope prediction {skip_msg.lower()}")
                print(f"├── 🎯 Epitope prediction {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.EPITOPE.value, skip_msg)
                protein_metrics.epitope = {"status": skip_msg}
            
            # Conservancy prediction
            if self.metrics_to_run is None or MetricCategory.CONSERVANCY in self.metrics_to_run:
                logging.info(f"├── 🌐 Analyzing sequence conservancy...")
                print("├── 🌐 Analyzing sequence conservancy...")
                try:
                    epitopes_for_conservancy = []
                    if isinstance(protein_metrics.epitope, dict) and protein_metrics.epitope.get('epitope_sequences_list'):
                        epitopes_for_conservancy = protein_metrics.epitope['epitope_sequences_list']
                    elif isinstance(protein_metrics.epitope, dict) and (protein_metrics.epitope.get("status") == "Skipped by user configuration." or "error" in protein_metrics.epitope) :
                        warn_skip_msg = "Skipped: dependent epitope prediction was skipped or failed."
                        self._add_warning(protein_metrics, MetricCategory.CONSERVANCY.value, warn_skip_msg)
                        raise ValueError(warn_skip_msg) # This will be caught and logged as error for conservancy

                    conservancy_result = predict_conservancy(sequence, epitopes_for_conservancy, output_dir = None, save_csv = False)
                    protein_metrics.conservancy = conservancy_result if conservancy_result else {"error": "No conservancy result"}
                    logging.info(f"│   └── ✓ Conservancy analysis complete")
                    print("│   └── ✓ Conservancy analysis complete")
                except Exception as e:
                    error_msg = f"Conservancy prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.CONSERVANCY.value, error_msg)
                    protein_metrics.conservancy = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🌐 Conservancy analysis {skip_msg.lower()}")
                print(f"├── 🌐 Conservancy analysis {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.CONSERVANCY.value, skip_msg)
                protein_metrics.conservancy = {"status": skip_msg}
            
            # Developability prediction
            if self.metrics_to_run is None or MetricCategory.DEVELOPABILITY in self.metrics_to_run:
                logging.info(f"├── 🔧 Assessing developability...")
                print("├── 🔧 Assessing developability...")
                try:
                    developability_result = predict_developability(sequence)
                    protein_metrics.developability = developability_result if developability_result else {"error": "No developability result"}
                    logging.info(f"│   └── ✓ Developability assessment complete")
                    print("│   └── ✓ Developability assessment complete")
                except Exception as e:
                    error_msg = f"Developability prediction failed: {str(e)}"
                    logging.warning(f"│   └── ⚠️ {error_msg}")
                    print(f"│   └── ⚠️ {error_msg}")
                    self._add_warning(protein_metrics, MetricCategory.DEVELOPABILITY.value, error_msg)
                    protein_metrics.developability = {"error": error_msg}
            else:
                skip_msg = "Skipped by user configuration."
                logging.info(f"├── 🔧 Developability assessment {skip_msg.lower()}")
                print(f"├── 🔧 Developability assessment {skip_msg}")
                self._add_warning(protein_metrics, MetricCategory.DEVELOPABILITY.value, skip_msg)
                protein_metrics.developability = {"status": skip_msg}
            
            # Check how many metrics categories have errors
            current_metrics_dict = protein_metrics.to_dict()['metrics']
            error_categories = [cat for cat, val in current_metrics_dict.items() if isinstance(val, dict) and "error" in val]
            
            if error_categories:
                self._add_warning(protein_metrics, "General", 
                                f"Failed to calculate results for {len(error_categories)} metric categories: {', '.join(error_categories)}")
                if len(error_categories) > len(MetricCategory) * 0.3:
                    self._add_warning(protein_metrics, "Critical", 
                                    f"More than 30% of metric categories have errors or were skipped. Results may be unreliable.")
            
            logging.info(f"└── 📊 Collection of metrics complete...")
            print("└── 📊 Collection of metrics complete")
            
            if protein_metrics.warnings:
                warning_count = len(protein_metrics.warnings)
                logging.info(f"⚠️ Protein has {warning_count} warnings/issues")
                print(f"⚠️ Protein has {warning_count} warnings/issues")
            
            logging.info(f"\n\n✅ Successfully processed protein {sequence[:20]}...\n\n")
            print(f"\n\n✅ Successfully processed protein {sequence[:20]}...\n\n")
            
            return ProcessingResult(
                sequence=sequence,
                metrics=protein_metrics,
                success=True
            )
            
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            error_msg = f"Critical error during processing of sequence {sequence[:20]}...: {str(e)}"
            logging.error(f"❌ {error_msg}{tb}")
            print(f"❌ {error_msg}")
            
            # Ensure protein_metrics is at least the initialized version if an early error occurred
            if protein_metrics is None and 'sequence' in locals() and 'molecular_weight' in locals() and 'molecular_formula' in locals():
                 protein_metrics = ProteinMetrics(
                    sequence=sequence, antigen="Error before antigen set", antigen_id="Error",
                    molecular_weight=molecular_weight, molecular_formula=molecular_formula,
                    blast={}, protparam={}, immunogenicity={}, stability={}, aggregation={},
                    glycosylation={}, binding_affinity={}, structure={}, epitope={},
                    developability={}, conservancy={}, warnings=[f"Critical Processing Error: {str(e)}"]
                )
            elif protein_metrics is not None:
                 self._add_warning(protein_metrics, "Critical Processing Error", str(e))


            return ProcessingResult(
                sequence=sequence,
                metrics=protein_metrics, 
                success=False,
                error=error_msg
            )

    def process_proteins(self, input_file: str, output_dir: Optional[str] = None) -> List[ProcessingResult]:
        """Process multiple protein sequences from a file (text or FASTA) and save results to output directory."""
        # output_dir for results of this batch run, not related to self.pdb_files_path (antibody) or self.antigen_pdb_download_path
        effective_output_dir = output_dir or os.path.join(os.getcwd(), "protein_validator_v2_results", datetime.now().strftime("%Y%m%d_%H%M%S"))
        try:
            os.makedirs(effective_output_dir, exist_ok=True)
            
            sequences = []
            current_sequence_lines = []
            with open(input_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("//"):
                        continue
                    if line.startswith(">"):
                        if current_sequence_lines:
                            sequences.append("".join(current_sequence_lines))
                            current_sequence_lines = []
                    else:
                        current_sequence_lines.append(line)
            if current_sequence_lines:
                sequences.append("".join(current_sequence_lines))

            if not sequences:
                logging.warning(f"No valid sequences found in {input_file}")
                return []
                
            logging.info(f"Found {len(sequences)} protein sequences in {input_file}")
            
            all_results: List[ProcessingResult] = []
            success_count = 0
            error_count = 0
            
            logging.info(f"🌲 Starting batch processing of {len(sequences)} proteins...")
            print(f"🌲 Processing {len(sequences)} proteins from {input_file} into {effective_output_dir}")
            
            for i, sequence in enumerate(sequences):
                try:
                    logging.info(f"Processing protein {i+1}/{len(sequences)}:")
                    # Ensure antigen context is set if needed by any metric in metrics_to_run
                    # For now, we assume if BINDING_AFFINITY is to be run, set_antigen_context MUST have been called by user.
                    # A more robust check could be added here or at the start of process_proteins
                    # to ensure antigen context is ready if MetricCategory.BINDING_AFFINITY is in self.metrics_to_run.
                    
                    result = self.process_protein(sequence)
                    all_results.append(result)
                    
                    if result.success:
                        success_count += 1
                        logging.info(f"✅ Successfully processed protein {i+1}")
                        print(f"│  ├── {i+1}/{len(sequences)}: ✅ Success")
                    else:
                        error_count += 1
                        logging.warning(f"❌ Failed to process protein {i+1}: {result.error}")
                        print(f"│  ├── {i+1}/{len(sequences)}: ❌ Error - {result.error}")
                
                except Exception as e: # Catch unexpected errors in the loop itself
                    error_count += 1
                    error_msg = f"Unexpected error in batch loop for protein {i+1}: {str(e)}"
                    logging.error(f"❌ {error_msg}")
                    print(f"│  ├── {i+1}/{len(sequences)}: ❌ Exception - {error_msg}")
                    all_results.append(ProcessingResult(
                        sequence=sequence, metrics=None, success=False, error=error_msg
                    ))
            
            print("│")
            print(f"└── 📊 Summary: {success_count} successful, {error_count} failed")
            
            for i, res in enumerate(all_results):
                status_symbol = "✅" if res.success else "❌"
                detail = f"Sequence: {res.sequence[:20]}..." if res.success and res.metrics else f"Error: {res.error}"
                print(f"    ├── {i+1}: {status_symbol} {detail}")
            
            logging.info(f"Completed processing {len(sequences)} proteins")
            logging.info(f"✅ Success: {success_count}, ❌ Failed: {error_count}")
            
            self.save_processing_results(all_results, effective_output_dir)
            logging.info(f"💾 Batch results saved to {effective_output_dir}")
            
            return all_results
            
        except Exception as e:
            logging.error(f"❌ Error in batch processing (process_proteins): {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
            return locals().get('all_results', []) # Return any results collected so far, or empty list
    
    def save_individual_metric_result(self, metrics: ProteinMetrics, output_dir: str):
        """Save individual protein metrics to a JSON file in a 'json' subdirectory."""
        metrics_dict = metrics.to_dict()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f") # Added microseconds for uniqueness
        
        json_subdir = os.path.join(output_dir, "json_metrics")
        os.makedirs(json_subdir, exist_ok=True)

        # Use a more unique filename, perhaps based on sequence hash or part of sequence if safe
        seq_part = metrics.sequence[:10].replace("/","_") # Basic sanitization
        output_file = os.path.join(json_subdir, f"{metrics_dict.get('molecular_formula','protein')}_{seq_part}_{timestamp}.json")
        
        try:
            with open(output_file, 'w') as f:
                import json
                json.dump(metrics_dict, f, indent=2)
            logging.info(f"Saved individual metrics to {output_file}")
        except Exception as e:
            logging.error(f"Error saving individual metrics to {output_file}: {e}")

    def save_batch_summary(self, results: List[ProcessingResult], output_dir: str):
        """Save summary of all protein processing results for the batch."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        summary_file = os.path.join(output_dir, f"batch_validation_summary_{timestamp}.txt")
        
        successful_results = [r for r in results if r.success and r.metrics is not None]
        failed_results = [r for r in results if not r.success or r.metrics is None]
        
        with open(summary_file, 'w') as f:
            f.write("Protein Batch Validation Summary")
            f.write("==============================")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            f.write(f"Output Directory: {os.path.abspath(output_dir)}")
            
            f.write("Overall Statistics:")
            f.write("-----------------")
            f.write(f"Total proteins attempted: {len(results)}")
            f.write(f"Successfully processed: {len(successful_results)}")
            f.write(f"Failed or no metrics: {len(failed_results)}")
            
            if successful_results:
                f.write("Successfully Processed Proteins (Summary - first 10):")
                f.write("--------------------------------------------------")
                for i, result in enumerate(successful_results[:10]):
                    f.write(f"{i+1}. Sequence: {result.sequence[:30]}...")
                    if result.metrics:
                        f.write(f"   Formula: {result.metrics.molecular_formula}")
                        f.write(f"   Weight: {result.metrics.molecular_weight:.2f}")
                        f.write(f"   Warnings: {len(result.metrics.warnings)}")
                    f.write("-" * 50 + "")
                if len(successful_results) > 10:
                    f.write(f"... and {len(successful_results) - 10} more successfully processed proteins.")
        
            if failed_results:
                f.write("Failed or Incomplete Proteins:")
                f.write("----------------------------")
                for result in failed_results:
                    f.write(f"Sequence: {result.sequence[:30]}...")
                    f.write(f"Error: {result.error or 'Metrics not generated'}")
                    f.write("-" * 50 + "")
        logging.info(f"Batch summary saved to {summary_file}")

    def _read_sequences_from_input(self, input_source: Union[str, List[str], Path]) -> List[str]:
        """
        Reads protein sequences from various input types.

        Args:
            input_source: Can be a single sequence string, a list of sequence strings,
                          or a file path (string or Path object) to a file containing sequences.
                          The file can be plain text (one sequence per line) or FASTA format.

        Returns:
            A list of protein sequence strings.
        """
        sequences: List[str] = []
        input_path: Optional[Path] = None

        if isinstance(input_source, list):
            if all(isinstance(item, str) for item in input_source):
                sequences = [item for item in input_source if item]  # Filter out empty strings
            else:
                logging.warning("Input list contains non-string items. Only non-empty string items will be treated as sequences.")
                sequences = [str(item) for item in input_source if isinstance(item, str) and item]
        elif isinstance(input_source, str):
            # Try to interpret as a file path first
            potential_path = Path(input_source)
            if potential_path.is_file():
                input_path = potential_path
            elif input_source:  # Not a file, treat as a single sequence string
                logging.info("Input string is not a file, treating as a single sequence.")
                sequences = [input_source]
            else: # Empty string
                logging.warning("Input string is empty and not a file path.")
        elif isinstance(input_source, Path):
            if input_source.is_file():
                input_path = input_source
            else:
                logging.warning(f"Input Path object {input_source} is not a valid file.")
        else:
            logging.warning(f"Unsupported input type for sequences: {type(input_source)}. Expected string, list of strings, or file path.")

        if input_path:
            logging.info(f"📖 Reading sequences from file: {input_path}")
            try:
                with open(input_path, 'r') as f:
                    lines_from_file = [line.strip() for line in f if line.strip()] # Read all non-empty stripped lines
                
                if not lines_from_file:
                    logging.warning(f"No non-empty lines found in file {input_path}.")
                    return []

                # Determine if the file is likely FASTA by checking for '>' in any relevant line
                is_fasta = any(line.startswith(">") for line in lines_from_file)

                if is_fasta:
                    logging.debug(f"Detected FASTA format for {input_path}")
                    current_sequence_lines: List[str] = []
                    for line in lines_from_file:
                        if line.startswith(">"):
                            if current_sequence_lines:  # Save previous sequence if any
                                sequences.append("".join(current_sequence_lines))
                                current_sequence_lines = []
                            # FASTA header itself is ignored for sequence collection
                        else:  # Sequence line, remove internal spaces/newlines just in case
                            current_sequence_lines.append(line.replace(" ", "").replace("\\n", ""))
                    if current_sequence_lines:  # Add the last sequence in the file
                        sequences.append("".join(current_sequence_lines))
                else:  # Plain text, one sequence per line (after initial strip)
                    logging.debug(f"Detected plain text format (one sequence per line) for {input_path}")
                    # Each line that is not a comment (already stripped) and not empty is a sequence
                    sequences.extend([line for line in lines_from_file if not line.startswith("#") and not line.startswith("//")])
                
                if not sequences:
                    logging.warning(f"No valid sequences extracted from file {input_path} (format detected: {'FASTA' if is_fasta else 'plain text'}).")
            except Exception as e:
                logging.error(f"Error reading sequence file {input_path}: {e}", exc_info=True)
                return [] # Return empty list on error

        # Final validation of sequences (e.g., remove very short or obviously invalid characters if needed)
        valid_sequences = [seq for seq in sequences if seq and len(seq) >= 3] # Example: ensure sequences are not empty and have min length
        if len(valid_sequences) != len(sequences):
            logging.debug(f"Removed {len(sequences) - len(valid_sequences)} empty or very short sequences from input.")
        
        logging.info(f"Successfully read {len(valid_sequences)} sequences from input_source.")
        return valid_sequences

    def _append_to_realtime_csv(self, result: ProcessingResult, csv_path: str):
        """Appends a single ProcessingResult to a CSV file, creating/writing headers if needed."""
        record_to_append = {}
        if result.success and result.metrics:
            record_to_append = result.metrics.to_dict() # This flattens metrics inside
            record_to_append['success'] = True
            record_to_append['error'] = None
             # Flatten nested metrics
            if 'metrics' in record_to_append and isinstance(record_to_append['metrics'], dict):
                for cat_key, cat_value in record_to_append['metrics'].items():
                    if isinstance(cat_value, dict):
                        for sub_key, sub_value in cat_value.items():
                            record_to_append[f'{cat_key}_{sub_key}'] = sub_value
                    else: # if metric category is not a dict (e.g. simple value or error string)
                        record_to_append[cat_key] = cat_value
                del record_to_append['metrics'] # remove original nested metrics

        else: # Failed or no metrics
            record_to_append = {
                'sequence': result.sequence,
                'success': False,
                'error': result.error or "No metrics generated",
                'antigen': "N/A", 'antigen_id': "N/A", 'antigen_pdb_chain_id': "N/A",
                'molecular_weight': None, 'molecular_formula': "N/A", 'warnings': [],
            }
        
        # Ensure warnings are stringified
        if 'warnings' in record_to_append and isinstance(record_to_append['warnings'], list):
            record_to_append['warnings'] = "; ".join(record_to_append['warnings'])

        # Convert all list/dict values in record to string for CSV compatibility
        for key, value in record_to_append.items():
            if isinstance(value, (list, dict)):
                record_to_append[key] = str(value)
        
        df_record = pd.DataFrame([record_to_append])
        
        file_exists = os.path.isfile(csv_path)
        try:
            df_record.to_csv(csv_path, mode='a', header=not file_exists, index=False)
        except Exception as e:
            logging.error(f"Error appending to realtime CSV {csv_path}: {e}")

    def get_successful_metrics(self, results: List[ProcessingResult]) -> List[ProteinMetrics]:
        """Extract only successful metrics from processing results"""
        # This method seems more general and could be a static method or utility if not using instance state
        successful = []
        failed_info = []
        
        for i, result in enumerate(results):
            if result.success and result.metrics is not None:
                successful.append(result.metrics)
            else:
                error_msg = result.error or "Unknown error or metrics not generated"
                seq_preview = result.sequence[:30] + "..." if result.sequence and len(result.sequence) > 30 else result.sequence
                failed_info.append((i + 1, seq_preview, error_msg))
        
        if failed_info:
            logging.warning(f"❌ {len(failed_info)} proteins had issues or failed processing in the batch:")
            for i, seq, error in failed_info:
                logging.warning(f"  - Protein Entry #{i} (seq: {seq}): {error}")
        
        logging.info(f"✅ Successfully extracted metrics for {len(successful)} out of {len(results)} attempted proteins in the batch.")
        return successful

    def validate_protein_list(self, 
                              input_source: Union[str, List[str], Path], 
                              output_csv_path: Optional[str] = None,
                              realtime_csv_backup_path: Optional[str] = None, # New argument
                              realtime_json_backup_path: Optional[str] = None # New argument
                              ) -> List[ProteinMetrics]:
        """
        Validate a list of protein sequences, a single protein sequence string, or sequences from a file.
        Optionally saves all (successful and failed) processed results to a CSV file.
        Also, incrementally saves all attempts to realtime_csv_backup_path and successful ProteinMetrics to realtime_json_backup_path.
        Antigen context for binding affinity must be set via `set_antigen_context` before calling this
        if binding affinity is among the metrics to be run.
        
        Args:
            input_source: A single protein sequence string, a list of protein sequences,
                          or a path (str or Path) to a file containing sequences.
            output_csv_path (Optional[str]): If provided, path to save the final validation results CSV.
            realtime_csv_backup_path (Optional[str]): Path to save all validation attempts incrementally (CSV).
            realtime_json_backup_path (Optional[str]): Path to save successful ProteinMetrics incrementally (JSON).
            
        Returns:
            List of ProteinMetrics objects for successfully validated proteins.
        """
        sequences_to_process = self._read_sequences_from_input(input_source)

        if not sequences_to_process:
            logging.warning("No sequences to validate from the provided input source.")
            if output_csv_path: 
                logging.info(f"Writing empty CSV to {output_csv_path} as no sequences were processed.")
                pd.DataFrame([]).to_csv(output_csv_path, index=False)
            if realtime_csv_backup_path and not os.path.exists(realtime_csv_backup_path): # Touch empty CSV
                 pd.DataFrame([]).to_csv(realtime_csv_backup_path, index=False)
            if realtime_json_backup_path and not os.path.exists(realtime_json_backup_path): # Touch empty JSON list
                with open(realtime_json_backup_path, 'w') as f:
                    json.dump([], f)
            return []

        binding_affinity_to_run = self.metrics_to_run is None or MetricCategory.BINDING_AFFINITY in self.metrics_to_run
        if binding_affinity_to_run and (not self.target_antigen_pdb_path): # Simplified check: PDB path is key
            logging.warning("⚠️ Binding affinity is set to run, but antigen PDB path is not set via set_antigen_context(). Binding affinity will likely fail or be skipped.")

        print(f"\nStarting validation of {len(sequences_to_process)} proteins (from input_source)...")
        
        all_processing_results: List[ProcessingResult] = []
        current_successful_metrics_list: List[ProteinMetrics] = []

        # Ensure backup directories exist
        if realtime_csv_backup_path:
            os.makedirs(os.path.dirname(realtime_csv_backup_path), exist_ok=True)
        if realtime_json_backup_path:
            os.makedirs(os.path.dirname(realtime_json_backup_path), exist_ok=True)
        
        for idx, sequence in enumerate(sequences_to_process, 1):
            print(f"\nValidating protein {idx}/{len(sequences_to_process)}: {sequence[:20]}...")
            result = self.process_protein(sequence)
            all_processing_results.append(result)

            # Real-time CSV backup of current result (attempt)
            if realtime_csv_backup_path:
                try:
                    self._append_to_realtime_csv(result, realtime_csv_backup_path)
                    logging.info(f"Appended to realtime CSV backup: {realtime_csv_backup_path} for protein {idx}")
                except Exception as e_csv:
                    logging.error(f"Failed to append to realtime CSV backup {realtime_csv_backup_path} for protein {idx}: {e_csv}")

            # Real-time JSON backup of successful metrics list
            if result.success and result.metrics:
                current_successful_metrics_list.append(result.metrics)
                if realtime_json_backup_path:
                    try:
                        # Convert list of ProteinMetrics objects to list of dicts
                        dict_list_to_save = [pm.to_dict() for pm in current_successful_metrics_list]
                        with open(realtime_json_backup_path, 'w') as f_json:
                            json.dump(dict_list_to_save, f_json, indent=2)
                        logging.info(f"Updated realtime JSON backup: {realtime_json_backup_path} with {len(current_successful_metrics_list)} successful proteins.")
                    except Exception as e_json:
                        logging.error(f"Failed to update realtime JSON backup {realtime_json_backup_path}: {e_json}")
            
        # After loop, successful_metrics can be derived from current_successful_metrics_list
        successful_metrics = current_successful_metrics_list 
        failed_count = len(sequences_to_process) - len(successful_metrics)

        print(f"\nValidation via validate_protein_list complete:")
        print(f"├── Total proteins attempted: {len(sequences_to_process)}")
        print(f"├── Successfully validated (metrics collected): {len(successful_metrics)}")
        print(f"└── Failed validation or no metrics: {failed_count}\n")

        if output_csv_path:
            logging.info(f"Saving validation results to CSV: {output_csv_path}")
            # save_metrics_to_csv expects List[ProcessingResult] which includes errors
            self.to_csv(all_processing_results, output_csv_path) # This is the final, complete CSV
        
        return successful_metrics

    def to_dict_list(self, results: List[ProcessingResult]) -> List[Dict[str, Any]]:
        """Converts ProteinMetrics from successful ProcessingResult objects into a list of dictionaries."""
        dict_list = []
        for result in results:
            if result.success and result.metrics:
                try:
                    dict_list.append(result.metrics.to_dict())
                except Exception as e:
                    logging.error(f"Error converting metrics to dict for sequence {result.sequence[:20]}...: {e}")
        return dict_list

    def to_csv(self, results: List[ProcessingResult], output_csv_path: str):
        """Saves ProteinMetrics from successful ProcessingResult objects to a CSV file.
        If self.metrics_to_run is specified, only columns related to those metrics
        (and essential non-metric columns) will be included in the CSV.
        """
        if not results:
            logging.warning("No results provided to save_metrics_to_csv.")
            print("Warning: No results to save to CSV.")
            return

        # Convert successful results to a list of dictionaries
        # Note: to_dict_list only includes successful metrics. We need ALL results for the CSV.
        all_records = []
        for result in results:
            if result.success and result.metrics:
                record = result.metrics.to_dict()
                record['success'] = True
                record['error'] = None
            else: # Failed or no metrics
                record = {
                    'sequence': result.sequence,
                    'success': False,
                    'error': result.error or "No metrics generated",
                    # Initialize other core fields as None/empty for consistent structure
                    'antigen': "N/A",
                    'antigen_id': "N/A",
                    'antigen_pdb_chain_id': "N/A",
                    'molecular_weight': None,
                    'molecular_formula': "N/A",
                    'warnings': [],
                    'metrics': {} # Empty metrics dict for failed ones
                }
            all_records.append(record)

        if not all_records:
            logging.warning("No records (successful or failed) to save to CSV.")
            print("Warning: No records to save to CSV.")
            return

        try:
            df = pd.DataFrame(all_records)

            # Flatten the 'metrics' nested dict if it exists and is not empty
            if 'metrics' in df.columns:
                # Normalize only rows where 'metrics' is a dict
                valid_metrics_mask = df['metrics'].apply(lambda x: isinstance(x, dict))
                if valid_metrics_mask.any():
                    metrics_expanded = pd.json_normalize(df.loc[valid_metrics_mask, 'metrics'], sep='_')
                    # df might not have 'metrics' if all failed very early
                    if 'metrics' in df.columns:
                         df = df.drop(columns=['metrics'])
                    df = pd.concat([df, metrics_expanded], axis=1) # Concat might add NaNs for rows without metrics

            # Define core columns that should always be present
            core_columns = ['sequence', 'success', 'error', 'antigen', 'antigen_id', 
                            'antigen_pdb_chain_id', 'molecular_weight', 'molecular_formula', 'warnings']
            
            # Explicitly drop weighted_scores and total_score columns (ranker fields)
            columns_to_drop_ranker = ['weighted_scores', 'total_score']
            for col in columns_to_drop_ranker:
                if col in df.columns:
                    df = df.drop(columns=[col], errors='ignore')

            # Filter metric columns based on self.metrics_to_run
            if self.metrics_to_run is not None:
                allowed_metric_prefixes = [m.value.lower().replace(" ", "_").replace("-", "") for m in self.metrics_to_run]
                
                metric_columns_to_keep = []
                all_current_columns = list(df.columns)

                for col_name in all_current_columns:
                    is_core = col_name in core_columns
                    is_allowed_metric = any(col_name.startswith(prefix + "_") or col_name == prefix for prefix in allowed_metric_prefixes)
                    
                    # Special case: if a category in metrics_to_run has no sub-metrics (e.g., just 'blast' if it's a simple field)
                    # then the prefix itself should be kept if it's a column name.
                    # Example: if 'blast' is in allowed_metric_prefixes and there's a column 'blast'

                    if is_core or is_allowed_metric:
                        metric_columns_to_keep.append(col_name)
                
                # Ensure core columns that might have been dropped if not present in any record are added back if necessary
                # (though DataFrame creation from list of dicts usually handles this by creating NaN columns)
                final_columns_for_df = [col for col in core_columns if col in df.columns] # Start with existing core
                for col in metric_columns_to_keep: # Add allowed metrics
                    if col not in final_columns_for_df and col in df.columns:
                        final_columns_for_df.append(col)
                
                df = df[final_columns_for_df]

            # Convert list-like warnings to string for CSV
            if 'warnings' in df.columns:
                 df['warnings'] = df['warnings'].apply(lambda x: "; ".join(x) if isinstance(x, list) else x)
                 
            # Define the desired column order based on user requirements
            # First part of columns for validation output
            desired_validation_columns = [
                # Basic Identifiers
                "sequence",
                "antigen",
                "antigen_id",
                "antigen_pdb_chain_id",
                
                # Molecular Properties
                "molecular_weight",
                "molecular_formula",

                # General Status / Logs
                "warnings",
                "success",
                "error",
                
                # BLAST
                "blast",
                
                # ProtParam Analysis
                "protparam_molecular_weight",
                "protparam_aromaticity",
                "protparam_instability_index",
                "protparam_isoelectric_point",
                "protparam_gravy",
                "protparam_hydrophobic_amino_acids_percentage",
                "protparam_hydrophilic_amino_acids_percentage",
                "protparam_predicted_solubility",
                "protparam_secondary_structure_fraction",
                
                # Immunogenicity
                'immunogenicity_immunogenic_score',
                'immunogenicity_strong_binding',
                'immunogenicity_moderate_binding',
                'immunogenicity_weak_binding',
                
                # Stability
                "stability_melting_temperature_celsius",
                "stability_normalized_stability_score",
                "stability_details_stabilizing_residues_count",
                "stability_details_destabilizing_residues_count",
                "stability_details_sequence_length",
                "stability_details_unrecognized_residues_count",
                
                
                # Aggregation
                "aggregation_aggregation_propensity",
                "aggregation_aggregation_prone_regions",
                "aggregation_average_aggregation_score",
                
                # Glycosylation
                'glycosylation_n_glycosylation_sites',
                'glycosylation_n_glycosylation_count',
                'glycosylation_o_glycosylation_sites',
                'glycosylation_o_glycosylation_count',
                
                # Structure (Basic)
                "structure_message",
                "structure_project_id",
                "structure_gmqe",
                "structure_pdb_file_path",

                # Structure Model Details
                "structure_model_details_model_id",
                "structure_model_details_status",
                "structure_model_details_qmean",
                "structure_model_details_coordinate_url",
                "structure_model_details_modelcif_url",
                "structure_date_created",
                "structure_project_title",
                "structure_view_url"
                
                
                # Epitope Prediction
                # "epitope_iedb_raw_results_df",
                "epitope_predicted_epitopes",
                "epitope_epitope_sequences_list",
                "epitope_epitope_count",
                "epitope_overall_average_score",
                "epitope_parameters_method",
                "epitope_parameters_sequence_length",
                "epitope_parameters_window_size",

                # Conservancy
                "conservancy_results",
                "conservancy_conservancy_score",
                
                # Developability
                "developability_matched_proteins_df",
                "developability_search_summary_has_active_matches",
                "developability_search_summary_threshold_used",
                "developability_search_summary_search_time_seconds",
                "developability_search_summary_matches_found_count",
                "developability_developability_score",
                "developability_performance_statistics_total_database_entries",
                "developability_performance_statistics_candidates_after_length_filter",
                "developability_performance_statistics_query_sequence_length",
                "developability_performance_statistics_final_matches_found",
                "developability_report_string",

                # Errors
                "immunogenicity_error",
                "binding_affinity_error",
                "structure_error",
                "epitope_error",
                "conservancy_error",
            ]
            
            # Create a list of available columns in the desired order
            ordered_columns = [col for col in desired_validation_columns if col in df.columns]
            
            # Add any remaining columns that might not be in the desired order list
            remaining_columns = [col for col in df.columns if col not in ordered_columns]
            ordered_columns.extend(remaining_columns)
            
            # Reorder the DataFrame columns
            df = df[ordered_columns]

            output_dir_for_csv = os.path.dirname(output_csv_path)
            if output_dir_for_csv: # Ensure the directory exists
                os.makedirs(output_dir_for_csv, exist_ok=True)
            
            df.to_csv(output_csv_path, index=False)
            logging.info(f"Successfully saved validation results for {len(df)} proteins to CSV: {output_csv_path}")
            print(f"INFO: Validation results for {len(df)} proteins saved to {output_csv_path}")
        except Exception as e:
            import traceback
            logging.error(f"Error saving validation results to CSV at {output_csv_path}: {e}\\n{traceback.format_exc()}")
            print(f"ERROR: Could not save validation results to CSV: {e}")

    def save_processing_results(self, results: List[ProcessingResult], output_dir: str):
        """
        Save all results from a processing batch to the output directory.
        This includes individual JSONs for successful metrics, a summary, and a failed details file.
        """
        try:
            os.makedirs(output_dir, exist_ok=True) # Ensure output_dir exists
            
            self.save_batch_summary(results, output_dir) # Save overall summary text file
            
            successful_count = 0
            json_subdir = os.path.join(output_dir, "successful_protein_metrics_json")
            os.makedirs(json_subdir, exist_ok=True)

            for result in results:
                if result.success and result.metrics:
                    # Use a slightly more specific function for saving individual metrics from a batch
                    self.save_individual_metric_result(result.metrics, json_subdir) 
                    successful_count += 1
            
            if successful_count > 0:
                 logging.info(f"💾 Saved {successful_count} individual successful metric JSONs to {json_subdir}")

            failed_details_list = []
            for r in results:
                if not r.success or r.metrics is None:
                    failed_details_list.append({
                        "sequence": r.sequence,
                        "error": r.error or "Metrics not generated"
                    })
            
            if failed_details_list:
                failed_file_path = os.path.join(output_dir, "failed_proteins_batch_details.json")
                with open(failed_file_path, 'w') as f:
                    import json
                    json.dump(failed_details_list, f, indent=2)
                logging.info(f"❌ Saved {len(failed_details_list)} failed protein details to {failed_file_path}")
                
        except Exception as e:
            logging.error(f"❌ Error saving batch processing results: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())

    def _add_warning(self, metrics: Optional[ProteinMetrics], category: str, message: str) -> None:
        """Add a warning to the metrics object if it exists and has a warnings list."""
        if metrics and hasattr(metrics, 'warnings') and isinstance(metrics.warnings, list):
            warning_msg = f"{category}: {message}"
            metrics.warnings.append(warning_msg)
            # Logging the warning here might be too verbose if called many times,
            # consider logging at a higher level summary.
            # logging.warning(f"Protein {metrics.sequence[:10]}... Warning: {warning_msg}") 
        elif metrics: # If metrics object exists but no warnings attribute or not a list
            logging.debug(f"Could not add warning to metrics object for category {category} - 'warnings' attribute missing or not a list.")


if __name__ == "__main__":
    # Example usage for ProteinValidatorV2
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # --- Test Case 1: Basic validation (no antigen-dependent metrics explicitly) ---
    print("--- Test Case 1: Basic Validation (No Antigen Context Set) ---")
    validator_v2_t1 = ProteinValidatorV2(
        pdb_files_path="validation_results_v2/generated_antibody_pdbs_t1",
        metrics_to_run=[MetricCategory.PROTPARAM, MetricCategory.STABILITY] 
    )
    dummy_sequence_1 = "MTQVPSNPPPVVGARHNFSLKECGFKGRYSPTLASARERGYRAVDLLARHGITVSEAFRA"
    print(f"Validating sequence: {dummy_sequence_1[:30]}...")
    results_t1 = validator_v2_t1.validate_protein_list([dummy_sequence_1])
    if results_t1:
        print(f"Metrics for T1: {list(results_t1[0].to_dict()['metrics'].keys())}")
        print(f"Warnings T1: {results_t1[0].warnings}")
        # Binding affinity should be empty or show an error/skip message if it was implicitly run
        print(f"Binding Affinity T1: {results_t1[0].binding_affinity}") 


    # --- Test Case 2: Validation WITH antigen context for Binding Affinity ---
    print("--- Test Case 2: Validation with Antigen Context for Binding Affinity ---")
    validator_v2_t2 = ProteinValidatorV2(
        pdb_files_path="validation_results_v2/generated_antibody_pdbs_t2",
        metrics_to_run=[MetricCategory.PROTPARAM, MetricCategory.BINDING_AFFINITY, MetricCategory.STRUCTURE] # Include STRUCTURE for antibody PDB
    )
    
    # Set antigen context: Using PDB ID and Chain
    print("Setting antigen context for 4R19_A...")
    antigen_set_success = validator_v2_t2.set_antigen_context(
        target_antigen_pdb_chain_id="4R19_A", # e.g. Spike protein part
        antigen_pdb_download_dir="validation_results_v2/downloaded_antigen_pdbs_t2"
    )
    
    if antigen_set_success:
        print(f"Antigen context set. PDB: {validator_v2_t2.target_antigen_pdb_path}, Seq len: {len(validator_v2_t2.target_antigen_sequence or '')}")
        dummy_sequence_2 = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSRYGMDVWGQGTTVTVSS" # Example antibody CDRs
        print(f"Validating antibody sequence: {dummy_sequence_2[:30]}...")
        results_t2 = validator_v2_t2.validate_protein_list([dummy_sequence_2])
        if results_t2:
            print(f"Metrics for T2: {list(results_t2[0].to_dict()['metrics'].keys())}")
            print(f"Binding Affinity T2: {results_t2[0].binding_affinity}")
            print(f"Antibody Structure T2: {results_t2[0].structure.get('pdb_file_path', 'Not generated')}")
            print(f"Warnings T2: {results_t2[0].warnings}")
    else:
        print("Failed to set antigen context for Test Case 2. Binding affinity will not be calculated.")

    # --- Test Case 3: Antigen from sequence only (structure prediction for antigen) ---
    print("--- Test Case 3: Antigen from Sequence (predict structure) + Binding Affinity ---")
    validator_v2_t3 = ProteinValidatorV2(
        pdb_files_path="validation_results_v2/generated_antibody_pdbs_t3",
        metrics_to_run=[MetricCategory.BINDING_AFFINITY, MetricCategory.STRUCTURE]
    )
    antigen_seq_for_pred = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYLGR" # Short example antigen
    print(f"Setting antigen context using sequence for prediction: {antigen_seq_for_pred[:30]}...")
    antigen_set_success_t3 = validator_v2_t3.set_antigen_context(
        target_antigen_sequence=antigen_seq_for_pred,
        antigen_pdb_download_dir="validation_results_v2/predicted_antigen_pdbs_t3"
    )
    if antigen_set_success_t3:
        print(f"Antigen context set. Predicted PDB: {validator_v2_t3.target_antigen_pdb_path}")
        dummy_sequence_3 = "QVQLQESGPGLVKPSQTLSLTCTVSGGSISSYYWSWIRQPPGKGLEWMGYISYSGITTYNPSLKSRVTMLVDTSKNQFSLRLSSVTAADTAVYYCARDRGYSAWFAYWGQGTLVTVSA"
        results_t3 = validator_v2_t3.validate_protein_list([dummy_sequence_3])
        if results_t3:
            print(f"Binding Affinity T3: {results_t3[0].binding_affinity}")
            print(f"Warnings T3: {results_t3[0].warnings}")
    else:
        print("Failed to set antigen context for Test Case 3 (prediction may have failed).")

    # --- Test Case 4: Batch processing from a file ---
    print("--- Test Case 4: Batch Processing from File ---")
    # Create a dummy sequences file for testing
    dummy_sequences_file = "dummy_protein_sequences_v2.txt"
    with open(dummy_sequences_file, "w") as f:
        f.write("MGSTSGSKSSDKGKPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSRYGMDVWGQGTTVTVSS\n")
        f.write("MEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSRYGMDVWGQGTTVTVSSGGGGSGGGGSGGGGSQVQLQESGPGLVKPSQTLSLTCTVSGGSISSYYWSWIRQPPGKGLEWMGYISYSGITTYNPSLKSRVTMLVDTSKNQFSLRLSSVTAADTAVYYCARDRGYSAWFAYWGQGTLVTVSA\n")
        f.write("INVALIDSEQUENCE\n") # Test invalid sequence handling
        f.write(dummy_sequence_1 + "\n") # Reuse a valid sequence

    validator_v2_t4 = ProteinValidatorV2(
        pdb_files_path="validation_results_v2/batch_antibody_pdbs_t4",
        metrics_to_run=None # Run all metrics
    )
    # Set antigen context that will apply to all proteins in the batch for binding affinity
    print("Setting antigen context for batch processing (using 1A2Y_A)...")
    validator_v2_t4.set_antigen_context(
        target_antigen_pdb_chain_id="1A2Y_A", # Example known antigen
        antigen_pdb_download_dir="validation_results_v2/batch_antigen_pdbs_t4"
    )
    
    batch_results_t4 = validator_v2_t4.process_proteins(
        dummy_sequences_file, 
        output_dir="validation_results_v2/batch_run_output_t4"
    )
    print(f"Batch processing T4 complete. {len(batch_results_t4)} results obtained.")
    successful_batch_metrics_t4 = validator_v2_t4.get_successful_metrics(batch_results_t4)
    print(f"{len(successful_batch_metrics_t4)} proteins successfully processed in batch T4.")
    if successful_batch_metrics_t4:
        # Example: Save successful metrics to CSV
        validator_v2_t4.to_csv(batch_results_t4, "validation_results_v2/batch_run_output_t4/successful_metrics.csv")

    # Clean up dummy file
    # os.remove(dummy_sequences_file)
    print(f"\nNote: Dummy file {dummy_sequences_file} was created for testing, please review/delete if necessary.")

    # --- Test Case 5: Using validate_protein_list with file input and CSV output ---
    print("\n--- Test Case 5: validate_protein_list with File Input & CSV Output ---")
    validator_v2_t5 = ProteinValidatorV2(
        pdb_files_path="validation_results_v2/generated_antibody_pdbs_t5",
        metrics_to_run=None # Run all applicable metrics
    )
    # Create a dummy sequences file specifically for this test case
    dummy_file_t5 = "dummy_sequences_for_validate_list_t5.txt"
    with open(dummy_file_t5, "w") as f:
        f.write(dummy_sequence_1 + "\n")
        f.write(dummy_sequence_2 + "\n") # From Test Case 2
        f.write("SHORT\n") # Invalid sequence to test error handling in CSV

    # Set antigen context if binding affinity is to be run seriously
    # For this example, let's assume binding affinity might be run and set a context
    validator_v2_t5.set_antigen_context(
        target_antigen_pdb_chain_id="4R19_A", 
        antigen_pdb_download_dir="validation_results_v2/antigen_pdbs_t5"
    )

    output_csv_t5 = "validation_results_v2/validate_list_output_t5.csv"
    print(f"Calling validate_protein_list with file '{dummy_file_t5}' to output to '{output_csv_t5}'")
    
    successful_metrics_t5 = validator_v2_t5.validate_protein_list(
        input_source=dummy_file_t5,
        output_csv_path=output_csv_t5
    )

    print(f"validate_protein_list (T5) returned {len(successful_metrics_t5)} successful ProteinMetrics objects.")
    if os.path.exists(output_csv_t5):
        print(f"CSV output successfully created at: {output_csv_t5}")
        # You can add a line here to print the contents of the CSV if desired for quick verification
        # import pandas as pd
        # print(pd.read_csv(output_csv_t5))
    else:
        print(f"CSV output was NOT created at: {output_csv_t5}")
    
    # Clean up dummy file for t5
    # os.remove(dummy_file_t5)
    print(f"Note: Dummy file {dummy_file_t5} was created for testing, please review/delete if necessary.")

    pass 