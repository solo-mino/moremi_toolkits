"""
Antibody Ranker V2 (adapted for ProteinValidatorV2)

This module implements a comprehensive scoring system for proteins based on
their physicochemical, structural, and functional properties. It provides normalized 
scoring (0-1) for each property and combines them with weighted scoring for final ranking.
This version is updated to work with ProteinValidatorV2 and respects metrics_to_run.

Key Features:
1. Weighted scoring for different metric categories (respecting skipped metrics).
2. Normalized scoring (0-1) for each property.
3. Comprehensive property evaluation.
4. Detailed scoring breakdown and reporting.
"""

from typing import List, Dict, Optional, Union, Any
import pandas as pd
import numpy as np
from datetime import datetime
import os
from dataclasses import dataclass, field
from .protein_validator_v2 import (
    ProteinMetrics, MetricCategory, ProteinValidatorV2
)
from pathlib import Path
import logging
from tqdm import tqdm
import importlib.util

# Setup basic logging if not already configured by a higher-level script
# This is a library, so ideally logging is configured by the application using it.
# However, for standalone execution or direct use, some basic config is helpful.
if not logging.getLogger().hasHandlers():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Try to import report generators, but don't fail if they're not available
try:
    from .reports.protein_enhanced_report_generator import generate_enhanced_report
    enhanced_report_available = True
except ImportError:
    logging.warning("protein_enhanced_report_generator not found. Enhanced individual PDF/CSV reports may not be generated.")
    enhanced_report_available = False

try:
    from .reports.protein_report_generator import generate_ranking_report as gen_ranking_report
    ranking_report_available = True
except ImportError:
    logging.warning("protein_report_generator (for overall ranking PDF) not found. Overall ranking PDF report may not be generated.")
    ranking_report_available = False

@dataclass
class ScoringConfig:
    """Configuration for protein property scoring calculations"""
    
    category_weights: Dict[MetricCategory, float] = field(default_factory=lambda: {
        MetricCategory.BINDING_AFFINITY: 0.30,
        MetricCategory.STRUCTURE: 0.20,
        MetricCategory.GLYCOSYLATION: 0.10,
        MetricCategory.AGGREGATION: 0.10,
        MetricCategory.PROTPARAM: 0.10,
        MetricCategory.IMMUNOGENICITY: 0.05,
        MetricCategory.CONSERVANCY: 0.05,
        MetricCategory.STABILITY: 0.05,
        MetricCategory.EPITOPE: 0.03,
        MetricCategory.DEVELOPABILITY: 0.02,
        # MetricCategory.BLAST is informational, not typically weighted for ranking.
    })
    
    property_configs: Dict[str, Dict[str, Any]] = field(default_factory=lambda: {
        'dissociation_constant': {
            'optimal_range': (1e-12, 1e-3), 
            'score_func': lambda x: 0.0 if not isinstance(x, (float, int)) or x < 1e-12 or x > 1e-3 else min(1.0, max(0.0, 1 - (np.log10(x) + 3) / 9))
        },
        'gmqe': {
            'optimal_range': (0.6, 1.0),
            'score_func': lambda x: (x - 0.6) / 0.4 if isinstance(x, (float, int)) and 0.6 <= x <= 1.0 else (0.0 if isinstance(x, (float, int)) and x < 0.6 else (1.0 if isinstance(x, (float, int)) and x > 1.0 else 0.0))
        },
        'n_glyc_sites_count': {
            'optimal_range': (1, 4),
            'score_func': lambda x: (x - 1) / 3 if isinstance(x, (float, int)) and 1 <= x <= 4 else (0.0 if isinstance(x, (float, int)) and x < 1 else (1.0 if isinstance(x, (float, int)) and x > 4 else 0.0))
        },
        'aggregation_propensity': {
            'score_func': lambda x: 1.0 if x == 'Low' else (0.5 if x == 'Medium' else 0.0)
        },
        'aggregation_regions': {
            'optimal_range': (0, 8),
            'score_func': lambda x: 1 - (x / 8) if isinstance(x, (float, int)) and 0 <= x <= 8 else (0.0 if isinstance(x, (float, int)) and x > 8 else (1.0 if isinstance(x, (float, int)) and x < 0 else 0.0))
        },
        'gravy': {
            'optimal_range': (-1.5, -0.5),
            'score_func': lambda x: (float(x) + 1.5) / 1.0 if isinstance(x, (float, int)) and -1.5 <= float(x) <= -0.5 else 0.0
        },
        'solubility': {
            'score_func': lambda x: 1.0 if x == 'Soluble' else 0.0
        },
        'immunogenicity_score': {
            'optimal_range': (0, 1),
            'score_func': lambda x: x if isinstance(x, (float, int)) and 0 <= x <= 1 else (0.0 if isinstance(x, (float, int)) and x < 0 else (1.0 if isinstance(x, (float, int)) and x > 1 else 0.0))
        },
        'conservancy_score': {
            'optimal_range': (0, 1),
            'score_func': lambda x: x if isinstance(x, (float, int)) and 0 <= x <= 1 else (0.0 if isinstance(x, (float, int)) and x < 0 else (1.0 if isinstance(x, (float, int)) and x > 1 else 0.0))
        },
        'melting_temperature_celsius': {
            'optimal_range': (65, 90),
            'score_func': lambda x: 0.0 if not isinstance(x, (float, int)) or x < 65 or x > 90 else (x - 65) / 25
        },
        'epitope_score': {
            'optimal_range': (0, 1),
            'score_func': lambda x: x if isinstance(x, (float, int)) and 0 <= x <= 1 else (0.0 if isinstance(x, (float, int)) and x < 0 else (1.0 if isinstance(x, (float, int)) and x > 1 else 0.0))
        },
        'developability_score': {
            'optimal_range': (0, 1),
            'score_func': lambda x: x if isinstance(x, (float, int)) and 0 <= x <= 1 else (0.0 if isinstance(x, (float, int)) and x < 0 else (1.0 if isinstance(x, (float, int)) and x > 1 else 0.0))
        }
    })

class ProteinRanker:
    """Enhanced protein ranker implementing comprehensive scoring system, respecting skipped metrics."""
    
    def __init__(self, config: Optional[ScoringConfig] = None, generate_pdf: bool = False, generate_csv: bool = True):
        self.config = config or ScoringConfig()
        self.df: Optional[pd.DataFrame] = None
        self.reports_dir: Optional[Path] = None
        self.generate_pdf = generate_pdf
        self.generate_csv = generate_csv
        
        try:
            importlib.import_module('.reports.protein_enhanced_report_generator', package='moremi_biokit.proteins')
            importlib.import_module('.reports.protein_report_generator', package='moremi_biokit.proteins')
        except ImportError as e:
            logging.warning(f"Could not import one or more report generators: {e}. Corresponding reports may not be generated.")

    def set_output_directory(self, output_dir: str):
        self.reports_dir = Path(output_dir)
        (self.reports_dir / "protein_reports").mkdir(parents=True, exist_ok=True)
        (self.reports_dir / "rankings").mkdir(parents=True, exist_ok=True)

    def _add_warning_to_metrics(self, metrics: ProteinMetrics, category: str, message: str) -> None:
        if not hasattr(metrics, 'warnings') or metrics.warnings is None:
            metrics.warnings = []
        warning_msg = f"{category}: {message}"
        if warning_msg not in metrics.warnings:
            metrics.warnings.append(warning_msg)
        # Limit logging verbosity here, main warnings will be in ProteinMetrics object

    def _parse_metric_value(self, value: Any, metric_name: str) -> Union[float, str, list, dict]:
        # (Assuming this helper is mostly correct from previous version, minor adjustments for robustness)
        numeric_metrics = [
            'dissociation_constant', 
            'gmqe',
            'n_glyc_sites_count',
            'aggregation_regions',
            'gravy',
            'immunogenicity_score',
            'conservancy_score',
            'melting_temperature_celsius',
            'epitope_score',
            'developability_score'
        ]
        if metric_name in numeric_metrics:
            if isinstance(value, (float, int)):
                return float(value)
            if isinstance(value, str):
                try:
                    return float(value)
                except (ValueError, TypeError):
                    logging.debug(f"Could not convert numeric metric {metric_name} value '{value}' to float, using 0.0")
                    return 0.0
            logging.debug(f"Unexpected type for numeric metric {metric_name} value '{value}', using 0.0")
            return 0.0
        
        if metric_name in ['aggregation_propensity', 'solubility']:
            return str(value) if value is not None else ""
        if isinstance(value, (list, dict)): return value
        if isinstance(value, (float, int)): return float(value)
        if isinstance(value, str):
            try: 
                return float(value)
            except (ValueError, TypeError): 
                return value
            return str(value) if value is not None else ""

    def _normalize_score(self, value: Any, metric_name: str) -> float:
        if isinstance(value, dict):
            if value.get("status") == "Skipped by user configuration.":
                logging.debug(f"Metric {metric_name} was skipped by validator, normalizing to 0.0.")
                return 0.0
            if "error" in value:
                logging.debug(f"Metric {metric_name} has an error from validator: {value.get('error')}, normalizing to 0.0.")
                return 0.0
        
        config_entry = self.config.property_configs.get(metric_name)
        if not config_entry or 'score_func' not in config_entry:
            if isinstance(value, (int, float)): return min(1.0, max(0.0, float(value)))
            # Further handling for unconfigured non-numeric types can be added if necessary
            logging.debug(f"No scoring config for {metric_name} or value type {type(value)} not directly scorable.")
            return 0.0 
            
        parsed_value = self._parse_metric_value(value, metric_name)
        try:
            return round(max(0.0, min(1.0, config_entry['score_func'](parsed_value))), 4)
        except Exception as e:
            logging.warning(f"Error applying score_func for metric {metric_name} with value {parsed_value} (original: {value}): {e}")
            return 0.0

    def calculate_overall_score(self, metrics: ProteinMetrics) -> Dict[str, Any]:
        logging.debug(f"Calculating scores for protein {metrics.sequence[:20]}...")
        
        category_scores_final = {}
        metric_scores_detailed = {}
        
        for cat_enum in MetricCategory:
            cat_name = cat_enum.value
            cat_attr = cat_name.lower().replace(' ', '_').replace('-', '')
            
            metric_data = getattr(metrics, cat_attr, None)
            current_cat_sub_scores = {}
            cat_final_score_value = 0.0

            if metric_data is None or (isinstance(metric_data, dict) and metric_data.get("status") == "Data not found in metrics object"):
                status_msg = "Data not found in metrics object"
                metric_scores_detailed[cat_name] = {"status": status_msg}
                category_scores_final[cat_name] = 0.0
                self._add_warning_to_metrics(metrics, cat_name, status_msg)
                continue
                
            if isinstance(metric_data, dict) and (metric_data.get("status") == "Skipped by user configuration." or "error" in metric_data):
                status_msg = metric_data.get("status", metric_data.get("error", "Skipped/Error by Validator"))
                metric_scores_detailed[cat_name] = {"status": status_msg}
                category_scores_final[cat_name] = 0.0 # Skipped/Errored by validator = 0 score for ranking
                if "error" in status_msg.lower(): 
                    self._add_warning_to_metrics(metrics, cat_name, status_msg)
                    logging.info(f"Category '{cat_name}' for {metrics.sequence[:20]} was '{status_msg}'. Score set to 0.")
                    continue
                    
            if not metric_data and cat_enum != MetricCategory.BLAST: # Empty data for non-BLAST category
                status_msg = "No data available from validator"
                metric_scores_detailed[cat_name] = {"status": status_msg}
                category_scores_final[cat_name] = 0.0
                self._add_warning_to_metrics(metrics, cat_name, status_msg)
                continue

            # --- Score active categories --- 
            try:
                if cat_enum == MetricCategory.PROTPARAM and isinstance(metric_data, dict):
                    current_cat_sub_scores['gravy'] = self._normalize_score(metric_data.get('gravy'), 'gravy')
                    current_cat_sub_scores['solubility'] = self._normalize_score(metric_data.get('predicted_solubility'), 'solubility')
                elif cat_enum == MetricCategory.STABILITY and isinstance(metric_data, dict):
                    current_cat_sub_scores['melting_temperature_celsius'] = self._normalize_score(metric_data.get('melting_temperature_celsius'), 'melting_temperature_celsius')
                elif cat_enum == MetricCategory.AGGREGATION and isinstance(metric_data, dict):
                    current_cat_sub_scores['aggregation_propensity'] = self._normalize_score(metric_data.get('aggregation_propensity'), 'aggregation_propensity')
                    num_regions = len(metric_data.get('aggregation_prone_regions', []))
                    current_cat_sub_scores['aggregation_regions'] = self._normalize_score(num_regions, 'aggregation_regions')
                elif cat_enum == MetricCategory.GLYCOSYLATION and isinstance(metric_data, dict):
                    n_glyc_data = metric_data.get('n_glycosylation')
                    if isinstance(n_glyc_data, dict):
                        n_glyc_count = n_glyc_data.get('count', 0)
                    elif isinstance(n_glyc_data, list):
                        n_glyc_count = len(n_glyc_data)
                    elif isinstance(n_glyc_data, (int, float)):
                        n_glyc_count = n_glyc_data
                    else:
                        n_glyc_count = 0
                        if n_glyc_data is not None:
                             logging.debug(f"Unexpected format for n_glycosylation data: {n_glyc_data}. Using count 0.")
                    current_cat_sub_scores['n_glyc_sites_count'] = self._normalize_score(n_glyc_count, 'n_glyc_sites_count')
                elif cat_enum == MetricCategory.STRUCTURE and isinstance(metric_data, dict):
                    current_cat_sub_scores['gmqe'] = self._normalize_score(metric_data.get('gmqe'), 'gmqe') # ProteinValidatorV2 uses 'gmqe' now
                elif cat_enum == MetricCategory.BINDING_AFFINITY and isinstance(metric_data, dict):
                    current_cat_sub_scores['dissociation_constant'] = self._normalize_score(metric_data.get('dissociation_constant'), 'dissociation_constant')
                elif cat_enum == MetricCategory.IMMUNOGENICITY and isinstance(metric_data, dict):
                    current_cat_sub_scores['score'] = self._normalize_score(metric_data.get('immunogenic_score'), 'immunogenicity_score')
                elif cat_enum == MetricCategory.CONSERVANCY and isinstance(metric_data, dict):
                    current_cat_sub_scores['score'] = self._normalize_score(metric_data.get('conservancy_score'), 'conservancy_score')
                elif cat_enum == MetricCategory.EPITOPE and isinstance(metric_data, dict):
                    current_cat_sub_scores['score'] = self._normalize_score(metric_data.get('overall_average_score'), 'epitope_score')
                elif cat_enum == MetricCategory.DEVELOPABILITY and isinstance(metric_data, dict):
                    current_cat_sub_scores['score'] = self._normalize_score(metric_data.get('developability_score'), 'developability_score')
                elif cat_enum == MetricCategory.BLAST:
                    metric_scores_detailed[cat_name] = {"status": "Collected, not directly scored for ranking"}
                    category_scores_final[cat_name] = 0.0
                    continue # Handled BLAST, skip averaging for it
                else: # Default for other dicts - try to average known sub-metrics or log
                    if isinstance(metric_data, dict):
                        for k, v in metric_data.items():
                            if k in self.config.property_configs: # If sub-metric is known
                                current_cat_sub_scores[k] = self._normalize_score(v, k)
                        if not current_cat_sub_scores:
                             logging.debug(f"Category {cat_name} is a dict but no known sub-metrics scored.")
                             current_cat_sub_scores = {"status": "No scorable sub-metrics"}
                    else: # Non-dict data, non-BLAST, not specifically handled above
                        logging.debug(f"Category {cat_name} data type {type(metric_data)} not specifically handled for scoring.")
                        current_cat_sub_scores={"status": "Unhandled data type"}

                # Calculate category's final score from its sub_scores
                valid_numeric_sub_scores = [s for s in current_cat_sub_scores.values() if isinstance(s, (float, int))]
                if valid_numeric_sub_scores:
                    cat_final_score_value = np.mean(valid_numeric_sub_scores)
                elif "status" not in current_cat_sub_scores:
                     self._add_warning_to_metrics(metrics, cat_name, "No valid numeric sub-metric scores after processing.")
                
                metric_scores_detailed[cat_name] = current_cat_sub_scores
                category_scores_final[cat_name] = round(max(0.0, min(1.0, cat_final_score_value)), 4)

            except Exception as e_score:
                logging.error(f"Error during scoring details of category {cat_name} for {metrics.sequence[:20]}: {e_score}", exc_info=True)
                metric_scores_detailed[cat_name] = {"error_in_ranker_scoring": str(e_score)}
                category_scores_final[cat_name] = 0.0
                self._add_warning_to_metrics(metrics, cat_name, f"Ranker error: {e_score}")

        # Calculate overall weighted score
        overall_score_num = 0.0
        total_weight_den = 0.0
        for cat_enum_member in MetricCategory:
            cat_name_for_weight = cat_enum_member.value
            # Check if this category was SKIPPED by the VALIDATOR
            cat_status_info = metric_scores_detailed.get(cat_name_for_weight, {})
            validator_skipped_or_errored = False
            if isinstance(cat_status_info, dict):
                status_msg = str(cat_status_info.get("status", "")).lower()
                if "skipped by user configuration" in status_msg or "error" in status_msg or "error" in cat_status_info:
                    validator_skipped_or_errored = True
            
            if not validator_skipped_or_errored: # Only consider categories not skipped/errored by validator for weighting
                current_cat_score = category_scores_final.get(cat_name_for_weight, 0.0)
                weight = self.config.category_weights.get(cat_enum_member, 0.0)
                if weight > 0: # And category has a defined weight
                    overall_score_num += current_cat_score * weight
                    total_weight_den += weight
        
        final_overall_score_val = 0.0
        if total_weight_den > 0:
            final_overall_score_val = overall_score_num / total_weight_den
        elif any(s > 0 for s in category_scores_final.values()):
            logging.warning(f"Protein {metrics.sequence[:20]} has positive category scores but total_weight_den is 0. Overall score is 0.")

        final_overall_score_val = round(max(0.0, min(1.0, final_overall_score_val)), 4)

        summary_to_return = {
            'overall_score': final_overall_score_val,
            'category_scores': category_scores_final,
            'metric_scores': metric_scores_detailed,
            'missing_categories': [k for k, v in metric_scores_detailed.items() if isinstance(v, dict) and v.get("status")]
        }
        # self._validate_scores(summary_to_return) # Placeholder if a validation func is needed
        return summary_to_return

    def rank_proteins(self, metrics_list: Union[List[ProteinMetrics], ProteinMetrics]) -> pd.DataFrame:
        if isinstance(metrics_list, ProteinMetrics):
            metrics_list = [metrics_list]
        if not metrics_list: return pd.DataFrame()

        logging.info(f"🎯 Ranking {len(metrics_list)} proteins (Ranker V2 logic)...")
        all_protein_rows = []
        protein_scores_cache = {}

        for i, protein_metric_obj in enumerate(tqdm(metrics_list, desc="📊 Calculating scores (Ranker V2)", unit="protein")):
            try:
                scores = self.calculate_overall_score(protein_metric_obj)
                protein_scores_cache[protein_metric_obj.sequence] = scores
                
                # Store calculated scores back onto the ProteinMetrics object if needed for other reports
                protein_metric_obj.total_score = scores.get('overall_score', 0.0)
                # The 'category_scores' from calculate_overall_score are the UNWEIGHTED category scores
                protein_metric_obj.weighted_scores = {} # This will store the actual weighted scores

                row = {}

                # 1. Basic Information
                row['sequence'] = protein_metric_obj.sequence if protein_metric_obj.sequence else "N/A"
                row['antigen'] = protein_metric_obj.antigen if protein_metric_obj.antigen else "N/A"
                row['antigen_pdb_chain_id'] = protein_metric_obj.antigen_pdb_chain_id if protein_metric_obj.antigen_pdb_chain_id else "N/A"
                row['molecular_formula'] = protein_metric_obj.molecular_formula if protein_metric_obj.molecular_formula else "N/A"
                row['molecular_weight'] = protein_metric_obj.molecular_weight if protein_metric_obj.molecular_weight else "N/A"

                # 2. Scores from `calculate_overall_score` and ProteinMetrics
                row['total_score'] = scores.get('overall_score', 0.0)
                
                # Iterate through MetricCategory to populate scores and raw values
                for cat_enum in MetricCategory:
                    cat_name_for_column = cat_enum.value.lower().replace(" ", "_").replace("-", "") # e.g., "protparam"
                    cat_name_proper = cat_enum.value # e.g., "ProtParam"

                    # Get category data from ProteinMetrics object
                    metric_category_data = getattr(protein_metric_obj, cat_name_for_column, {})
                    if not isinstance(metric_category_data, dict): # Ensure it's a dict for .get()
                        metric_category_data = {}

                    # a. Category Score (unweighted)
                    category_score_value = scores.get('category_scores', {}).get(cat_name_proper, 0.0)
                    row[f'{cat_name_for_column}_score'] = category_score_value

                    # b. Weighted Category Score
                    weight = self.config.category_weights.get(cat_enum, 0.0)
                    weighted_category_score_value = category_score_value * weight
                    row[f'weighted_{cat_name_for_column}_score'] = weighted_category_score_value
                    protein_metric_obj.weighted_scores[cat_name_proper] = weighted_category_score_value

                    # c. Normalized Metric Scores (from scores['metric_scores'])
                    #    And Raw Metric Values (from protein_metric_obj.category_data)
                    normalized_metrics_for_cat = scores.get('metric_scores', {}).get(cat_name_proper, {})
                    if isinstance(normalized_metrics_for_cat, dict) and "status" in normalized_metrics_for_cat:
                        # Category was skipped or had an error, so no individual normalized scores
                        pass
                    elif isinstance(normalized_metrics_for_cat, dict):
                        for norm_metric_key, norm_metric_val in normalized_metrics_for_cat.items():
                            row[f'norm_{norm_metric_key}_score'] = norm_metric_val
                    
                    # Specific Raw Values based on column_naming_guide.md and common practice
                    if cat_enum == MetricCategory.PROTPARAM:
                        row['gravy_value'] = metric_category_data.get('gravy', np.nan)
                        row['solubility'] = metric_category_data.get('predicted_solubility', "Unknown")
                    elif cat_enum == MetricCategory.IMMUNOGENICITY:
                        row['immunogenicity_value'] = metric_category_data.get('immunogenicity_score', np.nan)
                    elif cat_enum == MetricCategory.STABILITY:
                        row['melting_temperature'] = metric_category_data.get('melting_temperature_celsius', np.nan)
                    elif cat_enum == MetricCategory.AGGREGATION:
                        row['aggregation_propensity'] = metric_category_data.get('aggregation_propensity', "Unknown")
                        regions = metric_category_data.get('aggregation_prone_regions', [])
                        row['aggregation_regions_count'] = len(regions) if isinstance(regions, list) else 0
                    elif cat_enum == MetricCategory.GLYCOSYLATION:
                        row['n_glyc_sites_count'] = metric_category_data.get('n_glycosylation', {}).get('count', 0)
                    elif cat_enum == MetricCategory.STRUCTURE:
                        row['gmqe_value'] = metric_category_data.get('gmqe', np.nan)
                    elif cat_enum == MetricCategory.BINDING_AFFINITY:
                        row['dissociation_constant_value'] = metric_category_data.get('dissociation_constant', np.nan)
                    elif cat_enum == MetricCategory.EPITOPE:
                        row['epitope_score_value'] = metric_category_data.get('overall_epitope_score', metric_category_data.get('epitope_score', np.nan))
                    elif cat_enum == MetricCategory.CONSERVANCY:
                        row['conservancy_score_value'] = metric_category_data.get('overall_conservancy_score', np.nan)
                    elif cat_enum == MetricCategory.DEVELOPABILITY:
                        row['developability_value'] = metric_category_data.get('developability_score', np.nan)
                    
                    # For BLAST, the raw data is often a dictionary. We might choose to serialize it or pick key parts.
                    # For now, let's skip detailed BLAST raw values in the main ranking CSV to avoid excessive width,
                    # unless specific raw values are requested by the guide.
                    # The 'blast_score' (category score) will be 0 if not weighted.

                # 3. Warnings and Missing Metrics
                row['warnings_count'] = len(protein_metric_obj.warnings)
                row['warning_details'] = "; ".join(protein_metric_obj.warnings) if protein_metric_obj.warnings else ""
                row['missing_metrics_count'] = len(scores.get('missing_categories', []))
                row['missing_metrics_details'] = "; ".join(scores.get('missing_categories', []))


                all_protein_rows.append(row)
            except Exception as e_rank_loop:
                logging.error(f"Error processing protein {protein_metric_obj.sequence[:20]} for ranking DataFrame: {e_rank_loop}", exc_info=True)
                # Add a basic row with error if this protein fails catastrophically during row construction
                all_protein_rows.append({
                    'sequence': protein_metric_obj.sequence,
                    'total_score': 0.0,
                    'warning_details': f"Ranking Error: {e_rank_loop}"
                })
        
        if not all_protein_rows: return pd.DataFrame()
        self.df = pd.DataFrame(all_protein_rows)
        
        # Drop BLAST related columns as they are not used for ranking
        blast_columns_to_drop = ['blast_score', 'weighted_blast_score']
        # Also, if any raw blast.* columns were inadvertently added, identify and drop them
        # For now, assuming only the score columns are relevant based on current logic.
        # If raw blast data columns (e.g. 'blast_some_value') are present, they would need to be added here.
        self.df = self.df.drop(columns=[col for col in blast_columns_to_drop if col in self.df.columns], errors='ignore')
        
        if 'total_score' in self.df.columns:
            self.df.sort_values('total_score', ascending=False, inplace=True)
            self.df.reset_index(drop=True, inplace=True)
        else:
            logging.error("'total_score' column missing, cannot sort ranks.")
            return self.df # Return unranked
        
        # --- Reporting (condensed, assuming it's similar to previous logic but uses self.df) ---
        if self.reports_dir and self.df is not None and not self.df.empty:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            # Generate individual reports
            # ... (loop through self.df, find matching ProteinMetrics from metrics_list, call generate_protein_report)
            # Save overall ranking CSV
            csv_path = self.save_rankings(self.reports_dir, timestamp)
            if csv_path:
                # Generate overall ranking PDF report
                self.generate_ranking_report(str(csv_path), str(self.reports_dir), timestamp)
        return self.df
    
    # (Keep _safe_serialize, save_rankings, generate_ranking_report, generate_basic_ranking_report, generate_protein_report as before, 
    #  they operate on self.df or passed data and should be largely compatible if column names are maintained)
    def _safe_serialize(self, obj):
        if isinstance(obj, list):
            return [self._safe_serialize(item) for item in obj]
        elif isinstance(obj, dict):
            return {k: self._safe_serialize(v_item) for k, v_item in obj.items()}
        elif isinstance(obj, (np.integer, np.floating, np.bool_)):
            return obj.item() 
        return obj

    def save_rankings(self, output_dir_path: Path, timestamp: str) -> Optional[Path]:
        if self.df is None or self.df.empty: return None
        rankings_subdir = output_dir_path / 'rankings'
        rankings_subdir.mkdir(parents=True, exist_ok=True)
        df_csv = self.df.copy()
        
        # Define the desired column order for ranking output
        desired_ranking_columns = [
            # Basic Identifiers
            "sequence",
            "antigen",
            "antigen_pdb_chain_id",

            # Physical Properties
            "molecular_formula",
            "molecular_weight",

            # Scores (general)
            "total_score",

            # ProtParam-related
            "protparam_score",
            "weighted_protparam_score",
            "norm_gravy_score",
            "gravy_value",
            "norm_solubility_score",
            "solubility",

            # Immunogenicity
            "immunogenicity_score",
            "weighted_immunogenicity_score",
            "norm_score_score",
            "immunogenicity_value",

            # Stability / Melting Temperature
            "stability_score",
            "weighted_stability_score",
            "norm_melting_temperature_celsius_score",
            "melting_temperature",

            # Aggregation
            "aggregation_score",
            "weighted_aggregation_score",
            "norm_aggregation_propensity_score",
            "aggregation_propensity",
            "norm_aggregation_regions_score",
            "aggregation_regions_count",

            # Glycosylation
            "glycosylation_score",
            "weighted_glycosylation_score",
            "norm_n_glyc_sites_count_score",
            "n_glyc_sites_count",

            # Structural Quality
            "structure_score",
            "weighted_structure_score",
            "norm_gmqe_score",
            "gmqe_value",
            
            # Binding Affinity
            "binding_affinity_score",
            "weighted_binding_affinity_score",
            "dissociation_constant_value",
            "norm_dissociation_constant_score",

            # Epitope Score
            "epitope_score",
            "weighted_epitope_score",
            "epitope_score_value",

            # Conservancy
            "conservancy_score",
            "weighted_conservancy_score",
            "conservancy_score_value",
            
            # Developability
            "developability_score",
            "weighted_developability_score",
            "developability_value",

            # Warnings and Missing Data
            "warnings_count",
            "warning_details",
            "missing_metrics_count",
            "missing_metrics_details"
        ]
        
        # Create ordered columns list from those available in the dataframe
        ordered_columns = [col for col in desired_ranking_columns if col in df_csv.columns]
        
        # Add any remaining columns that are not in the desired list
        remaining_columns = [col for col in df_csv.columns if col not in ordered_columns]
        ordered_columns.extend(remaining_columns)
        
        # Reorder dataframe columns
        df_csv = df_csv[ordered_columns]
        
        # Serialize complex objects for CSV output
        for col in df_csv.columns: 
            if df_csv[col].dtype == 'object':
                 df_csv[col] = df_csv[col].apply(lambda x: str(x) if isinstance(x, (dict, list)) else x)
                 
        out_file = rankings_subdir / f'rankings_{timestamp}.csv'
        df_csv.to_csv(out_file, index=False)
        logging.info(f"Rankings saved to {out_file}")
        return out_file

    # generate_protein_report, generate_ranking_report, generate_basic_ranking_report remain largely same
    # but ensure they use self.df and updated data structures correctly
    def generate_protein_report(self, metrics: ProteinMetrics, scores: Dict, output_dir_str: str, rank: int = 0) -> Dict[str, Optional[Path]]:
        # (Code from previous version, ensure it works with current scores dict structure from calculate_overall_score)
        # This function needs to correctly access scores['category_scores'] and scores['overall_score']
        # and use MetricCategory enum keys for self.config.category_weights if accessing directly.
        # For protein_data["category_scores"]["weighted"], ensure keys are strings if report generator expects that.
        # ... (ensure this method is robust)
        # Simplified for brevity, assume previous implementation is adapted for new scores structure
        if not enhanced_report_available: return {"pdf": None, "csv": None}
        protein_reports_dir = Path(output_dir_str) / "protein_reports"
        protein_reports_dir.mkdir(parents=True, exist_ok=True)
        # protein_data prep needs careful check with new `scores` structure
        protein_data = { "sequence": metrics.sequence, "rank": rank, "total_score": scores.get('overall_score',0.0), 
                           "category_scores_raw": scores.get('category_scores', {}), "warnings": metrics.warnings }
        # ... (populate other fields for report from metrics and scores['metric_scores']) ... 
        try:
            return generate_enhanced_report(protein_data, str(protein_reports_dir), self.generate_pdf, self.generate_csv)
        except Exception as e_report: 
            logging.error(f"Enhanced report gen failed for {metrics.sequence[:10]}: {e_report}")
            return {"pdf": None, "csv": None}

    def generate_ranking_report(self, csv_file_str: str, output_dir_str: str, timestamp: str) -> Optional[Path]:
        # PDF generation is commented out as per user request.
        # if not self.generate_pdf:
        #     logging.info("PDF generation for overall ranking report is disabled.")
        #     return None

        # if not ranking_report_available and not hasattr(self, 'generate_basic_ranking_report'): return None
        # rankings_dir_path = Path(output_dir_str) / "rankings"
        # rankings_dir_path.mkdir(parents=True, exist_ok=True)
        # pdf_out_path = rankings_dir_path / f"ranking_report_{timestamp}.pdf"
        # try:
        #     if ranking_report_available:
        #         gen_ranking_report(csv_file_str, str(pdf_out_path), top_n=None)
        #     else:
        #         self.generate_basic_ranking_report(csv_file_str, str(pdf_out_path))
        #     return pdf_out_path
        # except Exception as e_rank_report:
        #     logging.error(f"Ranking PDF report failed: {e_rank_report}")
        logging.info("Overall ranking PDF report generation is currently disabled.")
        return None # Ensure it returns None as no PDF path is generated

    def generate_basic_ranking_report(self, csv_file: str, output_path_str: str):
        # (Implementation from previous version, assumed mostly correct)
        try:
            from fpdf import FPDF # Ensure fpdf is a dependency
            df_rep = pd.read_csv(csv_file)
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font('Arial', 'B', 12)
            pdf.cell(0,10, 'Basic Protein Ranking Report',0,1,'C')
            # ... (Add more content to basic report from df_rep)
            pdf.output(output_path_str)
            return True
        except Exception as e_basic_rep:
            logging.error(f"Basic PDF report failed: {e_basic_rep}")
            return False

def rank_proteins_from_metrics(
    protein_sequences_input: Union[str, List[str], Path],
    output_dir: str,
    config: Optional[ScoringConfig] = None,
    generate_pdf: bool = False,
    generate_csv: bool = True,
    metrics_to_run: Optional[List[MetricCategory]] = None,
    target_antigen_sequence: Optional[str] = None,
    target_antigen_pdb_file_path: Optional[str] = None,
    target_antigen_pdb_chain_id: Optional[str] = None,
    antigen_pdb_download_path: Optional[str] = None 
) -> pd.DataFrame:
    """
    Main function to rank proteins using ProteinValidatorV2 and ProteinRanker.
    """
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    main_output_dir = Path(output_dir) / f"ranking_run_{run_timestamp}"
    main_output_dir.mkdir(parents=True, exist_ok=True)

    log_file = main_output_dir / "ranking_process.log"
    # Reconfigure logging for this run to go to the run-specific directory
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]: root_logger.removeHandler(handler)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file)]
    )

    logging.info(f"Ranking process started. Output Dir: {main_output_dir}")
    logging.info(f"Enhanced report generator available: {enhanced_report_available}")
    logging.info(f"Ranking report generator available: {ranking_report_available}")
    
    # Prepare config for dynamic weight assignment if metrics_to_run is specified
    processed_config = config
    if metrics_to_run:
        processed_config = config or ScoringConfig() # Ensure we have a config object
        num_metrics = len(metrics_to_run)
        if num_metrics > 0:
            equal_weight = 1.0 / num_metrics
            new_category_weights = {category_enum: 0.0 for category_enum in MetricCategory} # Initialize all to 0
            for metric_category_enum in metrics_to_run:
                if isinstance(metric_category_enum, MetricCategory):
                    new_category_weights[metric_category_enum] = equal_weight
                else:
                    # Attempt to convert string to MetricCategory enum if needed (robustness)
                    try:
                        # This assumes MetricCategory values are simple strings like "Structure", "Binding Affinity"
                        # This might need adjustment based on how MetricCategory enum is defined and how metrics_to_run strings are formatted
                        found_enum = False
                        for enum_member in MetricCategory:
                            if enum_member.value.lower() == str(metric_category_enum).lower():
                                new_category_weights[enum_member] = equal_weight
                                found_enum = True
                                break
                        if not found_enum:
                            logging.warning(f"Could not map '{metric_category_enum}' to a MetricCategory enum for dynamic weighting. It will be ignored.")
                    except Exception as e:
                        logging.warning(f"Error mapping '{metric_category_enum}' to MetricCategory: {e}. It will be ignored.")
            
            # Log the dynamically assigned weights
            logging.info(f"Dynamically assigning equal weights to {num_metrics} specified metrics:")
            for cat, wt in new_category_weights.items():
                if wt > 0:
                    logging.info(f"  - {cat.value}: {wt:.4f}")

            processed_config.category_weights = new_category_weights
        else: # metrics_to_run is an empty list
            logging.warning("'metrics_to_run' was provided as an empty list. Defaulting to zero weights for all categories unless overridden in a custom ScoringConfig.")
            # Ensure category_weights is at least an empty dict or all zeros if we want to enforce it
            if processed_config: # Should exist due to line above
                 processed_config.category_weights = {category_enum: 0.0 for category_enum in MetricCategory}

    ranker_instance = ProteinRanker(processed_config, generate_pdf=generate_pdf, generate_csv=generate_csv)
    ranker_instance.set_output_directory(str(main_output_dir)) # Ranker outputs to subfolders here

    validator_antibody_pdb_output_path = main_output_dir / "validator_antibody_pdbs"
    validator_antibody_pdb_output_path.mkdir(parents=True, exist_ok=True)

    effective_antigen_pdb_path = Path(antigen_pdb_download_path) if antigen_pdb_download_path else main_output_dir / "validator_antigen_pdbs"
    effective_antigen_pdb_path.mkdir(parents=True, exist_ok=True)

    validator = ProteinValidatorV2(
        pdb_files_path=str(validator_antibody_pdb_output_path),
        metrics_to_run=metrics_to_run
    )

    if target_antigen_sequence or target_antigen_pdb_file_path or target_antigen_pdb_chain_id:
        logging.info("Attempting to set antigen context for validator...")
        antigen_set = validator.set_antigen_context(
        target_antigen_sequence=target_antigen_sequence,
        target_antigen_pdb_file_path=target_antigen_pdb_file_path,
        target_antigen_pdb_chain_id=target_antigen_pdb_chain_id,
            antigen_pdb_download_dir=str(effective_antigen_pdb_path)
        )
        logging.info(f"Antigen context set: {antigen_set}")
    else:
        logging.info("No antigen parameters provided for validator.")

    validator_overall_csv_path = main_output_dir / "validation_attempts_summary.csv"
    # ProteinValidatorV2.validate_protein_list returns List[ProteinMetrics] of successful ones.
    # It also internally calls its own `to_csv` with all ProcessingResult if output_csv_path is given.
    metrics_for_ranking = validator.validate_protein_list(
        input_source=protein_sequences_input,
        output_csv_path=str(validator_overall_csv_path) 
    )

    # Reporting on validation phase (summary based on validator's output CSV)
    total_attempted_validation = 0
    failures_in_validation = 0
    if validator_overall_csv_path.exists():
        try:
            val_df = pd.read_csv(validator_overall_csv_path)
            total_attempted_validation = len(val_df)
            failures_in_validation = len(val_df[val_df['success'] == False]) if 'success' in val_df.columns else total_attempted_validation - len(metrics_for_ranking)
        except Exception as e_val_csv:
            logging.warning(f"Could not parse validation summary CSV {validator_overall_csv_path}: {e_val_csv}")
            # Estimate from input if possible
            if isinstance(protein_sequences_input, list): total_attempted_validation = len(protein_sequences_input)
            # else: difficult to estimate from file path here without re-reading
            failures_in_validation = total_attempted_validation - len(metrics_for_ranking) if total_attempted_validation >= len(metrics_for_ranking) else 0

    logging.info(f"Validation Phase: Total Attempted: {total_attempted_validation}, Successful for Ranking: {len(metrics_for_ranking)}, Failed/Skipped Validation: {failures_in_validation}")

    if not metrics_for_ranking:
        logging.warning("No proteins successfully validated. Ranking cannot proceed.")
        return pd.DataFrame()

    final_ranked_df = ranker_instance.rank_proteins(metrics_for_ranking)

    return final_ranked_df

if __name__ == "__main__":
    # Create dummy file for testing
    dummy_protein_file = "dummy_proteins_for_ranker_main.txt"
    with open(dummy_protein_file, "w") as f:
        f.write(">seq1\nMTQVPSNPPPVVGARHNFSLKECGFKGRYSPTLASARERGYRAVDLLARHGITVSEAFRA\n")
        f.write(">seq2_antibody_example\nEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSRYGMDVWGQGTTVTVSS\n")
        f.write("SHORT\n") # Invalid seq
    
    test_output_dir = "./ranking_tool_test_outputs"

    # Test Case 1: Run with a subset of metrics, and antigen context
    logging.info("--- Running Ranker Test Case 1: Specific Metrics + Antigen ---")
    rank_df_1 = rank_proteins_from_metrics(
        protein_sequences_input=dummy_protein_file,
        output_dir=os.path.join(test_output_dir, "test_run_specific_metrics"),
        metrics_to_run=[
            MetricCategory.PROTPARAM,
            MetricCategory.STABILITY,
            MetricCategory.BINDING_AFFINITY, # This will require antigen context
            MetricCategory.STRUCTURE # For antibody structure for binding affinity
        ],
        target_antigen_pdb_chain_id="1A2Y_A", # Example, ensure this PDB exists or can be fetched
        antigen_pdb_download_path=os.path.join(test_output_dir, "test_antigens_cache"),
        generate_pdf=False, generate_csv=True
    )
    print("--- Test Case 1 Results ---")
    if not rank_df_1.empty:
        print(f"Ranked {len(rank_df_1)} proteins.")
        print(rank_df_1[['sequence', 'total_score']].head())
    else:
        print("Test Case 1: No proteins were ranked (check logs).")

    # Test Case 2: Run with all metrics (default), no antigen context (binding affinity should be skipped/low score)
    logging.info("--- Running Ranker Test Case 2: All Metrics, No Antigen ---")
    rank_df_2 = rank_proteins_from_metrics(
        protein_sequences_input=dummy_protein_file,
        output_dir=os.path.join(test_output_dir, "test_run_all_metrics_no_antigen"),
        metrics_to_run=None, # All metrics
        # No antigen context provided
    )
    print("--- Test Case 2 Results ---")
    if not rank_df_2.empty:
        print(f"Ranked {len(rank_df_2)} proteins.")
        print(rank_df_2[['sequence', 'total_score']].head())
    else:
        print("Test Case 2: No proteins were ranked (check logs).")

    # os.remove(dummy_protein_file) # Clean up
    print(f"Test dummy file {dummy_protein_file} can be removed.")
