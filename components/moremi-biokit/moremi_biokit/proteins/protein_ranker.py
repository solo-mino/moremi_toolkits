"""
Antibody Ranker V1

This module implements a comprehensive scoring system for proteins based on
their physicochemical, structural, and functional properties. It provides normalized 
scoring (0-1) for each property and combines them with weighted scoring for final ranking.

Key Features:
1. Weighted scoring for different metric categories
2. Normalized scoring (0-1) for each property
3. No binary pass/fail criteria
4. Comprehensive property evaluation
5. Detailed scoring breakdown
"""

from typing import List, Dict, Optional, Union, Any
import pandas as pd
import numpy as np
from datetime import datetime
import os
from dataclasses import dataclass
from .protein_validator import (
    ProteinMetrics, MetricCategory, MetricRanges, ProteinValidator
)
from pathlib import Path
import logging
from tqdm import tqdm
import importlib.util
import sys

# Try to import report generators, but don't fail if they're not available
try:
    from .reports.protein_enhanced_report_generator import generate_enhanced_report
    enhanced_report_available = True
except ImportError:
    logging.warning("protein_enhanced_report_generator not found")
    enhanced_report_available = False

try:
    from .reports.protein_report_generator import generate_ranking_report as gen_ranking_report
    ranking_report_available = True
except ImportError:
    logging.warning("protein_report_generator not found")
    ranking_report_available = False

@dataclass
class ScoringConfig:
    """Configuration for protein property scoring calculations"""
    
    # Weights for different categories based on importance
    category_weights = {
        MetricCategory.BINDING_AFFINITY: 0.30,  # Most important
        MetricCategory.STRUCTURE: 0.20,         # Second most important
        MetricCategory.GLYCOSYLATION: 0.10,     # Medium importance
        MetricCategory.AGGREGATION: 0.10,       # Medium importance
        MetricCategory.PROTPARAM: 0.10,         # Medium importance
        MetricCategory.IMMUNOGENICITY: 0.05,    # Lower importance
        MetricCategory.CONSERVANCY: 0.05,       # Lower importance
        MetricCategory.STABILITY: 0.05,         # Lower importance
        MetricCategory.EPITOPE: 0.03,           # Lower importance
        MetricCategory.DEVELOPABILITY: 0.02     # Lowest importance
    }
    
    # Property scoring functions and ranges
    property_configs = {
        # Binding Affinity (0.30)
        # 'binding_affinity': {
        #     'optimal_range': (1e-12, 1e-3),  # e-12 to e-3, closer to e-12 is better
        #     'score_func': lambda x: 0.0 if eval(x) < 1e-12 or eval(x) > 1e-3 else 1 - (np.log10(eval(x)) + 3) / 9  # Normalize log(Kd)
        # },
        'dissociation_constant': {
            'optimal_range': (1e-12, 1e-3),  # e-12 to e-3, closer to e-12 is better
            'score_func': lambda x: 0.0 if eval(x) < 1e-12 or eval(x) > 1e-3 else 1 - (np.log10(eval(x)) + 3) / 9  # Normalize log(Kd)
        },
        
        # Structure (0.20)
        'gmqe_score': {
            'optimal_range': (0.6, 1.0),  # 0.6 to 1.0, closer to 1 is better
            'score_func': lambda x: (x - 0.6) / 0.4 if 0.6 <= x <= 1.0 else 0.0 if x < 0.6 else 1.0  # Normalize to 0-1
        },
        
        # Glycosylation (0.10)
        'n_glyc_sites_count': {
            'optimal_range': (1, 4),  # 1 to 4, closer to 4 is better
            'score_func': lambda x: (x - 1) / 3 if 1 <= x <= 4 else 0.0 if x < 1 else 1.0  # Normalize to 0-1
        },
        
        # Aggregation (0.10)
        'aggregation_propensity': {
            # Low, Medium, High - Low is better
            'score_func': lambda x: 1.0 if x == 'Low' else (0.5 if x == 'Medium' else 0.0)
        },
        'aggregation_regions': {
            'optimal_range': (0, 8),  # 0 to 8, closer to 0 is better
            'score_func': lambda x: 1 - (x / 8) if 0 <= x <= 8 else 0.0 if x > 8 else 1.0  # Normalize to 0-1
        },
        
        # ProtParam (0.10)
        'gravy': {
            'optimal_range': (-1.5, -0.5),  # -1.5 to -0.5, closer to -0.5 is better
            # 'score_func': lambda x: (float(x) + 1.5) / 1.0 if -1.5 <= float(x) <= -0.5 else 0.0 if float(x) < -1.5 else 1.0  # Normalize to 0-1
            'score_func': lambda x: (float(x) + 1.5) / 1.0 if -1.5 <= float(x) <= -0.5 else 0.0  # Normalize to 0-1
        },
        'solubility': {
            # Soluble/Not Soluble - Soluble is better
            'score_func': lambda x: 1.0 if x == 'Soluble' else 0.0
        },
        
        # Immunogenicity (0.05)
        'immunogenicity_score': {
            'optimal_range': (0, 1),  # 0 to 1, closer to 1 is better
            'score_func': lambda x: x if 0 <= x <= 1 else 0.0 if x < 0 else 1.0  # Higher is better
        },
        
        # Conservancy (0.05)
        'conservancy_score': {
            'optimal_range': (0, 1),  # 0 to 1, closer to 1 is better
            'score_func': lambda x: x if 0 <= x <= 1 else 0.0 if x < 0 else 1.0  # Higher is better
        },
        
        # Stability (0.05)
        'melting_temperature': {
            'optimal_range': (65, 90),  # 65°C to 90°C, closer to 90°C is better
            'score_func': lambda x: 0.0 if x < 65 or x > 90 else (x - 65) / 25  # Normalize to 0-1
        },
        
        # Epitope (0.03)
        'epitope_score': {
            'optimal_range': (0, 1),  # 0 to 1, closer to 1 is better
            'score_func': lambda x: x if 0 <= x <= 1 else 0.0 if x < 0 else 1.0  # Higher is better
        },
        
        # Developability (0.02)
        'developability_score': {
            'optimal_range': (0, 1),  # 0 to 1, closer to 1 is better
            'score_func': lambda x: x if 0 <= x <= 1 else 0.0 if x < 0 else 1.0  # Higher is better
        }
    }

class ProteinRanker:
    """Enhanced protein ranker implementing comprehensive scoring system"""
    
    def __init__(self, config: Optional[ScoringConfig] = None, generate_pdf: bool = True, generate_csv: bool = True):
        """Initialize ranker with configuration and report generation flags.

        Args:
            config (Optional[ScoringConfig]): Scoring configuration.
            generate_pdf (bool): Whether to generate PDF reports for individual protein.
            generate_csv (bool): Whether to generate CSV reports for individual protein.
        """
        self.config = config or ScoringConfig()
        self.ranges = MetricRanges()
        self.df = None
        self.reports_dir = None
        self.generate_pdf = generate_pdf
        self.generate_csv = generate_csv
        
        # Verify that report generators can be imported
        try:
            importlib.import_module('.reports.protein_enhanced_report_generator', package='moremi_biokit.proteins')
            importlib.import_module('.reports.protein_report_generator', package='moremi_biokit.proteins')
            logging.debug("Report generator modules successfully imported")
        except ImportError as e:
            logging.warning(f"Could not import report generators: {str(e)}")
            logging.warning("PDF report generation may not work correctly")

    def set_output_directory(self, output_dir: str):
        """Set the output directory for reports"""
        self.reports_dir = Path(output_dir)
        # Create directories for individual protein reports and rankings
        (self.reports_dir / "protein_reports").mkdir(parents=True, exist_ok=True)
        (self.reports_dir / "rankings").mkdir(parents=True, exist_ok=True)

    def generate_protein_report(self, metrics: ProteinMetrics, scores: Dict, output_dir: str, rank: int = 0) -> Dict[str, Optional[Path]]:
        """Generate enhanced report for a single protein based on instance flags.

        Returns:
            Dict[str, Optional[Path]]: Dictionary of generated report paths ('pdf', 'csv').
        """
        try:
            if not metrics or not metrics.sequence:
                raise ValueError("Invalid metrics object or missing sequence")

            # Create the main output directories
            protein_reports_dir = os.path.join(output_dir, "protein_reports")
            os.makedirs(protein_reports_dir, exist_ok=True)
            
            # Check if enhanced report generator is available
            if not enhanced_report_available:
                logging.warning("Enhanced report generator not available. Skipping individual report generation.")
                return {"pdf": None, "csv": None}

            # Prepare data for report - safely serialize all complex data types
            protein_data = {
                "sequence": metrics.sequence,
                "antigen": metrics.antigen,
                "molecular_formula": metrics.molecular_formula,
                "molecular_weight": metrics.molecular_weight,
                "rank": rank,
                "blast": self._safe_serialize(metrics.blast) if hasattr(metrics, 'blast') else [],
                "protparam": self._safe_serialize(metrics.protparam) if hasattr(metrics, 'protparam') else {},
                "immunogenicity": self._safe_serialize(metrics.immunogenicity) if hasattr(metrics, 'immunogenicity') else {},
                "stability": self._safe_serialize(metrics.stability) if hasattr(metrics, 'stability') else {},
                "aggregation": self._safe_serialize(metrics.aggregation) if hasattr(metrics, 'aggregation') else {},
                "glycosylation": self._safe_serialize(metrics.glycosylation) if hasattr(metrics, 'glycosylation') else {},
                "binding_affinity": self._safe_serialize(metrics.binding_affinity) if hasattr(metrics, 'binding_affinity') else {},
                "structure": self._safe_serialize(metrics.structure) if hasattr(metrics, 'structure') else {},
                "epitope": self._safe_serialize(metrics.epitope) if hasattr(metrics, 'epitope') else {},
                "developability": self._safe_serialize(metrics.developability) if hasattr(metrics, 'developability') else {},
                "conservancy": self._safe_serialize(metrics.conservancy) if hasattr(metrics, 'conservancy') else {},
                "warnings": []  # Add warnings list for any issues
            }
            
            # Add raw and weighted scores
            protein_data["category_scores"] = {
                "raw": {
                    key.replace(' ', '_').lower(): value 
                    for key, value in scores.get('category_scores', {}).items()
                },
                "weighted": {
                    key.replace(' ', '_').lower(): value * self.config.category_weights.get(MetricCategory(key), 0.0)
                    for key, value in scores.get('category_scores', {}).items()
                }
            }
            protein_data["total_score"] = scores.get('overall_score', 0.0)

            # Generate enhanced report for individual protein
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            try:
                # First attempt with full data
                report_paths = generate_enhanced_report(
                    protein_data, 
                    protein_reports_dir, 
                    generate_pdf=self.generate_pdf, 
                    generate_csv=self.generate_csv
                )
                logging.info(f"Generated enhanced report for protein {metrics.sequence[:20]}... in {protein_reports_dir}")
                return report_paths
            except Exception as e:
                logging.error(f"Error in enhanced report generation: {str(e)}")
                
                # Second attempt with simplified data
                try:
                    # Create a simplified version with only essential data
                    simplified_data = {
                        "sequence": metrics.sequence,
                        "antigen": metrics.antigen,
                        "molecular_formula": metrics.molecular_formula,
                        "molecular_weight": metrics.molecular_weight,
                        "rank": rank,
                        "category_scores": protein_data["category_scores"],
                        "total_score": protein_data["total_score"],
                        "warnings": [f"Report simplified due to error: {str(e)}"]
                    }
                    
                    report_paths = generate_enhanced_report(
                        simplified_data, 
                        protein_reports_dir,
                        generate_pdf=self.generate_pdf, 
                        generate_csv=self.generate_csv
                    )
                    logging.warning(f"Generated simplified report for protein {metrics.sequence[:20]}...")
                    return report_paths
                except Exception as e2:
                    logging.error(f"Failed to generate simplified report: {str(e2)}")
                    return {"pdf": None, "csv": None}
            
        except Exception as e:
            logging.error(f"Failed to generate report for protein {metrics.sequence[:20]}...")
            logging.error(f"Error details: {str(e)}")
            return {"pdf": None, "csv": None}

    def _normalize_score(self, value: Union[float, int, str, dict, list], metric_name: str) -> float:
        """Normalize a metric value to a score between 0 and 1"""
        # Handle the case when the metric is not configured for scoring
        if metric_name not in self.config.property_configs:
            # Handle different value types
            if isinstance(value, (int, float)):
                return min(1.0, max(0.0, float(value)))
            elif isinstance(value, str):
                try:
                    return min(1.0, max(0.0, float(value)))
                except ValueError:
                    return 0.0  # Default score for strings that can't be converted
            elif isinstance(value, dict):
                # For dictionaries, use the average of values if they're numeric
                numeric_values = []
                for k, v in value.items():
                    try:
                        if isinstance(v, (int, float)):
                            numeric_values.append(float(v))
                        elif isinstance(v, str):
                            numeric_values.append(float(v))
                    except (ValueError, TypeError):
                        pass
                if numeric_values:
                    return min(1.0, max(0.0, sum(numeric_values) / len(numeric_values)))
                return 0.0  # Default score for dictionaries with no numeric values
            elif isinstance(value, list):
                # For lists, give a score based on presence
                return 1.0 if value else 0.0
            else:
                return 0.0  # Default score for other types
            
        # Use the configured scoring function if available
        config = self.config.property_configs[metric_name]
        
        try:
            # Handle string values (like 'Soluble', 'Low', etc.)
            if isinstance(value, str):
                # For string values representing numeric metrics, try to convert
                if metric_name in ['gravy', 'gmqe_score', 'melting_temperature', 'n_glyc_sites_count', 'immunogenicity_score', 
                                   'conservancy_score', 'epitope_score', 'developability_score', 'aggregation_regions']:
                    try:
                        value = float(value)
                    except (ValueError, TypeError):
                        # If conversion fails, it's probably a non-numeric string like 'Soluble'
                        pass
                return config['score_func'](value)
            # Handle numeric values and ensure they're within range
            elif isinstance(value, (int, float)):
                # Check if the value is numeric and apply scoring function
                return config['score_func'](float(value))
            else:
                # Default fallback
                return 0.0  # Neutral score
        except (ValueError, TypeError, AttributeError) as e:
            # If the scoring function fails, return a default score
            logging.warning(f"Error scoring metric {metric_name} with value {value}: {str(e)}")
            return 0.0  # Neutral score

    def _calculate_category_score(self, scores: Dict[str, float]) -> float:
        """Calculate the average score for a category"""
        if not scores:
            return 0.0
        return np.mean(list(scores.values()))

    def _validate_scores(self, scores: Dict[str, Any]) -> None:
        """Validate that all scores are within the valid range [0,1]"""
        # Check overall score
        if not 0 <= scores['overall_score'] <= 1:
            logging.warning(f"Overall score {scores['overall_score']} outside valid range [0,1]")
            scores['overall_score'] = max(0.0, min(1.0, scores['overall_score']))
        
        # Check category scores
        for category, score in scores['category_scores'].items():
            if not 0 <= score <= 1:
                logging.warning(f"Category {category} score {score} outside valid range [0,1]")
                scores['category_scores'][category] = max(0.0, min(1.0, score))
        
        # Check metric scores
        for category, metrics in scores['metric_scores'].items():
            for metric, score in metrics.items():
                if not 0 <= score <= 1:
                    logging.warning(f"Metric {metric} score {score} outside valid range [0,1]")
                    scores['metric_scores'][category][metric] = max(0.0, min(1.0, score))

    def calculate_overall_score(self, metrics: ProteinMetrics) -> Dict[str, Any]:
        """Calculate overall score and category scores for an protein"""
        logging.debug(f"Calculating scores for protein {metrics.sequence[:20]}...")
        
        category_scores = {}
        metric_scores = {}
        
        # Special handling for categories with multiple metrics
        # Aggregation - has 'aggregation_propensity' and 'aggregation_regions'
        aggregation_metrics = {}
        if hasattr(metrics, 'aggregation') and metrics.aggregation:
            agg_propensity = metrics.aggregation.get('aggregation_propensity', 'High')  # Default to worst
            agg_regions = len(metrics.aggregation.get('aggregation_prone_regions', []))
            
            # Score for aggregation propensity (Low/Medium/High)
            aggregation_metrics['aggregation_propensity'] = self._normalize_score(agg_propensity, 'aggregation_propensity')
            
            # Score for aggregation regions count (0-8, fewer is better)
            aggregation_metrics['aggregation_regions'] = self._normalize_score(agg_regions, 'aggregation_regions')
            
            # Calculate average for Aggregation category
            if aggregation_metrics:
                category_scores['Aggregation'] = sum(aggregation_metrics.values()) / len(aggregation_metrics)
                metric_scores['Aggregation'] = aggregation_metrics

        # ProtParam - has 'gravy' and 'solubility'
        protparam_metrics = {}
        if hasattr(metrics, 'protparam') and metrics.protparam:
            # Convert gravy to float to avoid comparison issues
            gravy_str = metrics.protparam.get('gravy', '-1.0')  # Default to middle
            try:
                gravy = float(gravy_str)
            except (ValueError, TypeError):
                gravy = 0.0  # Default if conversion fails
                
            solubility = metrics.protparam.get('predicted_solubility', 'Not Soluble')  # Default to worst
            
            # Score for GRAVY (-1.5 to -0.5, closer to -0.5 is better)
            protparam_metrics['gravy'] = self._normalize_score(gravy, 'gravy')
            
            # Score for solubility (Soluble/Not Soluble)
            protparam_metrics['solubility'] = self._normalize_score(solubility, 'solubility')
            
            # Calculate average for ProtParam category
            if protparam_metrics:
                category_scores['ProtParam'] = sum(protparam_metrics.values()) / len(protparam_metrics)
                metric_scores['ProtParam'] = protparam_metrics
        
        # Calculate scores for each remaining category
        for category in MetricCategory:
            # Skip already processed categories
            if category.value in category_scores:
                continue
                
            try:
                category_name = category.value
                category_attr = category_name.lower().replace(' ', '_').replace('-', '')
                
                # Check if the category attribute exists
                if not hasattr(metrics, category_attr):
                    logging.warning(f"Category {category_name} not found in metrics")
                    continue
                    
                category_metrics = getattr(metrics, category_attr)
                
                # Skip if no metrics for this category
                if not category_metrics:
                    continue
                    
                metric_scores[category_name] = {}
                
                # Handle different types of metrics storage
                if isinstance(category_metrics, dict):
                    # Process dictionary metrics
                    for metric_name, value in category_metrics.items():
                        try:
                            # Try to use the primary metric directly if possible
                            if category_name == 'Binding Affinity' and 'dissociation_constant' in category_metrics:
                                kd_value = category_metrics['dissociation_constant']
                                # For Binding Affinity, we need to ensure we handle very small values (like 4.5e-53) correctly
                                if not isinstance(kd_value, (int, float)):
                                    try:
                                        kd_value = float(kd_value)
                                    except (ValueError, TypeError):
                                        kd_value = 0.0  # Default to 0 if invalid value
                                
                                # Apply the scoring function following the criteria in val_rank.md
                                if kd_value < 1e-12 or kd_value > 1e-3:
                                    score = 0.0  # Outside the valid range
                                else:
                                    # For values in range, normalize: closer to 1e-12 is better
                                    score = 1 - (np.log10(kd_value) + 3) / 9
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                
                                metric_scores[category_name]['dissociation_constant'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Structure' and 'gmqe_score' in category_metrics:
                                gmqe = category_metrics['gmqe_score']
                                score = self._normalize_score(gmqe, 'gmqe_score')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['gmqe_score'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Glycosylation' and 'n_glyc_sites_count' in category_metrics:
                                nglyc = category_metrics['n_glyc_sites_count']
                                score = self._normalize_score(nglyc, 'n_glyc_sites_count')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['n_glyc_sites_count'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Stability' and 'melting_temperature' in category_metrics:
                                tm = category_metrics['melting_temperature']
                                score = self._normalize_score(tm, 'melting_temperature')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['melting_temperature'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Immunogenicity' and 'score' in category_metrics:
                                imm_score = category_metrics['score']
                                score = self._normalize_score(imm_score, 'immunogenicity_score')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['score'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Conservancy' and 'conservancy_score' in category_metrics:
                                cons_score = category_metrics['conservancy_score']
                                score = self._normalize_score(cons_score, 'conservancy_score')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['conservancy_score'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Epitope' and 'score' in category_metrics:
                                epi_score = category_metrics['score']
                                score = self._normalize_score(epi_score, 'epitope_score')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['score'] = score
                                category_scores[category_name] = score
                                break
                            elif category_name == 'Developability' and 'developability_score' in category_metrics:
                                dev_score = category_metrics['developability_score']
                                score = self._normalize_score(dev_score, 'developability_score')
                                score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                                metric_scores[category_name]['developability_score'] = score
                                category_scores[category_name] = score
                                break
                                
                            # General case - process each metric
                            score = self._normalize_score(value, metric_name)
                            score = max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
                            metric_scores[category_name][metric_name] = score
                        except Exception as e:
                            logging.warning(f"Error scoring metric {metric_name}: {str(e)}")
                            metric_scores[category_name][metric_name] = 0.0  # Default score on error
                elif isinstance(category_metrics, list):
                    # For list metrics, assign score based on presence/length
                    metric_scores[category_name]["presence"] = 1.0 if category_metrics else 0.0
                    metric_scores[category_name]["length"] = min(1.0, len(category_metrics) / 10.0)  # Normalize based on length
                else:
                    # For other types, just assign a simple presence score
                    metric_scores[category_name]["value"] = 1.0 if category_metrics else 0.0
                
                # Calculate category score as average of metric scores if not already set
                if category_name not in category_scores and metric_scores[category_name]:
                    category_scores[category_name] = min(1.0, sum(metric_scores[category_name].values()) / len(metric_scores[category_name]))
            except Exception as e:
                logging.warning(f"Error processing category {category.value}: {str(e)}")
                category_scores[category.value] = 0.0  # Default score on error
        
        # Calculate overall score as weighted sum of category scores
        overall_score = 0.0
        
        for category, score in category_scores.items():
            try:
                weight = self.config.category_weights.get(MetricCategory(category), 0.0)
                # print(f"Category: {category}, Score: {score}, Weight: {weight}")
                overall_score += score * weight
                category_scores[category] = score * weight
            except Exception as e:
                logging.warning(f"Error applying weight for {category}: {str(e)}")
    
        # Ensure overall score is between 0 and 1
        overall_score = round(max(0.0, min(1.0, overall_score)), 4)
        
        scores = {
            'overall_score': overall_score,
            'category_scores': category_scores,
            'metric_scores': metric_scores
        }
        self._validate_scores(scores)
        return scores

    def rank_proteins(self, metrics_list: Union[List[ProteinMetrics], ProteinMetrics]) -> pd.DataFrame:
        """
        Rank proteins based on their metrics and generate individual reports.
        """
        if not metrics_list:
            logging.warning("⚠️ No proteins to rank!")
            return pd.DataFrame()

        logging.info(f"🎯 Ranking {len(metrics_list)} proteins...")
        
        # Process each protein
        all_proteins = []
        protein_scores = {}  # Store scores for later use
        failed_proteins = []  # Track proteins that failed during ranking
        
        for i, metrics in enumerate(tqdm(metrics_list, desc="📊 Calculating scores", unit="ab")):
            try:
                # Calculate scores
                scores = self.calculate_overall_score(metrics)
                protein_scores[metrics.sequence] = scores  # Store scores by sequence
                
                # Create row with all required columns
                row = {
                    'sequence': metrics.sequence,
                    'antigen': metrics.antigen,
                    'molecular_formula': metrics.molecular_formula,
                    'molecular_weight': metrics.molecular_weight,
                    'total_score': scores['overall_score'],
                    
                    #  Individual Metric Scores, Raw Category Scores, Weighted Category Scores
                    'gravy_score': metrics.protparam.get('gravy', 0.0),
                    'raw_gravy_score': scores['metric_scores'].get('ProtParam', {}).get('gravy', 0.0),
                    'raw_protparam_score': scores['category_scores'].get('ProtParam', 0.0),
                    'weighted_protparam_score': scores['category_scores'].get('ProtParam', 0.0) * self.config.category_weights[MetricCategory.PROTPARAM],
                    
                    'solubility': metrics.protparam.get('predicted_solubility', 'N/A'),
                    'raw_solubility_score': scores['metric_scores'].get('ProtParam', {}).get('solubility', 0.0),
                    
                    'immunogenicity_score': metrics.immunogenicity.get('immunogenic_score', 0.0),
                    'raw_immunogenicity_score': scores['category_scores'].get('Immunogenicity', 0.0),
                    'weighted_immunogenicity_score': scores['category_scores'].get('Immunogenicity', 0.0) * self.config.category_weights[MetricCategory.IMMUNOGENICITY],
                    
                    'melting_temperature': metrics.stability.get('melting_temperature', 0.0),
                    'raw_stability_score': scores['category_scores'].get('Stability', 0.0),
                    'weighted_stability_score': scores['category_scores'].get('Stability', 0.0) * self.config.category_weights[MetricCategory.STABILITY],
                    
                    'aggregation_propensity': metrics.aggregation.get('aggregation_propensity', 'N/A'),
                    'raw_aggregation_propensity_score': scores['metric_scores'].get('Aggregation', {}).get('aggregation_propensity', 0.0),
                    'aggregation_regions': len(metrics.aggregation.get('aggregation_prone_regions', [])),
                    'raw_aggregation_regions_score': scores['metric_scores'].get('Aggregation', {}).get('aggregation_regions', 0.0),
                    'raw_aggregation_score': scores['category_scores'].get('Aggregation', 0.0),
                    'weighted_aggregation_score': scores['category_scores'].get('Aggregation', 0.0) * self.config.category_weights[MetricCategory.AGGREGATION],
                    
                    'n_glyc_sites': metrics.glycosylation.get('n_glyc_sites_count', 0),
                    'raw_glycosylation_score': scores['category_scores'].get('Glycosylation', 0.0),
                    'weighted_glycosylation_score': scores['category_scores'].get('Glycosylation', 0.0) * self.config.category_weights[MetricCategory.GLYCOSYLATION],
                    
                    'gmqe_score': metrics.structure.get('gmqe_score', 0.0),
                    'raw_structure_score': scores['category_scores'].get('Structure', 0.0),
                    'weighted_structure_score': scores['category_scores'].get('Structure', 0.0) * self.config.category_weights[MetricCategory.STRUCTURE],
                    
                    'binding_affinity_kd': metrics.binding_affinity.get('dissociation_constant', 'N/A'),
                    'raw_binding_affinity_score': scores['category_scores'].get('Binding Affinity', 0.0),
                    'weighted_binding_affinity_score': scores['category_scores'].get('Binding Affinity', 0.0) * self.config.category_weights[MetricCategory.BINDING_AFFINITY],
                    
                    'epitope_score': metrics.epitope.get('score', 0.0),
                    'raw_epitope_score': scores['category_scores'].get('Epitope', 0.0),
                    'weighted_epitope_score': scores['category_scores'].get('Epitope', 0.0) * self.config.category_weights[MetricCategory.EPITOPE],
                    
                    'conservancy_score': metrics.conservancy.get('conservancy_score', 0.0),
                    'raw_conservancy_score': scores['category_scores'].get('Conservancy', 0.0),
                    'weighted_conservancy_score': scores['category_scores'].get('Conservancy', 0.0) * self.config.category_weights[MetricCategory.CONSERVANCY],
                    
                    'developability_score': metrics.developability.get('developability_score', 0.0),
                    'raw_developability_score': scores['category_scores'].get('Developability', 0.0),
                    'weighted_developability_score': scores['category_scores'].get('Developability', 0.0) * self.config.category_weights[MetricCategory.DEVELOPABILITY],
                }
                
                all_proteins.append(row)
                logging.debug(f"✅ Successfully ranked protein #{i+1}: {metrics.sequence[:20]}...")
                
            except Exception as e:
                error_msg = f"Error calculating scores: {str(e)}"
                logging.error(f"❌ Failed to rank protein #{i+1} ({metrics.sequence[:20]}...): {error_msg}")
                failed_proteins.append({
                    'sequence': metrics.sequence[:50] + '...' if len(metrics.sequence) > 50 else metrics.sequence,
                    'error': error_msg
                })
        
        # Log failed proteins during ranking
        if failed_proteins:
            logging.warning(f"⚠️ {len(failed_proteins)} proteins failed during ranking")
            if self.reports_dir:
                failed_ranking_path = os.path.join(self.reports_dir, "failed_during_ranking.json")
                try:
                    import json
                    with open(failed_ranking_path, 'w') as f:
                        json.dump(failed_proteins, f, indent=2)
                    logging.info(f"❌ Failed ranking details saved to {failed_ranking_path}")
                except Exception as e:
                    logging.error(f"❌ Could not save failed ranking details: {str(e)}")
        
        # Create DataFrame with results
        if not all_proteins:
            logging.warning("⚠️ No proteins could be successfully ranked!")
            return pd.DataFrame()
            
        self.df = pd.DataFrame(all_proteins)
        
        # Sort by total score
        self.df.sort_values('total_score', ascending=False, inplace=True)
        self.df.reset_index(drop=True, inplace=True)
        
        # Generate timestamp for file naming
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if self.reports_dir:
            output_dir = Path(self.reports_dir)
            logging.info("💾 Saving ranking results...")
            
            # Generate individual enhanced reports with correct ranks after sorting
            report_success_count = 0
            report_fail_count = 0
            
            for idx, row in enumerate(self.df.itertuples(), 1):
                try:
                    # Find the matching metrics object
                    matching_metrics = None
                    for m in metrics_list:
                        if m.sequence == row.sequence:
                            matching_metrics = m
                            break
                    
                    if matching_metrics:
                        scores = protein_scores.get(row.sequence)
                        if scores:
                            logging.info(f"📝 Generating enhanced report for protein {matching_metrics.sequence[:20]}... (Rank {idx})")
                            report_paths_dict = self.generate_protein_report(matching_metrics, scores, str(self.reports_dir), idx)
                            if report_paths_dict and any(report_paths_dict.values()): # Check if any report was generated
                                report_success_count += 1
                            else:
                                report_fail_count += 1
                        else:
                            logging.warning(f"⚠️ No scores found for protein {row.sequence[:20]}...")
                            report_fail_count += 1
                    else:
                        logging.warning(f"⚠️ No matching metrics found for protein in row {idx}")
                        report_fail_count += 1
                except Exception as e:
                    logging.error(f"Error generating report for protein {idx}: {str(e)}")
                    report_fail_count += 1
            
            logging.info(f"❌ Report generation: {report_success_count} successful, {report_fail_count} failed")
            
            # Save rankings CSV
            try:
                csv_file = self.save_rankings(output_dir, timestamp)
                if csv_file:
                    logging.info(f"💾 Rankings saved to CSV: {csv_file}")
                    
                    # Generate ranking report
                    try:
                        pdf_file = self.generate_ranking_report(str(csv_file), str(output_dir), timestamp)
                        if pdf_file:
                            logging.info(f"📂 Ranking report generated: {pdf_file}")
                        else:
                            logging.warning("⚠️ Could not generate ranking report PDF")
                    except Exception as e:
                        logging.error(f"❌ Error generating ranking report: {str(e)}")
                else:
                    logging.warning("⚠️ Could not save rankings to CSV")
            except Exception as e:
                logging.error(f"Error saving rankings to CSV: {str(e)}")
            
            logging.info(f"✅ Ranking complete! Results saved to {output_dir}")
        
        # Return a simplified version of the DataFrame for display
        display_df = self.df.copy()
        
        # Limit sequence length for display to prevent overflow
        if 'sequence' in display_df.columns:
            display_df['sequence'] = display_df['sequence'].apply(lambda s: s[:50] + '...' if len(s) > 50 else s)
        
        return display_df

    def get_ranking_results_as_dict(self) -> Optional[List[Dict[str, Any]]]:
        """Return the ranking results as a list of dictionaries.

        If rankings have been calculated (self.df exists), this method converts
        the DataFrame to a list of dictionaries, where each dictionary 
        represents an protein and its scores/metrics.

        Returns:
            Optional[List[Dict[str, Any]]]: A list of dictionaries representing 
                                            the ranked proteins, or None if 
                                            no ranking has been performed yet.
        """
        if self.df is not None:
            try:
                # Ensure all data is serializable for broader compatibility, though to_dict handles most.
                # Pandas to_dict(orient='records') handles numpy types well generally.
                return self.df.to_dict(orient='records')
            except Exception as e:
                logging.error(f"Error converting ranking DataFrame to dict: {e}")
                return None
        return None

    def save_rankings(self, output_dir: Path, timestamp: str):
        """Save rankings to CSV file"""
        if self.df is not None:
            # Create rankings directory
            rankings_dir = output_dir / 'rankings'
            rankings_dir.mkdir(parents=True, exist_ok=True)
            
            # Convert complex objects to string representation for CSV
            df_for_csv = self.df.copy()
            for col in df_for_csv.columns:
                # Check if the column contains complex objects
                if df_for_csv[col].apply(lambda x: isinstance(x, (dict, list))).any():
                    # Convert complex objects to string
                    df_for_csv[col] = df_for_csv[col].apply(lambda x: str(x) if isinstance(x, (dict, list)) else x)
            
            # Save rankings with timestamp
            rankings_file = rankings_dir / f'rankings_{timestamp}.csv'
            df_for_csv.to_csv(rankings_file, index=False)
            return rankings_file
        return None

    def generate_basic_ranking_report(self, csv_file: str, output_path: str):
        """
        Generate a basic PDF ranking report without relying on protein_report_generator.
        This is a fallback method used when the import fails.
        """
        try:
            import pandas as pd
            from fpdf import FPDF
            from datetime import datetime
            
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Create PDF
            pdf = FPDF()
            pdf.add_page()
            
            # Add title
            pdf.set_font('Arial', 'B', 16)
            pdf.cell(0, 10, 'Antibody Ranking Report', 0, 1, 'C')
            pdf.cell(0, 5, f'Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', 0, 1, 'C')
            pdf.ln(10)
            
            # Add summary
            pdf.set_font('Arial', 'B', 14)
            pdf.cell(0, 10, 'Summary', 0, 1, 'L')
            pdf.set_font('Arial', '', 10)
            pdf.cell(0, 7, f'Number of proteins: {len(df)}', 0, 1, 'L')
            
            if len(df) > 0:
                pdf.cell(0, 7, f'Average score: {df["total_score"].mean():.3f}', 0, 1, 'L')
                pdf.cell(0, 7, f'Score range: {df["total_score"].min():.3f} - {df["total_score"].max():.3f}', 0, 1, 'L')
            pdf.ln(10)
            
            # Add ranking table - Top 10 only
            pdf.set_font('Arial', 'B', 14)
            pdf.cell(0, 10, 'Top Ranked proteins', 0, 1, 'L')
            
            # Define columns to show
            display_cols = ['total_score', 'molecular_formula', 'molecular_weight']
            
            # Add header
            pdf.set_font('Arial', 'B', 10)
            col_width = 50
            pdf.cell(20, 7, 'Rank', 1, 0, 'C')
            pdf.cell(30, 7, 'Score', 1, 0, 'C')
            pdf.cell(70, 7, 'Formula', 1, 0, 'C')
            pdf.cell(30, 7, 'Weight', 1, 0, 'C')
            pdf.cell(40, 7, 'Antigen', 1, 1, 'C')
            
            # Add data rows - top 10 only
            pdf.set_font('Arial', '', 10)
            for i, row in df.head(10).iterrows():
                pdf.cell(20, 7, str(i+1), 1, 0, 'C')
                pdf.cell(30, 7, f"{row['total_score']:.3f}", 1, 0, 'C')
                formula = str(row['molecular_formula'])
                if len(formula) > 15:
                    formula = formula[:12] + '...'
                pdf.cell(70, 7, formula, 1, 0, 'C')
                pdf.cell(30, 7, f"{row['molecular_weight']:.1f}", 1, 0, 'C')
                
                antigen = str(row.get('antigen', 'Unknown'))
                if len(antigen) > 20:
                    antigen = antigen[:17] + '...'
                pdf.cell(40, 7, antigen, 1, 1, 'C')
            
            # Add note about CSV file
            pdf.ln(10)
            pdf.set_font('Arial', '', 10)
            pdf.cell(0, 7, f'For complete data, refer to: {os.path.basename(csv_file)}', 0, 1, 'L')
            
            # Save the PDF
            pdf.output(output_path)
            logging.info(f"Generated basic ranking report: {output_path}")
            return True
        except Exception as e:
            logging.error(f"Failed to generate basic ranking report: {str(e)}")
            return False
            
    def generate_ranking_report(self, csv_file: str, output_dir: str, timestamp: str):
        """Generate ranking report for all proteins"""
        try:
            # Create rankings directory
            rankings_dir = Path(output_dir) / "rankings"
            rankings_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate PDF report
            pdf_file = rankings_dir / f"ranking_report_{timestamp}.pdf"
            logging.info(f"📊 Generating comprehensive ranking report: {pdf_file}")
            
            # First try using the imported function if available
            if ranking_report_available:
                try:
                    gen_ranking_report(csv_file, str(pdf_file), top_n=None)
                    logging.info(f"Ranking report generated successfully: {pdf_file}")
                    return pdf_file
                except Exception as e:
                    logging.warning(f"Failed to use standard report generator: {str(e)}")
                    logging.warning("Falling back to basic report generator")
            else:
                logging.warning("Ranking report generator not available. Using basic report generator.")
            
            # Fall back to basic report
            if self.generate_basic_ranking_report(csv_file, str(pdf_file)):
                logging.info(f"Basic ranking report generated successfully: {pdf_file}")
                return pdf_file
            else:
                return None
                
        except Exception as e:
            import traceback
            logging.error(f"Failed to generate ranking report: {str(e)}")
            logging.debug(traceback.format_exc())
            return None

    def _safe_serialize(self, obj):
        """Convert complex objects to serializable format for DataFrame storage"""
        if isinstance(obj, list):
            # Handle list of tuples/objects
            if obj and isinstance(obj[0], tuple):
                return [list(item) for item in obj]
            return obj
        elif isinstance(obj, dict):
            # Check for nested types that need conversion
            result = {}
            for k, v in obj.items():
                if isinstance(v, list) and v and isinstance(v[0], tuple):
                    result[k] = [list(item) for item in v]
                elif isinstance(v, (dict, list)):
                    result[k] = self._safe_serialize(v)
                else:
                    result[k] = v
            return result
        else:
            return obj

def rank_proteins_from_metrics(
    protein_sequences: str,
    output_dir: str,
    config: Optional[ScoringConfig] = None,
    generate_pdf: bool = True,  # Add flags here for the standalone function
    generate_csv: bool = True   # Add flags here for the standalone function
) -> pd.DataFrame:
    """
    Main function to rank proteins from metrics file.
    
    Args:
        protein_sequences: Path to file containing protein sequences
        output_dir: Directory to save ranking reports
        config: Optional custom scoring configuration
        generate_pdf: Whether to generate PDF reports for individual proteins.
        generate_csv: Whether to generate CSV reports for individual proteins.
        
    Returns:
        DataFrame with ranked proteins
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up logging
    log_file = os.path.join(output_dir, "ranking.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file)
        ]
    )
    
    # Log information about available modules
    logging.info(f"Enhanced report generator available: {enhanced_report_available}")
    logging.info(f"Ranking report generator available: {ranking_report_available}")
    
    # Try to add bio_v010 to path if needed
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    bio_v010_dir = os.path.join(parent_dir, 'bio_v010')
    
    if os.path.exists(bio_v010_dir) and bio_v010_dir not in sys.path:
        sys.path.append(bio_v010_dir)
        logging.info(f"Added {bio_v010_dir} to Python path")
    
    # Initialize validator and ranker
    ranker = ProteinRanker(config, generate_pdf=generate_pdf, generate_csv=generate_csv)
    ranker.set_output_directory(output_dir)
    validator = ProteinValidator(pdb_files_path=f"./{output_dir}/pdbs")
    
    # Process proteins and handle failures gracefully
    results = validator.process_proteins(protein_sequences, output_dir)
    metrics_list = validator.get_successful_metrics(results)
    
    # Track successful and failed proteins
    processed_proteins = len(results)
    successful_proteins = len(metrics_list)
    failed_proteins = processed_proteins - successful_proteins
    
    # Create a file with information about failed proteins
    failed_proteins_path = os.path.join(output_dir, "failed_proteins.txt")
    with open(failed_proteins_path, 'w') as f:
        f.write("======== Failed proteins ========\n")
        f.write(f"Total processed: {processed_proteins}\n")
        f.write(f"Failed: {failed_proteins}\n\n")
        f.write("Details of failed proteins:\n")
        f.write("--------------------------------\n\n")
        
        for result in results:
            if not result.success:
                f.write(f"Sequence: {result.sequence[:50]}...\n")
                f.write(f"Error: {result.error}\n\n")
    
    logging.info(f"Total proteins processed: {processed_proteins}")
    logging.info(f"Successful proteins: {successful_proteins}")
    logging.info(f"Failed proteins: {failed_proteins}")
    logging.info(f"Failed proteins information saved to: {failed_proteins_path}")
    
    if not metrics_list:
        logging.warning("No valid proteins found. Cannot proceed with ranking.")
        # Create an empty report to indicate processing completed but no valid proteins were found
        empty_report_path = os.path.join(output_dir, "no_valid_proteins.txt")
        with open(empty_report_path, 'w') as f:
            f.write("No valid proteins found to rank.\n")
            f.write(f"Total processed: {processed_proteins}\n")
            f.write(f"All processing failed. See {failed_proteins_path} for details.\n")
        
        return pd.DataFrame()  # Return empty DataFrame since no valid proteins
    
    # Rank proteins and generate reports
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    try:
        # Try to rank proteins - even if there's just one
        rankings = ranker.rank_proteins(metrics_list)
        
        # Save and report results
        csv_file = ranker.save_rankings(Path(output_dir), timestamp)
        pdf_file = ranker.generate_ranking_report(str(csv_file), output_dir, timestamp)
        
        # Generate a simple text report that includes info about failed proteins
        txt_report_path = os.path.join(output_dir, f"simple_report_{timestamp}.txt")
        with open(txt_report_path, 'w') as f:
            f.write("======================================\n")
            f.write("       ANTIBODY RANKING RESULTS       \n")
            f.write("======================================\n\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Number of proteins ranked: {len(rankings)}\n")
            f.write(f"Number of proteins that failed: {failed_proteins}\n\n")
            
            if len(rankings) > 0:
                f.write(f"Average score: {rankings['total_score'].mean():.3f}\n")
                f.write(f"Score range: {rankings['total_score'].min():.3f} - {rankings['total_score'].max():.3f}\n\n")
                
                f.write("TOP 5 proteins (or all if less than 5):\n")
                f.write("----------------\n")
                for i, row in rankings.head(min(5, len(rankings))).iterrows():
                    f.write(f"{i+1}. Score: {row['total_score']:.3f}, Formula: {row['molecular_formula']}\n")
                    f.write(f"   Antigen: {row.get('antigen', 'Unknown')}\n")
                    f.write(f"   Weight: {row['molecular_weight']:.1f}\n")
                    f.write("\n")
            
            f.write("\nFULL RESULTS:\n")
            f.write(f"- CSV file: {csv_file}\n")
            f.write(f"- Ranking report: {pdf_file}\n")
            f.write(f"- Log file: {log_file}\n")
            f.write(f"- Failed proteins: {failed_proteins_path}\n")
        
        logging.info(f"\nRanking results:")
        logging.info(f"- Ranked proteins: {len(rankings)}")
        logging.info(f"- Failed proteins: {failed_proteins}")
        logging.info(f"- CSV file: {csv_file}")
        logging.info(f"- Ranking report: {pdf_file}")
        logging.info(f"- Simple text report: {txt_report_path}")
        logging.info(f"- Individual reports: {output_dir}/protein_reports/")
        
        return rankings
        
    except Exception as e:
        logging.error(f"Error during ranking process: {str(e)}")
        
        # Create a fallback report with basic information
        fallback_report_path = os.path.join(output_dir, f"fallback_report_{timestamp}.txt")
        with open(fallback_report_path, 'w') as f:
            f.write("======================================\n")
            f.write("       ERROR IN RANKING PROCESS       \n")
            f.write("======================================\n\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"An error occurred during the ranking process: {str(e)}\n\n")
            f.write(f"Number of successful proteins that were processed: {successful_proteins}\n")
            f.write(f"Number of proteins that failed: {failed_proteins}\n\n")
            f.write(f"Failed proteins list: {failed_proteins_path}\n")
        
        logging.info(f"Fallback report generated: {fallback_report_path}")
        
        # Create a basic DataFrame with minimal information
        fallback_df = pd.DataFrame([{
            'sequence': m.sequence,
            'molecular_formula': m.molecular_formula,
            'antigen': m.antigen,
            'molecular_weight': m.molecular_weight
        } for m in metrics_list])
        
        return fallback_df

if __name__ == "__main__":
    # Example usage
    protein_sequences = "example_proteins.txt"
    output_dir = "ranking_results"
    # Pass flags to the main function
    rankings = rank_proteins_from_metrics(protein_sequences, output_dir, generate_pdf=True, generate_csv=True)
    print(f"Ranked {len(rankings)} proteins successfully")
