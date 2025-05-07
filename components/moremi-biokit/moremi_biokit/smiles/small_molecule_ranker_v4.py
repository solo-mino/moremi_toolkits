"""
Small Molecule Ranker V4

This module implements a comprehensive scoring system for chemical compounds based on
their physicochemical and ADMET properties. It provides normalized scoring (0-1) for
each property and combines them with equal weighting for final ranking.

Key Features:
1. Equal weighting for all metrics
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
from .small_molecule_validator_v3 import (
    MoleculeMetrics, MetricCategory, MetricRanges, SmallMoleculeValidator
)
from pathlib import Path
from .reports import generate_enhanced_report
import logging
from tqdm import tqdm

@dataclass
class ScoringConfig:
    """Configuration for property scoring calculations"""
    
    # Equal weights for all categories
    category_weights = {cat: 1.0 for cat in MetricCategory}
    
    # Property scoring functions and ranges
    property_configs = {
        # Physicochemical
        'tpsa': {
            'optimal_range': (0, 140),
            'score_func': lambda x: np.where(
                (x >= 0) & (x <= 140),
                1.0,
                np.divide(140, x, out=np.zeros_like(x), where=x != 0)  # Avoid ZeroDivisionError
            )
        },

        'logp': {
            'optimal_range': (1, 3),
            'score_func': lambda x: 1 - np.minimum(1, np.abs(x - 2) / 3)
        },
        
        # Druglikeness
        'lipinski_violations': {
            'score_func': lambda x: 1 - (len(x) if isinstance(x, list) else int(x)) / 4  # 4 is max violations
        },
        'qed': {
            'score_func': lambda x: np.where(x > 0.67, 1.0, x / 0.67)
        },
        
        'bioavailability': {
            'threshold': 0.7,
            'score_func': lambda x: np.minimum(1, x / 0.7)  # Values greater than 0.7 are better
        },
        
        # Medicinal Chemistry
        'synthetic_accessibility': {
            'optimal_range': (1, 5),
            'score_func': lambda x: 1 - np.minimum(1, (x - 1) / 9)  # Scale 1-10
        },
        
        # Toxicity
        'herg': {
            'optimal_range': (0, 0.5),
            'score_func': lambda x: 1 - x  # Lower is better
        },
        
        # Absorption
        'caco2': {
            'threshold': -4.7,
            'score_func': lambda x: 1 / (1 + np.exp(-0.5 * (x + 4.7))),  # Sigmoid around threshold -4.7
            'reference': 'Pires et al., J Med Chem, 2015, DOI:10.1021/acs.jmedchem.5b00104; Kus et al., ADMET & DMPK, 2023 :cite[3]:cite[9]'
        },
        'pampa': {
            'threshold': 0.7,  # Probability >0.7 = high confidence
            'score_func': lambda x: max(0.0, min(1.0, (x - 0.3) / 0.4)),  # Linear scaling: 0.3→0, 0.7→1
            'reference': 'Di et al., Eur J Pharm Sci, 2003 (PAMPA benchmarks)'
        },
        'mdck': {
            'threshold': -5.0,  # logPapp > -5.0 = high permeability
            'score_func': lambda x: 1 / (1 + np.exp(-3 * (x + 5.0))),  # Sigmoid centered at -5.0
            'reference': 'Pires et al., J Med Chem (2015) (DOI:10.1021/acs.jmedchem.5b00104).'
        },
        'hia': {
            'threshold': 0.7,
            'score_func': lambda x: x
        },
        'pgp_substrate': {
            'optimal_range': (0, 0.5),
            'score_func': lambda x: 1 - x
        },
        
        # PGP Inhibitor: Use sigmoid for smoother penalties
        'pgp_inhibitor': {
            'threshold': 0.4,
            'score_func': lambda x: 1 / (1 + np.exp(-5 * (x - 0.4)))  # Steep penalty outside 0.3–0.5
        },
        
        # Distribution
        # ==============================================================
        # DISTRIBUTION SCORING FUNCTIONS (Based on Table 6B and Section D) of the paper Design and In-Silico Evaluation of Pyridine-4-Carbohydrazide Derivatives for Potential Therapeutic Applications -> https://www.auctoresonline.org/uploads/articles/1742643482JSCRI-24-CR-235-Galley_Proof.pdf 
        # ==============================================================
        'ppb': {
            'threshold': 90,  # PPB (%) threshold for low free fraction
            'score_func': lambda x: max(0.0, min(1.0, 1 - (x / 100))),   # Inversely scaled to PPB
            'reference': 'pkCSM Theory, Biosig Lab, 2023 :cite[3]:cite[9]'
        },
        
        # ----------------------------------------------------------
        # Blood-Brain Barrier (BBB) Permeability - LogBB
        # Thresholds: LogBB > 0.3 (good), LogBB < -1.0 (poor)
        # Sigmoid rewards values > 0.3 (Pires et al., 2015)
        # ----------------------------------------------------------
        'bbb': {
            'threshold': 0.3,
            'score_func': lambda x: 1 / (1 + np.exp(-4 * (x - 0.3))),  # Steep penalty below 0.3
            'reference': 'Pires et al., J Med Chem, 2015, DOI:10.1021/acs.jmedchem.5b00104'
        },
        
        # ----------------------------------------------------------
        # Fraction Unbound (FU) - Fraction (0-1)
        # Higher FU = Better bioavailability (Watanabe et al., 2018)
        # ----------------------------------------------------------
        'fu': {
            'threshold': 30,  # Target FU ≥ 0.3 for therapeutic activity
            'score_func': lambda x: max(0.0, min(1.0, (x/100) / 0.3)),  # Linear scaling capped at 1.0
            'reference': 'Watanabe et al., Drug Metab Dispos, 2018, DOI:10.1124/dmd.117.079020'
    },
        
        # ==============================================================
        # VOLUME OF DISTRIBUTION (VDss) SCORING FUNCTION
        # --------------------------------------------------------------
        # Criteria (Pires et al., 2015):
        # - Good distribution: VDss > 2.81 L/kg (logVDss > 0.45)
        # - Poor distribution: VDss < 0.71 L/kg (logVDss < -0.15)
        # Scoring: Sigmoid centered at 3 L/kg, capped at 10 L/kg.
        # ==============================================================
        'vdss': {
            'threshold': 3.0,  # Center of sigmoid (L/kg)
            'score_func': lambda x: min(1.0, 1 / (1 + np.exp(-1.5 * (x - 3.0)))),  # Steepness tuned for 3 L/kg
            'reference': 'Pires et al., J Med Chem, 2015, DOI:10.1021/acs.jmedchem.5b00104'
        },
        
        # Metabolism
        'cyp2c9_inhibition': {
            'score_func': lambda x: 1 - x  # Lower is better
        },
        'cyp2d6_inhibition': {
            'score_func': lambda x: 1 - x
        },
        'cyp3a4_inhibition': {
            'score_func': lambda x: 1 - x
        },
        
        # Excretion
        'half_life': {
            'threshold': 4,
            'score_func': lambda x: 1 / (1 + np.exp(-0.5 * (x - 4)))
        },
        'clearance': {
            'threshold': 5,
            'score_func': lambda x: 1 / (1 + np.exp(-0.5 * (x - 5)))
        }
    }

class SmallMoleculeRankerV4:
    """Enhanced molecule ranker implementing comprehensive scoring system"""
    
    def __init__(self, config: Optional[ScoringConfig] = None, generate_pdf: Optional[bool] = False, generate_csv: Optional[bool] = True):
        """Initialize ranker with configuration"""
        self.config = config or ScoringConfig()
        self.ranges = MetricRanges()
        self.df = None
        # Add output directory for reports
        self.reports_dir = None
        self.generate_pdf = generate_pdf
        self.generate_csv = generate_csv

    def set_output_directory(self, output_dir: str):
        """Set the output directory for reports"""
        self.reports_dir = Path(output_dir)
        # Create directories for individual molecule reports and rankings
        (self.reports_dir / "molecule_reports").mkdir(parents=True, exist_ok=True)
        (self.reports_dir / "rankings").mkdir(parents=True, exist_ok=True)

    def generate_molecule_report(self, metrics: MoleculeMetrics, scores: Dict, output_dir: str, rank: int = 0):
        """Generate report for a single molecule"""
        try:
            if not metrics or not metrics.smiles:
                raise ValueError("Invalid metrics object or missing SMILES")

            # Import required property calculators using relative paths
            from .property_calculators import (
                PhysiochemicalProperties,
                SolubilityProperties,
                PharmacokineticsProperties,
                DruglikenessProperties,
                MedicinalChemistryProperties,
                ADMETPredictor
            )
            
            # Calculate additional properties with error handling
            try:
                phys_props = PhysiochemicalProperties(metrics.smiles).calculate_physiochemical_properties()
            except Exception as e:
                print(f"Warning: Error calculating physiochemical properties - {str(e)}")
                phys_props = {}

            try:
                sol_props = SolubilityProperties(metrics.smiles).calculate_solubility()
            except Exception as e:
                print(f"Warning: Error calculating solubility properties - {str(e)}")
                sol_props = {}

            try:
                pharm_props = PharmacokineticsProperties(metrics.smiles).calculate_pharmacokinetics()
            except Exception as e:
                print(f"Warning: Error calculating pharmacokinetics properties - {str(e)}")
                pharm_props = {}

            try:
                drug_props = DruglikenessProperties(metrics.smiles).calculate_druglikeness()
            except Exception as e:
                print(f"Warning: Error calculating druglikeness properties - {str(e)}")
                drug_props = {}

            try:
                med_props = MedicinalChemistryProperties(metrics.smiles).calculate_medicinal_chemistry()
            except Exception as e:
                print(f"Warning: Error calculating medicinal chemistry properties - {str(e)}")
                med_props = {}

            try:
                admet_props = ADMETPredictor(metrics.smiles).calculate_admet_properties()
            except Exception as e:
                print(f"Warning: Error calculating ADMET properties - {str(e)}")
                admet_props = {}
            
            # Prepare data for report with all property groups
            molecule_data = {
                "smiles": metrics.smiles,
                "molecular_formula": metrics.molecular_formula,
                "molecular_weight": metrics.molecular_weight,
                "rank": rank,  # Use the provided rank parameter
                "physicochemical": {
                    **(metrics.physicochemical or {}),
                    **(phys_props or {}),
                    'molecular_density': admet_props['molecular_density'],
                    'molecular_volume': admet_props['molecular_volume'],
                    'ring_count': admet_props['ring_count'],
                    'max_ring_size': admet_props['max_ring_size'],
                    'molecular_flexibility': admet_props['molecular_flexibility'],
                    'molecular_rigidity': admet_props['molecular_rigidity'],
                    'heteroatom_count': admet_props['heteroatom_count'],
                    'formal_charge': admet_props['formal_charge'],
                    'stereogenic_centers': admet_props['stereogenic_centers'],
                },
                "medicinal_chemistry": {
                    **(metrics.medicinal_chemistry or {}),
                    **(med_props or {})
                },
                "solubility": {
                    **(sol_props or {})
                },
                "lipophilicity": metrics.lipophilicity or {},
                "druglikeness": {
                    **(metrics.druglikeness or {}),
                    **(drug_props or {}),
                     
                },
                "absorption": {
                    **(metrics.absorption or {}),
                    "oral_bioavailability": admet_props['oral_bioavailability'],
                    "gi_absorption": pharm_props['gi_absorption'],
                    "log_kp": pharm_props['log_kp']
                },
                "distribution": {
                    **(metrics.distribution or {}),
                    'hydration_free_energy': admet_props['hydration_free_energy'],
                    'aqueous_solubility': admet_props['aqueous_solubility'],
                },
                "metabolism": {
                    **(metrics.metabolism or {}),
                    "hepatic_clearance": admet_props['hepatic_clearance'],
                    "microsomal_clearance": admet_props['microsomal_clearance'],
                    "plasma_clearance": admet_props['plasma_clearance'],
                    "cyp1a2_inhibition": admet_props['cyp1a2_inhibition'],
                    "cyp2c19_inhibition": admet_props['cyp2c19_inhibition'],
                    "cyp2c9_substrate": admet_props['cyp2c9_substrate'],
                    "cyp2d6_substrate": admet_props['cyp2d6_substrate'],
                    "cyp3a4_substrate": admet_props['cyp3a4_substrate']
                },
                "excretion": metrics.excretion or {},
                "toxicity": {
                    **(metrics.toxicity or {}),
                    "drug_induced_liver_injury": admet_props['drug_induced_liver_injury'],
                    "ames_mutagenicity": admet_props['ames_mutagenicity'],
                    "skin_sensitization": admet_props['skin_sensitization'],
                    "lethal_dose_50": admet_props['lethal_dose_50'],
                    "clinical_toxicity": admet_props['clinical_toxicity'],
                    # Nuclear Receptor Pathways
                    "androgen_receptor_binding": admet_props['androgen_receptor_binding'],
                    "androgen_receptor_lbd": admet_props['androgen_receptor_lbd'],
                    "aryl_hydrocarbon_receptor": admet_props['aryl_hydrocarbon_receptor'],
                    "aromatase_inhibition": admet_props['aromatase_inhibition'],
                    "estrogen_receptor_binding": admet_props['estrogen_receptor_binding'],
                    "estrogen_receptor_lbd": admet_props['estrogen_receptor_lbd'],
                    "ppar_gamma": admet_props['ppar_gamma'],
                    # Stress Response Pathways
                    "oxidative_stress_response": admet_props['oxidative_stress_response'],
                    "dna_damage_response": admet_props['dna_damage_response'],
                    "heat_shock_response": admet_props['heat_shock_response'],
                    "mitochondrial_toxicity": admet_props['mitochondrial_toxicity'],
                    "p53_pathway_activation": admet_props['p53_pathway_activation']
                },
                "scores": {
                    "overall_score": scores.get('overall_score', 0.0),
                    "category_scores": scores.get('category_scores', {}),
                    "metric_scores": scores.get('metric_scores', {})
                }
            }

            # Create the main output directories
            molecule_reports_dir = os.path.join(output_dir, "molecule_reports")
            rankings_dir = os.path.join(output_dir, "rankings")
            os.makedirs(molecule_reports_dir, exist_ok=True)
            os.makedirs(rankings_dir, exist_ok=True)

            # Generate molecule-specific reports with simplified names
            # pdf_path = os.path.join(molecule_reports_dir, f"{metrics.molecular_formula}_analysis.pdf")
            # csv_path = os.path.join(molecule_reports_dir, f"{metrics.molecular_formula}_analysis.csv")
            
            # Generate reports
            # TODO: ADDED FALSE HERE TO IGNORE PDF GENERTATION FOR THE SAKE OF TOXICITY PAPER - User wants to control this
            generate_enhanced_report(molecule_data, molecule_reports_dir, generate_pdf=self.generate_pdf, generate_csv=self.generate_csv)
            
            # Generate ranking reports with timestamps
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            # ranking_csv = os.path.join(rankings_dir, f"rankings_{timestamp}.csv")
            # ranking_pdf = os.path.join(rankings_dir, f"ranking_report_{timestamp}.pdf")
            
            # print(f"\n\n\nmolecular data ->: {molecule_data}\n\n\n")
            print(f"Generated reports for {metrics.molecular_formula} in {molecule_reports_dir}")
        except Exception as e:
            import traceback
            print(f"Warning: Failed to generate report for {metrics.molecular_formula}")
            print("Error details:")
            print(traceback.format_exc())

    # TODO: NOT USED YET - REMOVE
    def calculate_metric_score(self, value: float, metric_name: str) -> float:
        """Calculate score for a specific metric"""
        if metric_name not in self.config.property_configs:
            return 0.0
            
        config = self.config.property_configs[metric_name]
        return self._normalize_score(value, metric_name)

    def _normalize_score(self, value: Union[float, int, list], metric_name: str) -> float:
        """Normalize a metric value to a score between 0 and 1"""
        if metric_name not in self.config.property_configs:
            return float(value) if isinstance(value, (int, float)) else 0.0
            
        config = self.config.property_configs[metric_name]
        
        # Handle special case for lipinski violations
        if metric_name == 'lipinski_violations':
            if isinstance(value, list):
                num_violations = len(value)
            else:
                try:
                    num_violations = float(value)
                except (ValueError, TypeError):
                    num_violations = 0
            return max(0.0, min(1.0, 1 - (num_violations / 4)))  # 4 is max violations
            
        # Handle special case for bioavailability score
        if metric_name == 'bioavailability_score':
            try:
                score = float(value)
                return max(0.0, min(1.0, score))  # Ensure score is between 0 and 1
            except (ValueError, TypeError):
                return 0.0
            
        try:
            return config['score_func'](float(value))
        except (ValueError, TypeError):
            return 0.0
        
    def _calculate_category_score(self, scores: Dict[str, float]) -> float:
        """Calculate the average score for a category"""
        if not scores:
            return 0.0
        return np.mean(list(scores.values()))
        
    def _validate_scores(self, scores: Dict[str, Any]) -> None:
        """Validate that all scores are within the valid range [0,1]"""
        # Check overall score
        if not 0 <= scores['overall_score'] <= 1:
            print(f"WARNING: Overall score {scores['overall_score']} outside valid range [0,1]")
            scores['overall_score'] = max(0.0, min(1.0, scores['overall_score']))
        
        # Check category scores
        for category, score in scores['category_scores'].items():
            if not 0 <= score <= 1:
                print(f"WARNING: Category {category} score {score} outside valid range [0,1]")
                scores['category_scores'][category] = max(0.0, min(1.0, score))
        
        # Check metric scores
        for category, metrics in scores['metric_scores'].items():
            for metric, score in metrics.items():
                if not 0 <= score <= 1:
                    print(f"WARNING: Metric {metric} score {score} outside valid range [0,1]")
                    scores['metric_scores'][category][metric] = max(0.0, min(1.0, score))

    def calculate_overall_score(self, metrics: MoleculeMetrics) -> Dict[str, Any]:
        """Calculate overall score and category scores for a molecule"""
        try:
            print(f"\n🧮 Calculating scores for {metrics.molecular_formula}")
            logging.info(f"🧮 Calculating scores for {metrics.molecular_formula}")
            
            category_scores = {}
            metric_scores = {}
            
            # Calculate scores for each category
            for category in MetricCategory:
                category_name = category.value
                category_metrics = getattr(metrics, category_name.lower().replace(' ', '_').replace('-', ''))
                
                if category_metrics:
                    metric_scores[category_name] = {}
                    for metric_name, value in category_metrics.items():
                        # Special handling for druglikeness metrics
                        if category_name == "Druglikeness":
                            if metric_name == "lipinski_violations":
                                if isinstance(value, list):
                                    num_violations = len(value)
                                else:
                                    try:
                                        num_violations = float(value)
                                    except (ValueError, TypeError):
                                        num_violations = 0
                                score = max(0.0, min(1.0, 1 - (num_violations / 4)))
                                metric_scores[category_name][metric_name] = score
                                continue
                            elif metric_name == "bioavailability_score":
                                try:
                                    value = float(value)
                                    score = max(0.0, min(1.0, value))  # Direct normalization for bioavailability
                                    metric_scores[category_name][metric_name] = score
                                    continue
                                except (ValueError, TypeError):
                                    score = 0.0
                        
                        # Normal metric scoring
                        score = self._normalize_score(value, metric_name)
                        # Ensure score is between 0 and 1
                        score = max(0.0, min(1.0, score))
                        metric_scores[category_name][metric_name] = score
                    
                    # Calculate category score as average of metric scores
                    if metric_scores[category_name]:
                        # Ensure we don't exceed 1.0 for any category
                        category_scores[category_name] = min(1.0, sum(metric_scores[category_name].values()) / len(metric_scores[category_name]))
                    else:
                        category_scores[category_name] = 0.0
            
            # Calculate overall score as weighted average of category scores
            total_weight = sum(self.config.category_weights.values())
            overall_score = 0.0
            
            if total_weight > 0:
                weighted_sum = 0.0
                for category, score in category_scores.items():
                    weight = self.config.category_weights.get(MetricCategory(category), 1.0)
                    weighted_sum += score * weight
                overall_score = weighted_sum / total_weight
                
            # Ensure overall score is between 0 and 1
            overall_score = round(max(0.0, min(1.0, overall_score)), 4)
            
            scores = {
                'overall_score': overall_score,
                'category_scores': category_scores,
                'metric_scores': metric_scores
            }
            self._validate_scores(scores)
            return scores
        except Exception as e:
            print(f"Error calculating scores for {metrics.molecular_formula}: {str(e)}")
            logging.error(f"❌ Error calculating scores for {metrics.molecular_formula}: {str(e)}")
            return {
                'overall_score': 0.0,
                'category_scores': {},
                'metric_scores': {}
            }

    def rank_molecules(self, metrics_list: List[MoleculeMetrics]) -> pd.DataFrame:
        """
        Rank molecules based on their metrics and generate individual reports.
        """
        if not metrics_list:
            logging.warning("⚠️ No molecules to rank!")
            return pd.DataFrame()

        logging.info(f"🎯 Ranking {len(metrics_list)} molecules...")
        
        # Process each molecule
        all_molecules = []
        molecule_scores = {}  # Store scores for later use
        
        for metrics in tqdm(metrics_list, desc="📊 Calculating scores", unit="mol"):
            scores = self.calculate_overall_score(metrics)
            molecule_scores[metrics.smiles] = scores  # Store scores by SMILES
            
            # Create row with all required columns
            row = {
                'SMILES': metrics.smiles,
                'Molecular Formula': metrics.molecular_formula,
                'Molecular Weight': metrics.molecular_weight,
                'Overall Score': scores['overall_score'],
                
                # Category Scores
                'Physicochemical Score': scores['category_scores'].get('Physicochemical', 0.0),
                'Medicinal Chemistry Score': scores['category_scores'].get('Medicinal Chemistry', 0.0),
                'Druglikeness Score': scores['category_scores'].get('Druglikeness', 0.0),
                'Absorption Score': scores['category_scores'].get('Absorption', 0.0),
                'Distribution Score': scores['category_scores'].get('Distribution', 0.0),
                'Metabolism Score': scores['category_scores'].get('Metabolism', 0.0),
                'Excretion Score': scores['category_scores'].get('Excretion', 0.0),
                'Toxicity Score': scores['category_scores'].get('Toxicity', 0.0),
                
                # Individual Metric Scores and Values
                'tpsa score': scores['metric_scores'].get('Physicochemical', {}).get('tpsa', 0.0),
                'tpsa value': metrics.physicochemical.get('tpsa', 0.0),
                'logp score': scores['metric_scores'].get('Physicochemical', {}).get('logp', 0.0),
                'logp value': metrics.physicochemical.get('logp', 0.0),
                'qed score': scores['metric_scores'].get('Medicinal Chemistry', {}).get('qed', 0.0),
                'qed value': metrics.medicinal_chemistry.get('qed', 0.0),
                
                # Add all logP values
                'ilogp value': metrics.lipophilicity.get('ilogp', 0.0),
                'xlogp3 value': metrics.lipophilicity.get('xlogp3', 0.0),
                'wlogp value': metrics.lipophilicity.get('wlogp', 0.0),
                'mlogp value': metrics.lipophilicity.get('mlogp', 0.0),
                'silicos_it value': metrics.lipophilicity.get('silicos_it', 0.0),
                'consensus_logp value': metrics.lipophilicity.get('consensus_logp', 0.0),
                
                # Medicinal Chemistry
                'synthetic_accessibility score': scores['metric_scores'].get('Medicinal Chemistry', {}).get('synthetic_accessibility', 0.0),
                'synthetic_accessibility value': metrics.medicinal_chemistry.get('synthetic_accessibility', 0.0),
                
                # Drug-likeness
                'lipinski_violations score': scores['metric_scores'].get('Druglikeness', {}).get('lipinski_violations', 0.0),
                'lipinski_violations value': metrics.druglikeness.get('lipinski_violations', []),
                'bioavailability_score score': scores['metric_scores'].get('Druglikeness', {}).get('bioavailability_score', 0.0),
                'bioavailability_score value': metrics.druglikeness.get('bioavailability_score', 0.0),
                
                # Absorption
                'caco2 score': scores['metric_scores'].get('Absorption', {}).get('caco2', 0.0),
                'caco2 value': metrics.absorption.get('caco2', 0.0),
                'pampa score': scores['metric_scores'].get('Absorption', {}).get('pampa', 0.0),
                'pampa value': metrics.absorption.get('pampa', 0.0),
                'mdck score': scores['metric_scores'].get('Absorption', {}).get('mdck', 0.0),
                'mdck value': metrics.absorption.get('mdck', 0.0),
                'hia score': scores['metric_scores'].get('Absorption', {}).get('hia', 0.0),
                'hia value': metrics.absorption.get('hia', 0.0),
                'pgp_substrate score': scores['metric_scores'].get('Absorption', {}).get('pgp_substrate', 0.0),
                'pgp_substrate value': metrics.absorption.get('pgp_substrate', 0.0),
                'pgp_inhibitor score': scores['metric_scores'].get('Absorption', {}).get('pgp_inhibitor', 0.0),
                'pgp_inhibitor value': metrics.absorption.get('pgp_inhibitor', 0.0),
                
                # Distribution
                'vdss score': scores['metric_scores'].get('Distribution', {}).get('vdss', 0.0),
                'vdss value': metrics.distribution.get('vdss', 0.0),
                'ppb score': scores['metric_scores'].get('Distribution', {}).get('ppb', 0.0),
                'ppb value': metrics.distribution.get('ppb', 0.0),
                'bbb score': scores['metric_scores'].get('Distribution', {}).get('bbb', 0.0),
                'bbb value': metrics.distribution.get('bbb', 0.0),
                'fu score': scores['metric_scores'].get('Distribution', {}).get('fu', 0.0),
                'fu value': metrics.distribution.get('fu', 0.0),
                
                # Metabolism
                'cyp2c9_inhibition score': scores['metric_scores'].get('Metabolism', {}).get('cyp2c9_inhibition', 0.0),
                'cyp2c9_inhibition value': metrics.metabolism.get('cyp2c9_inhibition', 0.0),
                'cyp2d6_inhibition score': scores['metric_scores'].get('Metabolism', {}).get('cyp2d6_inhibition', 0.0),
                'cyp2d6_inhibition value': metrics.metabolism.get('cyp2d6_inhibition', 0.0),
                'cyp3a4_inhibition score': scores['metric_scores'].get('Metabolism', {}).get('cyp3a4_inhibition', 0.0),
                'cyp3a4_inhibition value': metrics.metabolism.get('cyp3a4_inhibition', 0.0),
                
                # Excretion
                'half_life score': scores['metric_scores'].get('Excretion', {}).get('half_life', 0.0),
                'half_life value': metrics.excretion.get('half_life', 0.0),
                'clearance score': scores['metric_scores'].get('Excretion', {}).get('clearance', 0.0),
                'clearance value': metrics.excretion.get('clearance', 0.0),
                
                # Toxicity
                'herg score': scores['metric_scores'].get('Toxicity', {}).get('herg', 0.0),
                'herg value': metrics.toxicity.get('herg', 0.0)
            }
            
            all_molecules.append(row)
        
        # Create DataFrame with results
        self.df = pd.DataFrame(all_molecules)
        self.df.sort_values('Overall Score', ascending=False, inplace=True)
        self.df.reset_index(drop=True, inplace=True)
        
        # Generate timestamp for file naming
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if self.reports_dir:
            output_dir = Path(self.reports_dir)
            logging.info("💾 Saving ranking results...")
            
            # Generate reports with correct ranks after sorting
            for idx, row in enumerate(self.df.itertuples(), 1):
                metrics = next(m for m in metrics_list if m.smiles == row.SMILES)
                scores = molecule_scores[row.SMILES]  # Get stored scores
                logging.info(f"📝 Generating report for {metrics.molecular_formula} (Rank {idx})")
                self.generate_molecule_report(metrics, scores, self.reports_dir, idx)
            
            # Save rankings
            self.save_rankings(output_dir, timestamp)
            
            # Generate final report
            self.generate_report(str(output_dir), timestamp)
            
            logging.info(f"✅ Ranking complete! Results saved to {output_dir}")
        
        return self.df

    def save_rankings(self, output_dir: Path, timestamp: str):
        """Save rankings to CSV file"""
        if self.df is not None:
            # Create rankings directory
            rankings_dir = output_dir / 'rankings'
            rankings_dir.mkdir(parents=True, exist_ok=True)
            
            # Save rankings with timestamp
            rankings_file = rankings_dir / f'rankings_{timestamp}.csv'
            self.df.to_csv(rankings_file, index=False)
            return rankings_file
        return None

    def generate_reports(self, output_dir: Path, timestamp: str):
        """Generate reports for each molecule"""
        if self.df is None:
            return []
        
        # Create molecule reports directory
        reports_dir = output_dir / 'molecule_reports'
        reports_dir.mkdir(parents=True, exist_ok=True)
        
        report_files = []
        for idx, (_, row) in enumerate(self.df.iterrows(), 1):
            # Generate report using molecular formula in filename
            mol_formula = row['Molecular Formula']
            report_base = f"{mol_formula}_{timestamp}"
            pdf_file = reports_dir / f"{report_base}.pdf"
            csv_file = reports_dir / f"{report_base}.csv"
            
            # Generate report content
            molecule_data = {
                "smiles": row['SMILES'],
                "molecular_formula": row['Molecular Formula'],
                "molecular_weight": row['Molecular Weight'],
                "rank": row.name + 1,  # Add rank to molecule data
                "physicochemical": {
                    'tpsa': row['tpsa value'],
                    'logp': row['logp value'],
                    'qed': row['qed value']
                },
                "medicinal_chemistry": {
                    'synthetic_accessibility': row['synthetic_accessibility value']
                },
                "lipophilicity": {
                    'ilogp': row['ilogp value'],
                    'xlogp3': row['xlogp3 value'],
                    'wlogp': row['wlogp value'],
                    'mlogp': row['mlogp value'],
                    'silicos_it': row['silicos_it value'],
                    'consensus_logp': row['consensus_logp value']
                },
                "druglikeness": {
                    'lipinski_violations': row['lipinski_violations value'],
                    'bioavailability_score': row['bioavailability_score value']
                },
                "absorption": {
                    'caco2': row['caco2 value'],
                    'pampa': row['pampa value'],
                    'mdck': row['mdck value'],
                    'hia': row['hia value'],
                    'pgp_substrate': row['pgp_substrate value'],
                    'pgp_inhibitor': row['pgp_inhibitor value']
                },
                "distribution": {
                    'vdss': row['vdss value'],
                    'ppb': row['ppb value'],
                    'bbb': row['bbb value'],
                    'fu': row['fu value']
                },
                "metabolism": {
                    'cyp2c9_inhibition': row['cyp2c9_inhibition value'],
                    'cyp2d6_inhibition': row['cyp2d6_inhibition value'],
                    'cyp3a4_inhibition': row['cyp3a4_inhibition value']
                },
                "excretion": {
                    'half_life': row['half_life value'],
                    'clearance': row['clearance value']
                },
                "toxicity": {
                    'herg': row['herg value']
                },
                "scores": {
                    'overall_score': row['Overall Score'],
                    'physicochemical_score': row['Physicochemical Score'],
                    'medicinal_chemistry_score': row['Medicinal Chemistry Score'],
                    'druglikeness_score': row['Druglikeness Score'],
                    'absorption_score': row['Absorption Score'],
                    'distribution_score': row['Distribution Score'],
                    'metabolism_score': row['Metabolism Score'],
                    'excretion_score': row['Excretion Score'],
                    'toxicity_score': row['Toxicity Score']
                }
            }

            # Generate reports
            # TODO: ADDED FALSE HERE TO IGNORE PDF GENERTATION FOR THE SAKE OF TOXICITY PAPER - User wants to control this
            generate_enhanced_report(molecule_data, str(reports_dir), generate_pdf=self.generate_pdf, generate_csv=self.generate_csv)
            report_files.extend([pdf_file, csv_file])
        
        return report_files

    def generate_report(self, output_dir: str, timestamp: str):
        """Generate final ranking report"""
        if self.df is None or self.df.empty:
            return

        # Create rankings directory
        rankings_dir = Path(output_dir) / "rankings"
        rankings_dir.mkdir(parents=True, exist_ok=True)

        # Save rankings to CSV
        csv_file = rankings_dir / f"rankings_{timestamp}.csv"
        self.df.to_csv(csv_file, index=False)
        
        # Generate PDF report using report generator
        from .reports import generate_ranking_report
        pdf_file = rankings_dir / f"ranking_report_{timestamp}.pdf"
        # TODO: MODIFIED THE TOP 10 T0 100
        generate_ranking_report(str(csv_file), str(pdf_file), top_n=100)
        
        print(f"\nRanking results saved to:")
        print(f"CSV: {csv_file}")
        print(f"PDF: {pdf_file}")
        
    def get_ranked_data_as_dict(self) -> List[Dict[str, Any]]:
        """Returns the ranked molecules data as a list of dictionaries.

        Each dictionary in the list represents a molecule and its associated data
        as present in the internal ranked DataFrame.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries, where each dictionary
                                 is a row from the ranked DataFrame. Returns an
                                 empty list if ranking has not been performed or
                                 resulted in an empty DataFrame.
        
        Example:
            >>> # Assuming 'ranker' is an instance of SmallMoleculeRankerV4
            >>> # and rank_molecules() has been called.
            >>> ranked_list_of_dicts = ranker.get_ranked_data_as_dict()
            >>> if ranked_list_of_dicts:
            ...     print(f"Top ranked molecule SMILES: {ranked_list_of_dicts[0]['SMILES']}")
        """
        if self.df is not None and not self.df.empty:
            return self.df.to_dict(orient="records")
        return []

def rank_molecules_from_metrics(
    metrics_input: Union[str, List[dict]],
    output_dir: str,
    config: Optional[ScoringConfig] = None
) -> pd.DataFrame:
    """
    Main function to rank molecules from metrics file or list.
    
    Args:
        metrics_input: Either a path to metrics file or a list of metrics dictionaries
        output_dir: Directory to save ranking reports
        config: Optional custom scoring configuration
        
    Returns:
        DataFrame with ranked molecules
    """
    # Initialize validator and ranker
    validator = SmallMoleculeValidator()
    ranker = SmallMoleculeRankerV4(config)
    ranker.set_output_directory(output_dir)
    
    # Process molecules based on input type
    if isinstance(metrics_input, str):
        # Process from file
        results = validator.process_molecules(metrics_input, output_dir)
        metrics_list = validator.get_successful_metrics(results)
    elif isinstance(metrics_input, list):
        # Directly use the provided metrics list
        metrics_list = metrics_input
    else:
        raise TypeError("metrics_input must be either a file path (str) or a list of metrics dictionaries")
    
    if not metrics_list:
        raise ValueError("No valid molecules found in input")
    
    # Rank molecules and generate reports
    rankings = ranker.rank_molecules(metrics_list)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    ranker.save_rankings(output_dir, timestamp)
    ranker.generate_reports(output_dir, timestamp)
    
    return rankings

if __name__ == "__main__":
    # Example usage
    metrics_file = "example_molecules.txt"
    output_dir = "ranking_results"
    rankings = rank_molecules_from_metrics(metrics_file, output_dir)
    print(f"Ranked {len(rankings)} molecules successfully")
