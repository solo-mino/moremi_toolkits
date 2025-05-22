"""
Enhanced Molecule Analysis Report Generator

This module generates comprehensive PDF and CSV reports for molecular analysis results,
including detailed property breakdowns, visualizations, and structure rendering.
"""

import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import importlib.resources as pkg_resources

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Define color scheme based on the image
COLOR_SCHEME = {
    'header_bg': '#4B89FF',  # Light blue for headers
    'text': '#000000',       # Black for text
    'grid': '#E5E5E5',       # Light gray for grid lines
    'chart': {
        'line': '#4B89FF',   # Light blue for radar chart line
        'fill': '#4B89FF33', # Semi-transparent blue for radar fill
        'grid': '#E5E5E5',   # Light gray for radar grid
        'text': '#333333'    # Dark gray for labels
    }
}

class EnhancedMoleculeReport(FPDF):
    """Enhanced PDF report generator for molecular analysis."""
    
    def __init__(self, output_dir: str):
        super().__init__('P', 'mm', 'A4')
        self.output_dir = Path(output_dir)
        self.set_auto_page_break(auto=True, margin=15)
        
        # Set up fonts using importlib.resources
        try:
            package_ref = pkg_resources.files('moremi_biokit')
            assets_dir_ref = package_ref.joinpath('assets')
            
            regular_font_ref = assets_dir_ref.joinpath('fonts', 'SpaceGrotesk-Regular.ttf')
            bold_font_ref = assets_dir_ref.joinpath('fonts', 'SpaceGrotesk-Bold.ttf')

            with pkg_resources.as_file(regular_font_ref) as regular_font_path:
                self.add_font('SpaceGrotesk', '', str(regular_font_path), uni=True)
            
            with pkg_resources.as_file(bold_font_ref) as bold_font_path:
                self.add_font('SpaceGrotesk', 'B', str(bold_font_path), uni=True)
                
        except FileNotFoundError:
            # Fallback or error handling if fonts are not found
            # For now, we'll let FPDF handle it or raise an error if fonts are critical
            print("Warning: Custom fonts not found. Default fonts will be used.")
        
        # Initialize section counter
        self._section_counter = 0
        
        # Set margins
        self.set_margins(15, 15, 15)

    def header(self):
        """Add header to pages"""
        # Add logo using importlib.resources
        try:
            package_ref = pkg_resources.files('moremi_biokit')
            assets_dir_ref = package_ref.joinpath('assets')
            logo_ref = assets_dir_ref.joinpath('minologo.png')
            with pkg_resources.as_file(logo_ref) as logo_path:
                if logo_path.exists(): # Check if path from context manager is valid
                    self.image(str(logo_path), 10, 15, 15)
        except FileNotFoundError:
            print("Warning: Logo image 'minologo.png' not found in package assets.")
        
        # Add title
        self.set_font('SpaceGrotesk', 'B', 15)
        title = "Comprehensive Small Molecule Analysis Report"
        self.cell(0, 10, title, 0, 1, 'C')
        
        # Add subtitle with molecular formula and date
        if hasattr(self, 'molecular_formula'):
            self.set_font('SpaceGrotesk', '', 11)
            self.cell(0, 5, f"Molecule: {self.molecular_formula}", 0, 1, 'C')
            self.cell(0, 5, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 0, 1, 'C')
        
        # Add line
        self.line(10, self.get_y() + 2, self.w - 10, self.get_y() + 2)
        self.ln(10)

    def footer(self):
        """Create footer with page number."""
        self.set_y(-15)
        self.set_font('SpaceGrotesk', '', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

    def add_molecule_structure(self, smiles: str):
        """Add molecular structure image generated from SMILES."""
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return
            
            # Generate 2D coordinates
            AllChem.Compute2DCoords(mol)
            
            # Create image
            img = Draw.MolToImage(mol, size=(300, 300))
            
            # Save temporary image
            temp_path = self.output_dir / 'temp_structure.png'
            img.save(str(temp_path))
            
            # Add to PDF
            self.image(str(temp_path), x=50, y=None, w=100)
            
            # Clean up
            os.remove(temp_path)
            
        except Exception as e:
            print(f"Error generating molecular structure: {e}")

    def add_section_title(self, title: str):
        """Add formatted section title."""
        self._section_counter += 1  # Increment counter
        self.set_font('SpaceGrotesk', 'B', 14)  # Increased font size
        self.set_fill_color(75, 137, 255)  # Bright blue background
        self.set_text_color(255, 255, 255)  # White text
        self.cell(0, 10, f"{self._section_counter}. {title}", 0, 1, 'L', True)
        self.set_text_color(0, 0, 0)  # Reset text color
        self.ln(4)

    def add_subsection_title(self, title: str):
        """Add formatted subsection title."""
        self.set_font('SpaceGrotesk', 'B', 12)  # Medium font size
        self.set_text_color(75, 137, 255)  # Bright blue text
        self.cell(0, 8, title, 'B', 1, 'L')  # Added bottom border
        self.set_text_color(0, 0, 0)  # Reset text color
        self.ln(3)

    def add_method_title(self, title: str):
        """Add formatted method title."""
        self.set_font('SpaceGrotesk', '', 10)  # Regular font, smaller size
        self.set_text_color(100, 100, 100)  # Gray color
        self.cell(0, 6, "" + title, 0, 1, 'L')  # Changed to arrow symbol
        self.set_text_color(0, 0, 0)  # Reset text color
        self.ln(2)

    def create_property_table(self, data: Dict[str, Any], headers: List[str]):
        """Create a formatted table for property data."""
        # Set up table formatting
        self.set_font('SpaceGrotesk', '', 9)
        col_width = self.get_string_width('X' * 30)
        line_height = 7
        
        # Add headers
        self.set_font('SpaceGrotesk', 'B', 9)
        self.set_fill_color(236, 240, 241)  # Light gray background
        for header in headers:
            self.cell(col_width, line_height, header, 1, 0, 'L', True)
        self.ln()
        
        # Add data
        self.set_font('SpaceGrotesk', '', 9)
        for key, value in data.items():
            if isinstance(value, dict):
                # Handle nested dictionaries
                for sub_key, sub_value in value.items():
                    self.cell(col_width, line_height, f"{key} ({sub_key})", 1)
                    self.cell(col_width, line_height, str(sub_value), 1)
                    self.ln()
            else:
                self.cell(col_width, line_height, str(key), 1)
                self.cell(col_width, line_height, str(value), 1)
                self.ln()

    def format_value(self, value, unit=None, decimals=3, scientific=False):
        """Format a value with its unit"""
        if value is None or value == "N/A":
            return "N/A"
            
        if isinstance(value, bool):
            return str(value)
            
        if isinstance(value, (int, float)):
            if scientific:
                formatted = f"{value:.{decimals}e}"
            else:
                if float(value).is_integer():
                    formatted = f"{int(value)}"
                else:
                    formatted = f"{value:.{decimals}f}"
        else:
            formatted = str(value)
            
        if unit:
            return f"{formatted} {unit}"
        return formatted

    def add_basic_properties(self, data):
        """Add basic molecular properties section"""
        self.add_section_title("Physicochemical Properties")
        
        basic_props = {
            # "SMILES": (data.get("smiles", "N/A"), None),
            # "Molecular Formula": (data.get("molecular_formula", "N/A"), None),
            "Molecular Weight": (data.get("molecular_weight", "N/A"), "g/mol"),
            "Molar Refractivity": (data['physicochemical'].get("molar_refractivity", "N/A"), None),
            "Molecular Density": (data['physicochemical'].get("molecular_density", "N/A"), "g/cm³"),
            "Molecular Volume": (data['physicochemical'].get("molecular_volume", "N/A"), "Å³, Van der Waals volume"),
            "Number of Rotatable Bonds": (data['physicochemical'].get("num_rotatable_bonds", "N/A"), None),
            "H-Bond Donors (HBD)": (data['physicochemical'].get("num_hbd", "N/A"), None),
            "H-Bond Acceptors (HBA)": (data['physicochemical'].get("num_hba", "N/A"), None),
            'Ring Count': (data['physicochemical'].get("ring_count", "N/A"), None),
            "Max Ring Size": (data['physicochemical'].get("max_ring_size", "N/A"), None),
            "Number of Heavy Atoms": (data['physicochemical'].get("num_heavy_atoms", "N/A"), None),
            "Number of Aromatic Heavy Atoms": (data['physicochemical'].get("num_aromatic_heavy_atoms", "N/A"), None),
            "Fraction CSP3": (data['physicochemical'].get("fraction_csp3", "N/A"), None),
            "LogP": (data['physicochemical'].get("logp", "log-ratio"), None),
            "TPSA": (data['physicochemical'].get("tpsa", "N/A"), "Å²"),
            "Molecular Flexibility": (data['physicochemical'].get("molecular_flexibility", "N/A"), None),
            "Molecular Rigidity": (data['physicochemical'].get("molecular_rigidity", "N/A"), None),
            "Number of non-C/H atoms": (data['physicochemical'].get("heteroatom_count", "N/A"), None),
            "Formal Charge": (data['physicochemical'].get("formal_charge", "N/A"), "e, Net Molecular Charge"),
            "Stereogenic Centers": (data['physicochemical'].get("stereogenic_centers", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in basic_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_lipophilicity_section(self, data):
        """Add lipophilicity section"""
        if 'lipophilicity' not in data:
            return
            
        self.add_section_title("Lipophilicity")
        
        lipophilicity_props = {
            "iLogP": (data['lipophilicity'].get("ilogp", "N/A"), None),
            "XLogP3": (data['lipophilicity'].get("xlogp3", "N/A"), None),
            "WLogP": (data['lipophilicity'].get("wlogp", "N/A"), None),
            "MLogP": (data['lipophilicity'].get("mlogp", "N/A"), None),
            "SILICOS-IT LogP": (data['lipophilicity'].get("silicos_it", "N/A"), None),
            "Consensus LogP": (data['lipophilicity'].get("consensus_logp", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in lipophilicity_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_solubility_section(self, data):
        """Add solubility section"""
        if 'solubility' not in data:
            return
            
        self.add_section_title("Solubility")
        
        # Add each solubility method
        for method in ['ESOL', 'Ali', 'SILICOS-IT']:
            if method in data['solubility']:
                method_data = data['solubility'][method]
                solubility_props = {
                    "LogS": (method_data.get("log_s", "N/A"), None),
                    "Solubility": (method_data.get("solubility_mg_ml", "N/A"), "mg/mL"),
                    "Molar Solubility": (method_data.get("solubility_mol_l", "N/A"), "mol/L"),
                    "Class": (method_data.get("class", "N/A"), None)
                }
                
                self.add_subsection_title(f"{method} Method")
                formatted_props = {k: self.format_value(v[0], v[1], scientific=True if k != "Class" else False) 
                                 for k, v in solubility_props.items()}
                self.create_property_table(formatted_props, ["Property", "Value"])
                self.ln(3)

    def add_druglikeness_section(self, data):
        """Add druglikeness section"""
        if 'druglikeness' not in data:
            return
            
        self.add_section_title("Drug-likeness Assessment")
        
        # Bioavailability score
        bioavailability_props = {
            "Bioavailability Score": (data['druglikeness'].get("bioavailability_score", "N/A"), None),
        }
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in bioavailability_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(3)
        
        # Rule-based evaluations
        rules = ['lipinski', 'ghose', 'veber', 'egan', 'muegge']
        for rule in rules:
            if rule in data['druglikeness']:
                rule_data = data['druglikeness'][rule]
                self.add_subsection_title(f"{rule.title()} Rules")
                
                rule_props = {
                    "Status": ("Passes" if rule_data.get("passes", False) else "Fails", None)
                }
                formatted_props = {k: self.format_value(v[0], v[1]) for k, v in rule_props.items()}
                self.create_property_table(formatted_props, ["Property", "Value"])
                
                violations = rule_data.get("violations", [])
                if violations:
                    self.set_font('SpaceGrotesk', '', 10)
                    self.cell(0, 6, "Violations:", 0, 1, 'L')
                    for violation in violations:
                        self.cell(0, 6, f"• {violation}", 0, 1, 'L')
                self.ln(3)

    def add_medicinal_chemistry_section(self, data):
        """Add medicinal chemistry section"""
        if 'medicinal_chemistry' not in data:
            return
            
        self.add_section_title("Medicinal Chemistry")
        
        # Synthetic Properties
        self.add_subsection_title("Synthetic Properties")
        med_chem_props = {
            "Synthetic Accessibility": (data['medicinal_chemistry'].get("synthetic_accessibility", "N/A"), None),
            "QED Score": (data['medicinal_chemistry'].get("qed", "N/A"), None)
        }
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in med_chem_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(3)
        
        # Structural Alerts
        if "pains" in data['medicinal_chemistry']:
            self.add_subsection_title("PAINS Analysis")
            pains_data = data['medicinal_chemistry']["pains"]
            pains_props = {
                "Alert Count": (pains_data.get("alert_count", 0), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in pains_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            
            alerts = pains_data.get("alerts", [])
            if alerts:
                self.set_font('SpaceGrotesk', '', 10)
                self.cell(0, 6, "Alert Details:", 0, 1, 'L')
                for alert in alerts:
                    self.cell(0, 6, f"• {alert}", 0, 1, 'L')
            self.ln(3)
        
        if "brenk" in data['medicinal_chemistry']:
            self.add_subsection_title("Brenk Analysis")
            brenk_data = data['medicinal_chemistry']["brenk"]
            brenk_props = {
                "Alert Count": (brenk_data.get("alert_count", 0), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in brenk_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            
            alerts = brenk_data.get("alerts", [])
            if alerts:
                self.set_font('SpaceGrotesk', '', 10)
                self.cell(0, 6, "Alert Details:", 0, 1, 'L')
                for alert in alerts:
                    self.cell(0, 6, f"• {alert}", 0, 1, 'L')
            self.ln(3)
        
        # Leadlikeness
        if "leadlikeness" in data['medicinal_chemistry']:
            self.add_subsection_title("Leadlikeness")
            leadlikeness_data = data['medicinal_chemistry']["leadlikeness"]
            leadlikeness_props = {
                "Pass/Fail Status": ("Pass" if leadlikeness_data.get("pass", False) else "Fail", None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in leadlikeness_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            
            violations = leadlikeness_data.get("violations", [])
            if violations:
                self.set_font('SpaceGrotesk', '', 10)
                self.cell(0, 6, "Violations:", 0, 1, 'L')
                for violation in violations:
                    self.cell(0, 6, f"• {violation}", 0, 1, 'L')
            self.ln(3)

    def add_admet_section(self, data):
        """Add ADMET properties section"""
        # Absorption
        if 'absorption' in data:
            self.add_section_title("Absorption")
            absorption_props = {
                "Caco2 Permeability": (data['absorption'].get("caco2", "N/A"), "log(10-⁶ cm/s)"),
                "PAMPA Permeability": (data['absorption'].get("pampa", "N/A"), None),
                "MDCK Permeability": (data['absorption'].get("mdck", "N/A"), "log(10-⁶ cm/s)"),
                "P-gp Substrate": (data['absorption'].get("pgp_substrate", "N/A"), None),
                "P-gp Inhibitor": (data['absorption'].get("pgp_inhibitor", "N/A"), None),
                "HIA": (data['absorption'].get("hia", "N/A"), None),
                "Oral Bioavailability": (data['absorption'].get("oral_bioavailability", "N/A"), None),
                "GI Absorption": (data['absorption'].get("gi_absorption", "N/A"), None),
                "Log Kp": (data['absorption'].get("log_kp", "N/A"), "cm/s")
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in absorption_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(5)

        # Distribution
        if 'distribution' in data:
            self.add_section_title("Distribution")
            distribution_props = {
                "VDss": (data['distribution'].get("vdss", "N/A"), "L/kg"),
                "Plasma Protein Binding Rate": (data['distribution'].get("ppb", "N/A"), "%"),
                "BBB Penetration": (data['distribution'].get("bbb", "N/A"), None),
                "Fraction Unbound": (data['distribution'].get("fu", "N/A"), "%"),
                'Hydration Free Energy': (data['distribution'].get("hydration_free_energy", "N/A"), "kcal/mol"),
                "Aqueous Solubility": (data['distribution'].get("aqueous_solubility", "N/A"), "log(mol/L)"),
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in distribution_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(5)

        # Metabolism
        if 'metabolism' in data:
            self.add_section_title("Metabolism")
            
            # CYP Inhibition
            self.add_subsection_title("CYP Inhibition")
            cyp_inhibition_props = {
                "CYP1A2": (data['metabolism'].get("cyp1a2_inhibition", "N/A"), None),
                "CYP2C19": (data['metabolism'].get("cyp2c19_inhibition", "N/A"), None),
                "CYP2C9": (data['metabolism'].get("cyp2c9_inhibition", "N/A"), None),
                "CYP2D6": (data['metabolism'].get("cyp2d6_inhibition", "N/A"), None),
                "CYP3A4": (data['metabolism'].get("cyp3a4_inhibition", "N/A"), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in cyp_inhibition_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(3)
            
            # CYP Substrate
            self.add_subsection_title("CYP Substrate")
            cyp_substrate_props = {
                "CYP2C9": (data['metabolism'].get("cyp2c9_substrate", "N/A"), None),
                "CYP2D6": (data['metabolism'].get("cyp2d6_substrate", "N/A"), None),
                "CYP3A4": (data['metabolism'].get("cyp3a4_substrate", "N/A"), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in cyp_substrate_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(3)
            
            # Clearance
            self.add_subsection_title("Clearance")
            clearance_props = {
                "Hepatic Clearance": (data['metabolism'].get("hepatic_clearance", "N/A"), "µL/min/106 cells"),
                "Microsomal Clearance": (data['metabolism'].get("microsomal_clearance", "N/A"), "µL/min/mg")
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in clearance_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(5)

        # Excretion
        if 'excretion' in data:
            self.add_section_title("Excretion")
            excretion_props = {
                "Half-life": (data['excretion'].get("half_life", "N/A"), "hr"),
                "Clearance": (data['excretion'].get("clearance", "N/A"), "µL/min/10⁶")
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in excretion_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(5)

        # Toxicity
        if 'toxicity' in data:
            self.add_section_title("Toxicity")
            
            # Basic toxicity properties
            basic_tox_props = {
                "hERG Inhibition": (data['toxicity'].get("herg", "N/A"), None),
                "Drug-Induced Liver Injury": (data['toxicity'].get("drug_induced_liver_injury", "N/A"), None),
                "Ames Mutagenicity": (data['toxicity'].get("ames_mutagenicity", "N/A"), None),
                "Skin Sensitization": (data['toxicity'].get("skin_sensitization", "N/A"), None),
                "Lethal Dose 50": (data['toxicity'].get("lethal_dose_50", "N/A"), "log(1/(mol/kg))"),
                "Clinical Toxicity": (data['toxicity'].get("clinical_toxicity", "N/A"), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in basic_tox_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(3)
            
            # Nuclear Receptor Pathways
            self.add_subsection_title("Nuclear Receptor (NR) Pathways")
            nr_props = {
                "NR-AR: Androgen Receptor Binding": (data['toxicity'].get("androgen_receptor_binding", "N/A"), None),
                "NR-AR-LBD: Androgen Receptor LBD": (data['toxicity'].get("androgen_receptor_lbd", "N/A"), None),
                "NR-AhR: Aryl Hydrocarbon Receptor": (data['toxicity'].get("aryl_hydrocarbon_receptor", "N/A"), None),
                "NR-Aromatase: Aromatase Inhibition": (data['toxicity'].get("aromatase_inhibition", "N/A"), None),
                "NR-ER: Estrogen Receptor Binding": (data['toxicity'].get("estrogen_receptor_binding", "N/A"), None),
                "NR-ER-LBD: Estrogen Receptor LBD": (data['toxicity'].get("estrogen_receptor_lbd", "N/A"), None),
                "NR-PPAR-gamma: PPAR Gamma": (data['toxicity'].get("ppar_gamma", "N/A"), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in nr_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(3)
            
            # Stress Response Pathways
            self.add_subsection_title("Stress Response (SR) Pathways")
            sr_props = {
                "SR-ARE: Oxidative Stress Response": (data['toxicity'].get("oxidative_stress_response", "N/A"), None),
                "SR-ATAD5: DNA Damage Response": (data['toxicity'].get("dna_damage_response", "N/A"), None),
                "SR-HSE: Heat Shock Response": (data['toxicity'].get("heat_shock_response", "N/A"), None),
                "SR-MMP: Mitochondrial Toxicity": (data['toxicity'].get("mitochondrial_toxicity", "N/A"), None),
                "SR-p53: p53 Pathway Activation": (data['toxicity'].get("p53_pathway_activation", "N/A"), None)
            }
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in sr_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(5)

    def add_radar_chart(self, scores):
        """Add radar chart from category scores"""
        if not scores:
            return
        
        # Remove lipophilicity from scores
        scores.pop('Lipophilicity', None)
        
        # Create radar chart
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
        
        # Get category names and values
        categories = list(scores.keys())
        values = list(scores.values())
        
        # Number of variables
        num_vars = len(categories)
        
        # Compute angle for each axis
        angles = [n / float(num_vars) * 2 * np.pi for n in range(num_vars)]
        angles += angles[:1]
        
        # Plot data
        values += values[:1]
        ax.plot(angles, values, color='#3498db', linewidth=1.5)
        ax.fill(angles, values, color='#3498db', alpha=0.25)
        
        # Set chart properties
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, size=10)
        ax.set_ylim(0, 1)
        
        # Add grid
        ax.grid(True, color='#95a5a6', alpha=0.2)
        
        # Add subtle background colors
        for i in range(5):
            ax.fill(angles, [i/4]*len(angles), color='gray', alpha=0.05)
        
        # Add value labels
        for angle, value in zip(angles[:-1], values[:-1]):
            if value > 0:  # Only show non-zero values
                ax.text(angle, value + 0.05, f'{value:.3f}', 
                       ha='center', va='center')
        
        # Save the chart
        self.temp_chart_path = Path(self.output_dir) / 'temp_radar_chart.png'
        plt.savefig(str(self.temp_chart_path), dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        # Add the chart to the PDF
        if self.temp_chart_path.exists():
            self.image(str(self.temp_chart_path), x=30, y=None, w=150)
            try:
                os.remove(self.temp_chart_path)
            except:
                pass

    def create_report(self, data: Dict):
        """Create the complete report"""
        # Store molecular formula for header
        self.molecular_formula = data.get("molecular_formula", "Unknown")
        
        # Add new page with header
        self.add_page()
        
        # Add molecule structure
        self.add_molecule_structure(data.get("smiles", ""))
        
        # Add SMILES and Molecular Formula
        self.set_font('SpaceGrotesk', '', 12)
        self.cell(0, 10, f"SMILES: {data.get('smiles', 'N/A')}", 0, 1, 'L')
        self.cell(0, 10, f"MOLECULAR FORMULA: {data.get('molecular_formula', 'N/A')}", 0, 1, 'L')
        self.set_font('SpaceGrotesk', '', 10)  # Reset font to default
        
        # Add all property sections (each on a new page with header)
        
        self.add_basic_properties(data)
        
        self.add_lipophilicity_section(data)
        
        self.add_solubility_section(data)
        
        self.add_medicinal_chemistry_section(data)

        self.add_druglikeness_section(data)
        
        self.add_admet_section(data)
        
        # Add radar chart from ranking scores if available
        if 'scores' in data and 'category_scores' in data['scores']:
            self.add_page()
            self.add_section_title("Molecular Properties Overview")
            self.add_radar_chart(data['scores']['category_scores'])
            self.ln(5)

    def generate_pdf_report(self, data: Dict, output_path: str):
        """Generate the complete PDF report."""
        self.create_report(data)
        
        # Save the report
        self.output(output_path)

def generate_csv_report(data: Dict[str, Any], output_path: str):
    """Generate CSV report with hierarchical structure."""
    # Flatten the nested dictionary while preserving hierarchy
    flat_data = {}
    
    def flatten_dict(d: Dict, prefix: str = ""):
        for key, value in d.items():
            if isinstance(value, dict):
                flatten_dict(value, f"{prefix}{key}_")
            else:
                flat_data[f"{prefix}{key}"] = value
    
    flatten_dict(data)
    
    # Convert to DataFrame and save
    df = pd.DataFrame([flat_data])
    df.to_csv(output_path, index=False)

def generate_enhanced_report(molecule_data: Dict[str, Any], output_dir: str, generate_pdf: Optional[bool] = True, generate_csv: Optional[bool] = True):
    """Generate enhanced PDF and CSV reports for molecular analysis."""
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize report
    report = EnhancedMoleculeReport(output_dir)
    
    # Generate timestamp for file names
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Get rank from molecule data, default to 0 if not provided
    rank = molecule_data.get('rank', 0)
    
    pdf_path = None
    csv_path = None

    # Create PDF report with rank in filename
    if generate_pdf:
        pdf_path = output_dir / f'r{rank}_{molecule_data["molecular_formula"]}_{timestamp}.pdf'
        report.generate_pdf_report(molecule_data, str(pdf_path))
    
    # Create CSV report with rank in filename
    if generate_csv:
        csv_path = output_dir / f'r{rank}_{molecule_data["molecular_formula"]}_{timestamp}.csv'
        generate_csv_report(molecule_data, str(csv_path))
    
    if generate_pdf and generate_csv:
        return pdf_path, csv_path
    elif generate_pdf:
        return pdf_path, None
    elif generate_csv:
        return None, csv_path
    return None, None

def generate_radar_chart(data: Dict[str, float], output_dir: str):
    """Create a radar chart for a category of properties."""
    # Create radar chart
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    
    # Prepare data
    categories = list(data.keys())
    values = list(data.values())
    
    # Number of variables
    num_vars = len(categories)
    
    # Compute angle for each axis
    angles = [n / float(num_vars) * 2 * np.pi for n in range(num_vars)]
    angles += angles[:1]
    
    # Plot data
    values += values[:1]
    ax.plot(angles, values, color='blue', linewidth=2)
    ax.fill(angles, values, color='#010175', alpha=0.05)
    
    # Set chart properties
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_title("Overall Category Scores Analysis", pad=20)
    
    # Save chart
    temp_path = output_dir / f'radar_chart.png'
    plt.savefig(str(temp_path), dpi=300, bbox_inches='tight')
    plt.close()

def generate_enhanced_report_v2(molecule_data: Dict[str, Any], output_dir: str):
    """Generate enhanced PDF and CSV reports for molecular analysis."""
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create PDF report
    pdf = FPDF()
    
    # Set up fonts
    pdf.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf', uni=True)
    pdf.add_font('DejaVu', 'B', 'DejaVuSansCondensed-Bold.ttf', uni=True)
    
    # Title Page
    pdf.add_page()
    pdf.set_font('DejaVu', 'B', 16)
    pdf.cell(0, 10, 'Molecular Analysis Report', 0, 1, 'C')
    pdf.ln(10)
    
    # Basic Information
    pdf.set_font('DejaVu', 'B', 14)
    pdf.cell(0, 10, 'Basic Information', 0, 1, 'L')
    pdf.set_font('DejaVu', '', 12)
    pdf.cell(0, 10, f'SMILES: {molecule_data["smiles"]}', 0, 1, 'L')
    pdf.cell(0, 10, f'Molecular Formula: {molecule_data["molecular_formula"]}', 0, 1, 'L')
    pdf.cell(0, 10, f'Molecular Weight: {molecule_data["molecular_weight"]:.2f}', 0, 1, 'L')
    pdf.ln(10)
    
    # Property Groups
    property_groups = [
        ("Physicochemical Properties", "physicochemical"),
        ("Medicinal Chemistry", "medicinal_chemistry"),
        ("Lipophilicity", "lipophilicity"),
        ("Drug-likeness", "druglikeness"),
        ("Absorption", "absorption"),
        ("Distribution", "distribution"),
        ("Metabolism", "metabolism"),
        ("Excretion", "excretion"),
        ("Toxicity", "toxicity")
    ]
    
    for title, key in property_groups:
        if key in molecule_data and molecule_data[key]:
            pdf.add_page()
            pdf.set_font('DejaVu', 'B', 14)
            pdf.cell(0, 10, title, 0, 1, 'L')
            pdf.ln(5)
            
            pdf.set_font('DejaVu', '', 12)
            for metric, value in molecule_data[key].items():
                if isinstance(value, (int, float)):
                    pdf.cell(0, 10, f'{metric}: {value:.4f}', 0, 1, 'L')
                else:
                    pdf.cell(0, 10, f'{metric}: {value}', 0, 1, 'L')
            pdf.ln(5)
    
    # Overall Analysis
    if "scores" in molecule_data:
        pdf.add_page()
        pdf.set_font('DejaVu', 'B', 14)
        pdf.cell(0, 10, 'Overall Analysis', 0, 1, 'L')
        pdf.ln(5)
        
        scores = molecule_data["scores"]
        pdf.set_font('DejaVu', '', 12)
        pdf.cell(0, 10, f'Overall Score: {scores["overall_score"]:.4f}', 0, 1, 'L')
        
        pdf.ln(5)
        pdf.set_font('DejaVu', 'B', 12)
        pdf.cell(0, 10, 'Category Scores:', 0, 1, 'L')
        pdf.set_font('DejaVu', '', 12)
        for category, score in scores["category_scores"].items():
            pdf.cell(0, 10, f'{category}: {score:.4f}', 0, 1, 'L')
        
        # Add radar chart
        if scores["category_scores"]:
            generate_radar_chart(scores["category_scores"], output_dir)
            pdf.image(f"{output_dir}/radar_chart.png", x=10, y=None, w=190)
    
    # Save the report
    pdf_path = f"{output_dir}/{molecule_data['molecular_formula']}_molecular_analysis.pdf"
    pdf.output(pdf_path)
    
    # Create CSV report
    csv_data = flatten_dict(molecule_data)
    pd.DataFrame([csv_data]).to_csv(f"{output_dir}/{molecule_data['molecular_formula']}_molecular_analysis.csv", index=False)
    
    return pdf_path, f"{output_dir}/{molecule_data['molecular_formula']}_molecular_analysis.csv"

def flatten_dict(d: Dict, prefix: str = ""):
    flat_data = {}
    for key, value in d.items():
        if isinstance(value, dict):
            flat_data.update(flatten_dict(value, f"{prefix}{key}_"))
        else:
            flat_data[f"{prefix}{key}"] = value
    return flat_data

if __name__ == "__main__":
    # Example usage
    example_data = {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "molecular_formula": "C9H8O4",
        # Add other data...
    }
    
    output_dir = "reports"
    pdf_path, csv_path = generate_enhanced_report_v2(example_data, output_dir)
    print(f"Generated reports:\nPDF: {pdf_path}\nCSV: {csv_path}")
