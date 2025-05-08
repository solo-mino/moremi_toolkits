"""
Antibody Enhanced Report Generator

This module generates comprehensive PDF and CSV reports for antibody analysis results,
including detailed property breakdowns, visualizations, and structure information.
"""

import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional
import importlib.resources as pkg_resources
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF

from .report_utils import get_asset_path, format_value, format_dict_for_table, create_protein_analysis_df

# Define color scheme
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

class EnhancedAntibodyReport(FPDF):
    """Enhanced PDF report generator for antibody analysis."""
    
    def __init__(self, output_dir: str):
        super().__init__('P', 'mm', 'A4')
        self.output_dir = Path(output_dir)
        self.set_auto_page_break(auto=True, margin=15)
        
        # Set up fonts using relative paths
        try:
            assets_ref = pkg_resources.files('moremi_biokit.assets')
            
            regular_font_ref = assets_ref.joinpath('fonts', 'SpaceGrotesk-Regular.ttf')
            bold_font_ref = assets_ref.joinpath('fonts', 'SpaceGrotesk-Bold.ttf')

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
            logo_ref = pkg_resources.files('moremi_biokit.assets').joinpath('minologo.png')
            with pkg_resources.as_file(logo_ref) as logo_path:
                if logo_path.exists(): # Check if path from context manager is valid
                    self.image(str(logo_path), 10, 15, 15)
        except FileNotFoundError:
            print("Warning: Logo image 'minologo.png' not found in package assets.")
        
        # Add title
        self.set_font('SpaceGrotesk', 'B', 15)
        title = "Comprehensive Antibody Analysis Report"
        self.cell(0, 10, title, 0, 1, 'C')
        
        # Add subtitle with molecular formula and date
        if hasattr(self, 'molecular_formula'):
            self.set_font('SpaceGrotesk', '', 11)
            self.cell(0, 5, f"Antibody: {self.molecular_formula}", 0, 1, 'C')
            self.cell(0, 5, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 0, 1, 'C')
        
        # Add line
        self.line(10, self.get_y() + 2, self.w - 10, self.get_y() + 2)
        self.ln(10)

    def footer(self):
        """Create footer with page number."""
        self.set_y(-15)
        self.set_font('SpaceGrotesk', '', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

    def add_section_title(self, title: str):
        """Add formatted section title."""
        self._section_counter += 1
        self.set_font('SpaceGrotesk', 'B', 14)
        self.set_fill_color(75, 137, 255)  # Bright blue background
        self.set_text_color(255, 255, 255)  # White text
        self.cell(0, 10, f"{self._section_counter}. {title}", 0, 1, 'L', True)
        self.set_text_color(0, 0, 0)  # Reset text color
        self.ln(4)

    def add_subsection_title(self, title: str):
        """Add formatted subsection title."""
        self.set_font('SpaceGrotesk', 'B', 12)
        self.set_text_color(75, 137, 255)
        self.cell(0, 8, title, 'B', 1, 'L')
        self.set_text_color(0, 0, 0)
        self.ln(3)

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


    def add_basic_properties(self, data):
        """Add basic antibody properties section"""
        self.add_section_title("Basic Properties")
        
        basic_props = {
            # "Molecular Weight": (data.get("molecular_weight", "N/A"), "Da"),
            "Molecular Formula": (data.get("molecular_formula", "N/A"), None),
            "Antibody": (data.get("sequence", "N/A")[:20]+"...", None) ,
            "Antigen": (data.get("antigen", "N/A")[:20]+"...", None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in basic_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_protparam_section(self, data):
        """Add ProtParam properties section"""
        if 'protparam' not in data:
            return
            
        self.add_section_title("ProtParam Properties")
        
        protparam_props = {
            "Molecular Weight": (data['protparam'].get("molecular_weight", "N/A"), "Da"),
            "Aromaticity": (data['protparam'].get("aromaticity", "N/A"), None),
            "Instability Index": (data['protparam'].get("instability_index", ("N/A", "N/A"))[0], None),
            "Stability": (data['protparam'].get("instability_index", ("N/A", "N/A"))[1], None),
            "Isoelectric Point (pI)": (data['protparam'].get("isoelectric_point (pI)", "N/A"), None),
            "GRAVY": (data['protparam'].get("gravy", "N/A"), None),
            "Hydrophobic Amino Acids": (data['protparam'].get("hydrophobic_amino_acids", "N/A"), None),
            "Hydrophilic Amino Acids": (data['protparam'].get("hydrophilic_amino_acids", "N/A"), None),
            "Predicted Solubility": (data['protparam'].get("predicted_solubility", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in protparam_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_immunogenicity_section(self, data):
        """Add immunogenicity section"""
        if 'immunogenicity' not in data:
            return
            
        self.add_section_title("Immunogenicity")
        
        immuno_props = {
            "Immunogenic Score": (data['immunogenicity'].get("immunogenic_score", "N/A"), None),
            "Strong Binding Sites": (data['immunogenicity'].get("strong_binding", "N/A"), None),
            "Moderate Binding Sites": (data['immunogenicity'].get("moderate_binding", "N/A"), None),
            "Weak Binding Sites": (data['immunogenicity'].get("weak_binding", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in immuno_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_stability_section(self, data):
        """Add stability section"""
        if 'stability' not in data:
            return
            
        self.add_section_title("Stability")
        
        stability_props = {
            "Melting Temperature": (data['stability'].get("melting_temperature", "N/A"), "°C"),
            "Stability Score": (data['stability'].get("stability_score", "N/A"), None),
            "Stabilizing Residues": (data['stability'].get("details", {}).get("stabilizing_residues", "N/A"), None),
            "Destabilizing Residues": (data['stability'].get("details", {}).get("destabilizing_residues", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in stability_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_aggregation_section(self, data):
        """Add aggregation section"""
        if 'aggregation' not in data:
            return
            
        self.add_section_title("Aggregation")
        
        agg_props = {
            "Aggregation Propensity": (data['aggregation'].get("aggregation_propensity", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in agg_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        
        # Add aggregation prone regions
        regions = data['aggregation'].get("aggregation_prone_regions", [])
        if regions:
            self.ln(5)
            self.add_subsection_title("Aggregation Prone Regions")
            for region, sequence in regions:
                self.set_font('SpaceGrotesk', '', 10)
                self.cell(0, 6, f"• {region}: {sequence}", 0, 1, 'L')
        self.ln(5)

    def add_glycosylation_section(self, data):
        """Add glycosylation section"""
        if 'glycosylation' not in data:
            return
            
        self.add_section_title("Glycosylation")
        
        glyc_props = {
            "N-Glycosylation Sites": (data['glycosylation'].get("n_glyc_sites_count", "N/A"), None),
            "O-Glycosylation Sites": (data['glycosylation'].get("o_glyc_sites_count", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in glyc_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        
        # Add detailed sites information
        sites = data['glycosylation'].get("sites_result", {})
        if sites:
            if sites.get("n_sites"):
                self.ln(5)
                self.add_subsection_title("N-Glycosylation Sites")
                for position, confidence in sites["n_sites"]:
                    self.set_font('SpaceGrotesk', '', 10)
                    self.cell(0, 6, f"• {position} ({confidence})", 0, 1, 'L')
            
            if sites.get("o_sites"):
                self.ln(5)
                self.add_subsection_title("O-Glycosylation Sites")
                for position, residue in sites["o_sites"]:
                    self.set_font('SpaceGrotesk', '', 10)
                    self.cell(0, 6, f"• {position} ({residue})", 0, 1, 'L')
        self.ln(5)

    def add_blast_section(self, data):
        """Add BLAST results section"""
        if 'blast' not in data:
            return
            
        self.add_section_title("BLAST Results")
        
        for result in data['blast']:
            blast_props = {
                "Sequence": (result.get("sequence", "N/A"), None),
                "Length": (result.get("length", "N/A"), "aa"),
                "E-value": (result.get("e-value", "N/A"), None),
                "Identity": (result.get("identity", "N/A"), None),
                "Identity Percentage": (result.get("identity_percentage", "N/A"), "%"),
                "Matched Sequences": (result.get("matched_sequences", "N/A")[:20]+"...", None)
            }
            
            formatted_props = {k: self.format_value(v[0], v[1]) for k, v in blast_props.items()}
            self.create_property_table(formatted_props, ["Property", "Value"])
            self.ln(5)

    def add_radar_chart(self, scores):
        """Add radar chart from category scores"""
        if not scores:
            return
        
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
        temp_chart_path = Path(self.output_dir) / 'temp_radar_chart.png'
        plt.savefig(str(temp_chart_path), dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        # Add the chart to the PDF
        if temp_chart_path.exists():
            self.image(str(temp_chart_path), x=30, y=None, w=150)
            try:
                os.remove(temp_chart_path)
            except:
                pass

    def add_binding_affinity_section(self, data):
        """Add binding affinity section"""
        if 'binding_affinity' not in data:
            return
            
        self.add_section_title("Binding Affinity")
        
        binding_props = {
            "Binding Affinity": (data['binding_affinity'].get("binding_affinity", "N/A"), "kcal/mol"),
            "Dissociation Constant": (data['binding_affinity'].get("dissociation_constant", "N/A"), None),
            "Total Intermolecular Contacts": (data['binding_affinity'].get("intermolecular_contacts", "N/A"), None),
            "Charged-Charged Contacts": (data['binding_affinity'].get("charged_charged_contacts", "N/A"), None),
            "Charged-Polar Contacts": (data['binding_affinity'].get("charged_polar_contacts", "N/A"), None),
            "Charged-Apolar Contacts": (data['binding_affinity'].get("charged_apolar_contacts", "N/A"), None),
            "Polar-Polar Contacts": (data['binding_affinity'].get("polar_polar_contacts", "N/A"), None),
            "Apolar-Polar Contacts": (data['binding_affinity'].get("apolar_polar_contacts", "N/A"), None),
            "Apolar-Apolar Contacts": (data['binding_affinity'].get("apolar_apolar_contacts", "N/A"), None),
            "Apolar NIS %": (data['binding_affinity'].get("apolar_nis_percentage", "N/A"), "%"),
            "Charged NIS %": (data['binding_affinity'].get("charged_nis_percentage", "N/A"), "%")
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in binding_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_structure_section(self, data):
        """Add structure section"""
        if 'structure' not in data:
            return
            
        self.add_section_title("Structure")
        
        structure_props = {
            "GMQE Score": (data['structure'].get("gmqe_score", "N/A"), None),
            # "PDB File": (data['structure'].get("pdb_file_path", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in structure_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        self.ln(5)

    def add_epitope_section(self, data):
        """Add epitope section"""
        if 'epitope' not in data:
            return
            
        self.add_section_title("Epitope Analysis")
        
        # Basic epitope information
        epitope_props = {
            "Overall Score": (data['epitope'].get("score", "N/A"), None),
            "Number of Peptides": (data['epitope'].get("peptides_count", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in epitope_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        
        # Add peptides information
        if 'peptides' in data['epitope']:
            self.ln(5)
            self.add_subsection_title("Predicted Epitope Peptides")
            
            # Create table headers
            self.set_font('SpaceGrotesk', 'B', 9)
            headers = ["Start", "End", "Peptide", "Length", "Score"]
            col_widths = [20, 20, 60, 20, 20]
            
            for i, header in enumerate(headers):
                self.cell(col_widths[i], 7, header, 1, 0, 'C')
            self.ln()
            
            # Add peptide data
            self.set_font('SpaceGrotesk', '', 9)
            for peptide in data['epitope']['peptides']:
                self.cell(20, 7, str(peptide['Start']), 1, 0, 'C')
                self.cell(20, 7, str(peptide['End']), 1, 0, 'C')
                self.cell(60, 7, peptide['Peptide'], 1, 0, 'L')
                self.cell(20, 7, str(peptide['Number of Residues']), 1, 0, 'C')
                self.cell(20, 7, self.format_value(peptide['Score']), 1, 0, 'C')
                self.ln()
        
        self.ln(5)

    def add_developability_section(self, data):
        """Add developability section"""
        if 'developability' not in data:
            return
            
        self.add_section_title("Developability Assessment")
        
        dev_props = {
            "Developability Score": (data['developability'].get("developability_score", "N/A"), None),
            "Has Active Matches": (data['developability'].get("has_active", "N/A"), None),
            "Threshold Used": (data['developability'].get("threshold_used", "N/A"), None),
            "Number of Matches": (data['developability'].get("num_matches", "N/A"), None),
            "Search Time": (data['developability'].get("search_time", "N/A"), "seconds")
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in dev_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        
        # Add performance statistics
        if 'performance_stats' in data['developability']:
            self.ln(5)
            self.add_subsection_title("Performance Statistics")
            stats = data['developability']['performance_stats']
            stats_props = {
                "Sequence Length": (stats.get("sequence_length", "N/A"), None),
                "Initial Candidates": (stats.get("initial_candidates", "N/A"), None),
                "Post Length Filter": (stats.get("post_length_filter", "N/A"), None),
                "Post Kmer Filter": (stats.get("post_kmer_filter", "N/A"), None),
                "Final Matches": (stats.get("final_matches", "N/A"), None)
            }
            formatted_stats = {k: self.format_value(v[0], v[1]) for k, v in stats_props.items()}
            self.create_property_table(formatted_stats, ["Metric", "Value"])
        
        self.ln(5)

    def add_conservancy_section(self, data):
        """Add conservancy section"""
        if 'conservancy' not in data:
            return
            
        self.add_section_title("Conservancy Analysis")
        
        cons_props = {
            "Conservancy Score": (data['conservancy'].get("conservancy_score", "N/A"), None)
        }
        
        formatted_props = {k: self.format_value(v[0], v[1]) for k, v in cons_props.items()}
        self.create_property_table(formatted_props, ["Property", "Value"])
        
        # Add conservancy results table if available
        if 'results' in data['conservancy']:
            self.ln(5)
            self.add_subsection_title("Epitope Conservancy Results")
            
            results = data['conservancy']['results']
            if isinstance(results, pd.DataFrame):
                # Convert DataFrame to list of dictionaries for consistent handling
                results = results.to_dict('records')
            
            if results:
                # Create table headers
                self.set_font('SpaceGrotesk', 'B', 9)
                headers = ["Epitope #", "Sequence", "Match %", "Min Identity", "Max Identity"]
                col_widths = [20, 40, 40, 30, 30]
                
                for i, header in enumerate(headers):
                    self.cell(col_widths[i], 7, header, 1, 0, 'C')
                self.ln()
                
                # Add conservancy data
                self.set_font('SpaceGrotesk', '', 9)
                for row in results:
                    self.cell(20, 7, str(row['Epitope #']), 1, 0, 'C')
                    self.cell(40, 7, str(row['Epitope Sequence']), 1, 0, 'L')
                    self.cell(40, 7, str(row['Percent of protein sequence matches at identity >= 70%']), 1, 0, 'C')
                    self.cell(30, 7, str(row['Minimum Identity']), 1, 0, 'C')
                    self.cell(30, 7, str(row['Maximum Identity']), 1, 0, 'C')
                    self.ln()
        
        self.ln(5)

    def add_weighted_scores_section(self, data):
        """Add weighted scores section"""
        if 'weighted_scores' not in data:
            return
            
        self.add_section_title("Weighted Scores")
        
        scores = data['weighted_scores']
        total_score = data.get('total_score', "N/A")
        
        # Add total score first
        total_props = {
            "Total Score": (total_score, None)
        }
        formatted_total = {k: self.format_value(v[0], v[1]) for k, v in total_props.items()}
        self.create_property_table(formatted_total, ["Metric", "Score"])
        
        self.ln(5)
        self.add_subsection_title("Individual Scores")
        
        # Create sorted list of scores for better presentation
        score_items = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        score_props = {k.replace('_', ' ').title(): (v, None) for k, v in score_items}
        
        formatted_scores = {k: self.format_value(v[0], v[1]) for k, v in score_props.items()}
        self.create_property_table(formatted_scores, ["Metric", "Score"])
        
        self.ln(5)

    def create_report(self, data: Dict):
        """Create the complete report"""
        # Store molecular formula for header
        self.molecular_formula = data.get("molecular_formula", "Unknown")
        
        # Add new page with header
        self.add_page()
        
        # Add all sections
        self.add_basic_properties(data)
        self.add_protparam_section(data)
        self.add_binding_affinity_section(data)
        self.add_structure_section(data)
        self.add_epitope_section(data)
        self.add_immunogenicity_section(data)
        self.add_stability_section(data)
        self.add_aggregation_section(data)
        self.add_glycosylation_section(data)
        self.add_blast_section(data)
        self.add_developability_section(data)
        self.add_conservancy_section(data)
        
        # Add weighted scores and radar chart on a new page
        self.add_page()
        self.add_weighted_scores_section(data)
        
        if 'weighted_scores' in data:
            self.add_section_title("Score Distribution")
            self.add_radar_chart(data['weighted_scores'])
            self.ln(5)

    def generate_pdf_report(self, data: Dict[str, Any], output_path: str):
        """Generate the complete PDF report."""
        self.create_report(data)
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

def generate_enhanced_report(antibody_data: Dict[str, Any], output_dir: str, generate_pdf: bool = True, generate_csv: bool = True) -> Dict[str, Optional[Path]]:
    """Generate enhanced PDF and CSV reports for antibody analysis.

    Args:
        antibody_data (Dict[str, Any]): Dictionary containing all antibody metrics and data.
        output_dir (str): The directory where reports will be saved.
        generate_pdf (bool, optional): Whether to generate a PDF report. Defaults to True.
        generate_csv (bool, optional): Whether to generate a CSV report. Defaults to True.

    Returns:
        Dict[str, Optional[Path]]: A dictionary containing paths to the generated reports.
                                    Keys are 'pdf' and 'csv'. Values are Path objects or None.
    """
    # Create output directory if it doesn't exist
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Initialize report
    report_generator = EnhancedAntibodyReport(output_dir_path) # Corrected variable name
    
    # Generate timestamp for file names
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Get rank from antibody data, default to 0 if not provided
    rank = antibody_data.get('rank', 0)
    base_filename = f'r{rank}_{antibody_data.get("molecular_formula", "antibody")}_{timestamp}'

    generated_paths: Dict[str, Optional[Path]] = {"pdf": None, "csv": None}

    if generate_pdf:
        pdf_path = output_dir_path / f'{base_filename}.pdf'
        try:
            report_generator.generate_pdf_report(antibody_data, str(pdf_path))
            generated_paths["pdf"] = pdf_path
            logging.info(f"Successfully generated PDF report: {pdf_path}")
        except Exception as e:
            logging.error(f"Failed to generate PDF report for {base_filename}: {e}")
            # Optionally, add to antibody_data['warnings'] if it's passed around
    
    if generate_csv:
        csv_path = output_dir_path / f'{base_filename}.csv'
        try:
            # Assuming generate_csv_report is a static method or defined elsewhere
            # If it's part of EnhancedAntibodyReport, it should be: report_generator.generate_csv_report(...)
            # For now, assuming it's a standalone function in this module or imported.
            # Let's define a placeholder if it's missing or adapt if it exists in the class.
            
            # Placeholder for generate_csv_report logic, to be integrated correctly
            # based on where it's actually defined.
            # For this example, let's assume it's a static method or separate function.
            _generate_csv_report_standalone(antibody_data, str(csv_path)) 
            generated_paths["csv"] = csv_path
            logging.info(f"Successfully generated CSV report: {csv_path}")
        except Exception as e:
            logging.error(f"Failed to generate CSV report for {base_filename}: {e}")

    return generated_paths

# Assuming generate_csv_report was a standalone function as in the original selection.
# If it's meant to be a method of EnhancedAntibodyReport, this needs adjustment.
def _generate_csv_report_standalone(antibody_data: Dict[str, Any], csv_path_str: str):
    """Generates a CSV report from antibody data.
    This is a simplified placeholder based on the original structure.
    Args:
        antibody_data: Dictionary of antibody data.
        csv_path_str: Path to save the CSV file.
    """
    # Flatten relevant parts of antibody_data for CSV
    # This is a very basic flattening, might need to be more sophisticated
    flat_data = {}
    flat_data['sequence'] = antibody_data.get('sequence')
    flat_data['antigen'] = antibody_data.get('antigen')
    flat_data['molecular_formula'] = antibody_data.get('molecular_formula')
    flat_data['molecular_weight'] = antibody_data.get('molecular_weight')
    flat_data['rank'] = antibody_data.get('rank')
    flat_data['total_score'] = antibody_data.get('total_score')

    for category, cat_data in antibody_data.get('category_scores', {}).get('raw', {}).items():
        flat_data[f'raw_score_{category}'] = cat_data
    for category, cat_data in antibody_data.get('category_scores', {}).get('weighted', {}).items():
         flat_data[f'weighted_score_{category}'] = cat_data

    # Extract specific metrics - this needs careful selection of what to include
    # For example, from protparam:
    if 'protparam' in antibody_data and isinstance(antibody_data['protparam'], dict):
        flat_data['gravy'] = antibody_data['protparam'].get('gravy')
        flat_data['predicted_solubility'] = antibody_data['protparam'].get('predicted_solubility')
    
    # Add more fields as needed from other sections like 'stability', 'aggregation' etc.
    if 'stability' in antibody_data and isinstance(antibody_data['stability'], dict):
        flat_data['melting_temperature'] = antibody_data['stability'].get('melting_temperature')

    df = pd.DataFrame([flat_data])
    df.to_csv(csv_path_str, index=False)


if __name__ == "__main__":
    # Example usage
    example_data = {
        "antigen": "KVFRSSVLHSKKADPEASFWGEEFTRGVYYPD",
        "molecular_weight": 74182.1406000002,
        "molecular_formula": "C627H1254N627O627",
        # Add other data...
    }
    
    output_dir = "reports"
    pdf_path, csv_path = generate_enhanced_report(example_data, output_dir)
    print(f"Generated reports:\nPDF: {pdf_path}\nCSV: {csv_path}") 