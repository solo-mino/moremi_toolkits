"""
Enhanced Molecule Ranking Report Generator

This module generates comprehensive PDF reports for molecule rankings including validation metrics,
category scores, and data visualizations in a professional format.
"""

import pandas as pd
from fpdf import FPDF
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
from math import ceil
import importlib.resources as pkg_resources

class MoleculeReportPDF(FPDF):
    def __init__(self):
        super().__init__('L')  # Landscape mode for wider tables
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
            
        self._current_section = 0
        
        # Set default font
        self.set_font('SpaceGrotesk', '', 10)

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
        
        # Add main title
        self.set_font('SpaceGrotesk', 'B', 16)
        self.cell(0, 10, 'Molecular Analysis and Ranking Report', 0, 1, 'C')
        
        # Add generated timestamp
        self.set_font('SpaceGrotesk', '', 10)
        self.cell(0, 5, f'Generated on: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', 0, 1, 'C')
        self.ln(10)  # Add some space after the header

    def footer(self):
        self.set_y(-15)
        self.set_font('SpaceGrotesk', '', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

    def section_title(self, title):
        """Add a section title with consistent formatting"""
        self.set_font('SpaceGrotesk', 'B', 12)
        self.set_fill_color(230, 230, 230)
        self.cell(0, 10, f'{title}', 0, 1, 'L', True)
        self.ln(5)

    def create_wrapped_cell(self, w, h, txt, border=0, align='L', fill=False, padding=0):
        """Create a cell with wrapped text and optional padding"""
        # Store current margins
        cw = w - (2 * padding)  # Reduce width by padding
        
        # Handle newlines and wrapping
        if '\n' in txt:
            lines = txt.split('\n')
        else:
            lines = self.multi_cell(cw, 5, txt, 0, align, fill, split_only=True)
        
        # Get height of all lines
        height = h
        
        # Draw Cell
        x_start = self.get_x()
        y_start = self.get_y()
        
        # Draw the border and fill first
        if fill:
            self.rect(x_start, y_start, w, height, 'F')
        if border:
            self.rect(x_start, y_start, w, height)
        
        # Calculate vertical centering
        total_line_height = len(lines) * 5
        y_offset = (height - total_line_height) / 2
        
        # Print each line
        self.set_xy(x_start + padding, y_start + y_offset)
        for line in lines:
            self.cell(cw, 5, line, 0, 0, align)
            self.set_xy(x_start + padding, self.get_y() + 5)
        
        # Reset position to the right of the cell
        self.set_xy(x_start + w, y_start)

    def create_table_header(self, headers, col_widths):
        """Create table header with improved formatting"""
        self._last_headers = headers.copy()
        self._last_col_widths = col_widths.copy()
        
        # Keep headers as single line
        formatted_headers = []
        for header in headers:
            # Remove any newlines and keep headers short
            formatted_header = header.strip()
            if header == 'Medicinal Chemistry':
                formatted_header = 'Med Chem'
            elif header == 'Physicochemical':
                formatted_header = 'Phys Chem'
            formatted_headers.append(formatted_header)
        
        self.set_font('SpaceGrotesk', 'B', 9)  # Slightly smaller font for headers
        self.set_fill_color(66, 133, 244)  # Blue header
        self.set_text_color(255, 255, 255)  # White text
        
        # Fixed height for header row
        header_height = 8  # Reduced height for single line
        
        # Draw header cells
        for header, width in zip(formatted_headers, col_widths):
            x_start = self.get_x()
            y_start = self.get_y()
            
            # Draw background
            self.rect(x_start, y_start, width, header_height, 'F')
            
            # Draw border
            self.rect(x_start, y_start, width, header_height)
            
            # Center text
            self.set_xy(x_start, y_start + (header_height - 4) / 2)
            self.cell(width, 4, header, 0, 0, 'C')
            
            # Move to next cell position
            self.set_xy(x_start + width, y_start)
        
        self.ln(header_height)
        self.set_font('SpaceGrotesk', '', 10)  # Reset font for content
        self.set_text_color(0, 0, 0)  # Reset to black text

    def create_header_cell(self, w, h, txt, border=0):
        """Create a header cell with proper alignment"""
        x_start = self.get_x()
        y_start = self.get_y()
        
        # Draw background first
        self.rect(x_start, y_start, w, h, 'F')
        
        # Draw border
        if border:
            self.rect(x_start, y_start, w, h)
        
        # Handle multi-line text
        if '\n' in txt:
            lines = txt.split('\n')
            line_height = h / len(lines)
            for i, line in enumerate(lines):
                y = y_start + (i * line_height)
                self.set_xy(x_start, y)
                self.cell(w, line_height, line, 0, 0, 'C')
        else:
            # Center single line text
            self.set_xy(x_start, y_start + (h - 5) / 2)  # 5 is approximate text height
            self.cell(w, 5, txt, 0, 0, 'C')
        
        # Move to next cell position
        self.set_xy(x_start + w, y_start)

    def create_table_row(self, data, col_widths, is_alternating=True):
        """Create table row with wrapped text and alternating colors"""
        if is_alternating and self.get_y() % 20 < 10:
            self.set_fill_color(245, 245, 245)  # Light gray
            fill = True
        else:
            self.set_fill_color(255, 255, 255)  # White
            fill = False
        
        self.set_text_color(0, 0, 0)
        
        # Calculate row height based on content
        max_height = 7  # Reduced minimum height
        for content, width in zip(data, col_widths):
            content_str = str(content)
            # Add extra height for SMILES column but with reduced multiplier
            if len(content_str) > 30:  # SMILES string
                needed_height = self.get_string_height(width - 4, content_str) * 1.2  # Reduced multiplier
            else:
                needed_height = self.get_string_height(width - 4, content_str)
            max_height = max(max_height, needed_height + 2)  # Reduced padding
        
        # Check for page break
        self.check_page_break(max_height)
        
        # Print the row with consistent spacing
        x_start = self.get_x()
        y_start = self.get_y()
        
        for content, width in zip(data, col_widths):
            self.set_xy(x_start, y_start)
            # Added cell padding of 2 points
            self.create_wrapped_cell(width, max_height, str(content), 1, 'C', fill, 2)
            x_start += width
        
        self.ln(max_height)

    def check_page_break(self, height):
        """Check if adding content of given height would exceed page boundary"""
        if self.get_y() + height > self.page_break_trigger:
            self.add_page()
            self.create_table_header(self._last_headers, self._last_col_widths)
            
    def get_string_height(self, width, txt):
        """Calculate the height needed for wrapped text"""
        cw = self.get_string_width(txt)
        if cw == 0:
            return 7
        nb_lines = max(1, ceil(cw / width))
        return 7 * nb_lines

def calculate_mol_weight(smiles):
    """Calculate molecular weight from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Descriptors.ExactMolWt(mol)
    except:
        pass
    return 0.0

def create_validation_tables(pdf, df):
    """Create multiple validation tables to show all parameters"""
    
    # Table 1: Basic Properties
    pdf.add_page()
    pdf.section_title('Basic Molecular Properties')
    headers = ['SMILES', 'Molecular Formula', 'Molecular Weight', 'TPSA', 'LogP']
    widths = [90, 40, 30, 30, 30]  # Reduced SMILES width
    pdf.create_table_header(headers, widths)
    for _, row in df.iterrows():
        mol_weight = calculate_mol_weight(row['SMILES'])  # Calculate molecular weight from SMILES
        data = [
            row['SMILES'],
            row['Molecular Formula'],
            format_numeric(mol_weight),
            format_numeric(row.get('tpsa value', 0.0)),
            format_numeric(row.get('logp value', 0.0))
        ]
        pdf.create_table_row(data, widths)
    
    # Table 2: Lipophilicity Parameters
    pdf.add_page()
    pdf.section_title('Lipophilicity Metrics')
    headers = ['SMILES', 'iLogP', 'xLogP3', 'wLogP', 'mLogP', 'SILICOS-IT','Consensus LogP']
    widths = [90, 30, 30, 30, 30, 30, 30]
    pdf.create_table_header(headers, widths)
    for _, row in df.iterrows():
        data = [
            row['SMILES'],
            format_numeric(row.get('ilogp value', 0.0)),
            format_numeric(row.get('xlogp3 value', 0.0)),
            format_numeric(row.get('wlogp value', 0.0)),
            format_numeric(row.get('mlogp value', 0.0)),
            format_numeric(row.get('silicos_it value', 0.0)),
            format_numeric(row.get('consensus_logp value', 0.0))
        ]
        pdf.create_table_row(data, widths)
    
    # Table 3: Druglikeness and Medicinal Chemistry
    pdf.add_page()
    pdf.section_title('Druglikeness and Medicinal Chemistry Metrics')
    headers = ['SMILES', 'QED', 'Synthetic Accessibility', 'Lipinski Violations', 'Bioavailability']
    widths = [90, 30, 35, 35, 35]  # Reduced SMILES width
    pdf.create_table_header(headers, widths)
    for _, row in df.iterrows():
        data = [
            row['SMILES'],
            format_numeric(row.get('qed value', 0.0), 3),
            format_numeric(row.get('synthetic_accessibility value', 0.0)),
            format_numeric(row.get('lipinski_violations value', 0), 0),
            format_numeric(row.get('bioavailability_score value', 0.0))
        ]
        pdf.create_table_row(data, widths)
    
    # Table 4: ADMETlab Absorption
    pdf.add_page()
    pdf.section_title('ADMETlab Absorption Metrics')
    headers = ['SMILES', 'Caco2', 'PAMPA', 'MDCK', 'HIA', 'PGP-sub', 'PGP-inh']
    widths = [90, 25, 25, 25, 25, 25, 25]  # Reduced SMILES width
    pdf.create_table_header(headers, widths)
    for _, row in df.iterrows():
        data = [
            row['SMILES'],
            format_numeric(row.get('caco2 value', 0.0)),
            format_numeric(row.get('pampa value', 0.0)),
            format_numeric(row.get('mdck value', 0.0)),
            format_numeric(row.get('hia value', 0.0)),
            format_numeric(row.get('pgp_substrate value', 0.0)),
            format_numeric(row.get('pgp_inhibitor value', 0.0))
        ]
        pdf.create_table_row(data, widths)
    
    # Table 5: Distribution and Toxicity
    pdf.add_page()
    pdf.section_title('Distribution and Toxicity Metrics')
    headers = ['SMILES', 'VDss', 'PPB', 'BBB', 'Fu', 'HERG']
    widths = [90, 30, 30, 30, 30, 30]  # Reduced SMILES width
    pdf.create_table_header(headers, widths)
    for _, row in df.iterrows():
        data = [
            row['SMILES'],
            format_numeric(row.get('vdss value', 0.0)),
            format_numeric(row.get('ppb value', 0.0)),
            format_numeric(row.get('bbb value', 0.0)),
            format_numeric(row.get('fu value', 0.0)),
            format_numeric(row.get('herg value', 0.0))
        ]
        pdf.create_table_row(data, widths)
    
    # Table 6: Metabolism and Excretion
    pdf.add_page()
    pdf.section_title('Metabolism and Excretion Metrics')
    headers = ['SMILES', 'CYP2C9', 'CYP2D6', 'CYP3A4', 'Half Life', 'Clearance']
    widths = [90, 30, 30, 30, 30, 30]  # Reduced SMILES width
    pdf.create_table_header(headers, widths)
    for _, row in df.iterrows():
        data = [
            row['SMILES'],
            format_numeric(row.get('cyp2c9_inhibition value', 0.0)),
            format_numeric(row.get('cyp2d6_inhibition value', 0.0)),
            format_numeric(row.get('cyp3a4_inhibition value', 0.0)),
            format_numeric(row.get('half_life value', 0.0)),
            format_numeric(row.get('clearance value', 0.0))
        ]
        pdf.create_table_row(data, widths)

def create_visualizations(pdf, df, top_n=None):
    """Create data visualizations"""
    # First visualization - Violin Plot
    pdf.add_page()
    pdf.section_title('Score Distributions')
    
    # Create violin plot with adjusted dimensions
    plt.figure(figsize=(17, 8))
    
    # Get score columns and prepare data
    score_cols = [col for col in df.columns if col.endswith(' Score') and col != 'Lipophilicity Score']
    
    if top_n is not None and top_n < len(df):
        df = df.head(top_n)
    
    data = [df[col].values for col in score_cols]
    labels = [col.replace(' Score', '') for col in score_cols]
    
    # Set color palette
    colors = ['#2ecc71', '#3498db', '#9b59b6', '#e74c3c', '#f1c40f', 
              '#1abc9c', '#e67e22', '#34495e', '#2980b9', '#8e44ad', '#c0392b']
    
    # Create violin plot with seaborn
    sns.violinplot(data=data, palette=colors, alpha=0.7)
    
    # Customize plot appearance
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.ylabel('Score', fontsize=12)
    plt.title('Distribution of Scores Across Categories', fontsize=14, pad=20)
    
    # Add grid for better readability
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    
    # Calculate y-axis limits with dynamic padding for both upper and lower bounds
    max_score = max(max(series) for series in data)
    min_score = min(min(series) for series in data)
    
    upper_padding = max(0.5, (1.5 - max_score))  # At least 0.5 padding, more if scores are low
    lower_padding = max(0.5, (min_score + 0.5))  # Similar padding for lower bound
    
    # Set y-axis limits with dynamic padding
    plt.ylim(min(-0.5, min_score - lower_padding), min(2.0, max_score + upper_padding))
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot to temporary file
    temp_plot = 'temp_violin_plot.png'
    plt.savefig(temp_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Add plot to PDF
    pdf.image(temp_plot, x=10, y=None, w=270)
    
    # Clean up
    os.remove(temp_plot)
    
    # Second visualization - Statistical Analysis Plot
    pdf.add_page()
    pdf.section_title('Statistical Analysis by Category Plot')
    
    # Calculate statistics for each category
    categories = score_cols
    means = [df[col].mean() for col in categories]
    stds = [df[col].std() for col in categories]
    mins = [df[col].min() for col in categories]
    maxs = [df[col].max() for col in categories]
    
    # Create figure and axis with larger size
    plt.figure(figsize=(17, 8))
    
    # Create positions for bars
    x = np.arange(len(categories))
    
    # Create bars
    plt.bar(x, means, yerr=stds, capsize=5, color='royalblue', alpha=0.7, 
            label='Mean ± Std Dev')
    
    # Add min/max points
    plt.scatter(x, mins, color='red', marker='v', s=100, label='Minimum')
    plt.scatter(x, maxs, color='green', marker='^', s=100, label='Maximum')
    
    # Customize the plot
    plt.xlabel('Categories', fontsize=12)
    plt.ylabel('Values', fontsize=12)
    plt.title('Molecular Analysis Statistics by Category', fontsize=14, pad=20)
    
    # Format category labels
    category_labels = [col.replace(' Score', '') for col in categories]
    plt.xticks(x, category_labels, rotation=45, ha='right')
    
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend()
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Add value labels on top of bars
    for i, v in enumerate(means):
        plt.text(i, v + stds[i] + 0.02, f'{v:.3f}', ha='center', fontsize=9)
    
    # Set y-axis limits with some padding
    plt.ylim(0, max(maxs) + 0.1)
    
    # Save plot to temporary file
    temp_plot = 'temp_stats_plot.png'
    plt.savefig(temp_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Add plot to PDF
    pdf.image(temp_plot, x=10, y=None, w=270)
    
    # Clean up
    os.remove(temp_plot)

def format_numeric(value, decimal_places=2):
    """Format a value as numeric with specified decimal places"""
    try:
        # Convert to float first to handle both string and numeric inputs
        num_value = float(value) if value is not None else 0.0
        return f"{num_value:.{decimal_places}f}"
    except (ValueError, TypeError):
        return str(value)  # Return as string if conversion fails

def analyze_data(df):
    """Analyze data to generate key findings"""
    findings = []
    
    # Keep only the essential findings
    findings.append(f"Number of molecules analyzed: {len(df)}")
    findings.append(f"Average overall score: {df['Overall Score'].mean():.3f}")
    findings.append(f"Score range: {df['Overall Score'].min():.3f} - {df['Overall Score'].max():.3f}")
    findings.append(f"Top performing molecule: {df.iloc[0]['SMILES'][:50]}")
    
    return findings

def generate_ranking_report(input_csv: str, output_pdf: str, top_n: int = None):
    """Generate a comprehensive PDF report from the ranking results CSV."""
    # Convert paths to Path objects
    input_csv = Path(input_csv)
    output_pdf = Path(output_pdf)
    
    # Read the data
    df = pd.read_csv(input_csv)
    
    # Create PDF
    pdf = MoleculeReportPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    
    # Executive Summary
    pdf.add_page()
    pdf.section_title('Executive Summary')
    pdf.set_font('SpaceGrotesk', '', 10)
    pdf.multi_cell(0, 7, "This comprehensive report presents a detailed analysis of molecule rankings based on multiple evaluation criteria including physicochemical properties, ADMET characteristics, and drug-likeness metrics. The analysis encompasses validation metrics, category-wise scores, and statistical analysis of the results.")
    pdf.ln(10)

    # Key Findings
    pdf.section_title('Key Findings')
    findings = analyze_data(df)
    for finding in findings:
        pdf.set_font('SpaceGrotesk', '', 10)
        pdf.cell(0, 7, "• " + finding, 0, 1)
    pdf.ln(10)

    # Category Scores and Rankings
    pdf.add_page()
    pdf.section_title('Category Scores and Rankings')
        
    # Define column widths for score table - adjusted for better fit
    score_widths = [15, 80]  # Rank, SMILES
    score_cols = [
        'Physicochemical',
        'Medicinal Chemistry',
        'Druglikeness',
        'Absorption',
        'Distribution',
        'Metabolism',
        'Excretion',
        'Toxicity',
        'Overall'
    ]
    # Add equal width for all score columns
    col_width = (pdf.w - sum(score_widths) - 20) / len(score_cols)  # -20 for margins
    score_widths.extend([col_width] * len(score_cols))
    
    
    # Add data rows
    if top_n is not None and top_n < len(df):
        df = df.head(top_n)
        pdf.ln(10)
        pdf.set_font('SpaceGrotesk', 'B', 10)
        pdf.cell(0, 10, f'Showing top {top_n} molecules', 0, 1, 'L')
    
    # Create header row with original column names
    headers = ['Rank', 'SMILES'] + score_cols
    pdf.create_table_header(headers, score_widths)    
        
    for idx, row in df.iterrows():
        data = [
            idx + 1,
            row['SMILES']
        ] + [format_numeric(row.get(col + ' Score', 0.0)) for col in score_cols]
        pdf.create_table_row(data, score_widths)
    
    # Validation Results Tables (multiple tables)
    create_validation_tables(pdf, df)

    # Data Visualizations with increased size
    create_visualizations(pdf, df, top_n)

    # Statistical Analysis
    pdf.add_page()
    pdf.section_title('Table of Statistical Analysis Values for Section 11: Category Plots') 
    stats_headers = ['Category', 'Mean', 'Std Dev', 'Min', 'Max']
    stats_widths = [80, 30, 30, 30, 30]
    pdf.create_table_header(stats_headers, stats_widths)
    
    score_columns = [col for col in df.columns if col.endswith('Score')]
    for col in score_columns:
        stats = df[col].describe()
        category = col.replace(' Score', '')
        data = [
            category,
            format_numeric(stats['mean'], 3),
            format_numeric(stats['std'], 3),
            format_numeric(stats['min'], 3),
            format_numeric(stats['max'], 3)
        ]
        pdf.create_table_row(data, stats_widths, is_alternating=True)

    # Save the PDF
    pdf.output(str(output_pdf))

def generate_report(smiles: str, mol_formula: str, output_dir: Path, timestamp: str):
    """
    Generate a PDF report for a molecule
    """
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filenames using molecular formula
    base_filename = f"{mol_formula}_{timestamp}"
    pdf_path = output_dir / f"{base_filename}.pdf"
    csv_path = output_dir / f"{base_filename}.csv"

if __name__ == "__main__":
    input_csv = "rankings_20250109_004826.csv"
    output_pdf = "enhanced_molecule_ranking_report.pdf"
    generate_ranking_report(input_csv, output_pdf, top_n=10)