"""
Antibody Ranking Report Generator

This module generates comprehensive PDF reports for protein rankings including validation metrics,
category scores, and data visualizations in a professional format.
"""

import pandas as pd
from fpdf import FPDF
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from math import ceil
import importlib.resources as pkg_resources

class ProteinReportPDF(FPDF):
    def __init__(self):
        super().__init__('L')  # Landscape mode for wider tables
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
        self._current_section = 0
        
        # Set default font
        self.set_font('SpaceGrotesk', '', 10)

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
        
        # Add main title
        self.set_font('SpaceGrotesk', 'B', 16)
        self.cell(0, 10, 'Antibody Analysis and Ranking Report', 0, 1, 'C')
        
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
            if header == 'Binding Affinity':
                formatted_header = 'Binding'
            elif header == 'Immunogenicity':
                formatted_header = 'Immuno'
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

def format_numeric(value, decimal_places=2):
    """Format a value as numeric with specified decimal places"""
    try:
        num_value = float(value) if value is not None else 0.0
        return f"{num_value:.{decimal_places}f}"
    except (ValueError, TypeError):
        return str(value)

def analyze_data(df):
    """Analyze data to generate key findings"""
    findings = []
    
    findings.append(f"Number of proteins analyzed: {len(df)}")
    
    if 'total_score' in df.columns:
        findings.append(f"Average overall score: {df['total_score'].mean():.3f}")
        findings.append(f"Score range: {df['total_score'].min():.3f} - {df['total_score'].max():.3f}")
    
    if len(df) > 0:
        findings.append(f"Top performing protein: {df.iloc[0]['antigen'][:50]}")
    
    return findings

def create_visualizations(pdf, df, top_n=None):
    """Create data visualizations"""
    # Create a subset of the data for visualization
    if top_n is not None and top_n < len(df):
        viz_df = df.head(top_n)
    else:
        viz_df = df
    
    # First visualization - Score Distribution
    pdf.add_page()
    pdf.section_title('Score Distributions by Category')
    
    # Get score columns (both raw and weighted)
    raw_score_cols = [col for col in df.columns if col.endswith('_Score') and not col.startswith('Weighted')]
    
    if not raw_score_cols:
        pdf.cell(0, 10, "No score data available for visualization", 0, 1, 'L')
        return
    
    # Create violin plot
    plt.figure(figsize=(12, 6))
    
    # Prepare data for violin plot
    data = []
    labels = []
    
    for col in raw_score_cols:
        if viz_df[col].notna().any():  # Only include columns with data
            data.append(viz_df[col].values)
            labels.append(col.replace('_Score', ''))
    
    if not data:
        pdf.cell(0, 10, "No valid score data available for visualization", 0, 1, 'L')
        return
    
    # Set color palette
    colors = sns.color_palette("viridis", len(data))
    
    # Create violin plot with seaborn
    violin_parts = sns.violinplot(data=data, palette=colors, alpha=0.7)
    
    # Add individual data points for better visibility
    for i, d in enumerate(data):
        plt.scatter([i] * len(d), d, color='white', alpha=0.5, s=5)
    
    # Customize plot appearance
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.ylabel('Score', fontsize=12)
    plt.title('Distribution of Raw Scores Across Categories', fontsize=14)
    
    # Add grid for better readability
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    
    # Ensure y-axis limits are appropriate
    plt.ylim(0, 1.1)
    
    plt.tight_layout()
    
    # Save plot to temporary file
    temp_plot = 'temp_violin_plot.png'
    plt.savefig(temp_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Add plot to PDF
    pdf.image(temp_plot, x=10, y=None, w=270)
    
    # Clean up
    if os.path.exists(temp_plot):
        os.remove(temp_plot)
    
    # Second visualization - Weighted Scores Bar Chart
    pdf.add_page()
    pdf.section_title('Weighted Score Comparison')
    
    # Get weighted score columns
    weighted_score_cols = [col for col in df.columns if col.startswith('Weighted_')]
    
    if not weighted_score_cols:
        pdf.cell(0, 10, "No weighted score data available for visualization", 0, 1, 'L')
        return
    
    # Create a bar chart for weighted scores
    plt.figure(figsize=(12, 6))
    
    # Get the top proteins
    top_proteins = viz_df.head(min(5, len(viz_df)))
    
    # Prepare data
    categories = [col.replace('Weighted_', '') for col in weighted_score_cols]
    protein_names = [f"Antibody {i+1}" for i in range(len(top_proteins))]
    
    # Create data matrix for bar chart
    data_matrix = []
    for _, row in top_proteins.iterrows():
        data_matrix.append([row[col] for col in weighted_score_cols])
    
    # Transpose data for grouped bar chart
    data_matrix = np.array(data_matrix).T
    
    # Set up the bar chart
    bar_width = 0.15
    x = np.arange(len(categories))
    
    # Plot bars for each protein
    for i in range(len(top_proteins)):
        plt.bar(x + i*bar_width - (len(top_proteins)*bar_width/2) + bar_width/2, 
                data_matrix[:, i], width=bar_width, label=protein_names[i])
    
    # Customize plot
    plt.xlabel('Categories', fontsize=12)
    plt.ylabel('Weighted Score', fontsize=12)
    plt.title('Weighted Scores by Category for Top proteins', fontsize=14)
    plt.xticks(x, categories, rotation=45, ha='right')
    plt.legend()
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    temp_plot = 'temp_bar_plot.png'
    plt.savefig(temp_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Add to PDF
    pdf.image(temp_plot, x=10, y=None, w=270)
    
    # Clean up
    if os.path.exists(temp_plot):
        os.remove(temp_plot)
    
    # Third visualization - Radar Chart of Weighted Scores for Top Antibody
    if len(viz_df) > 0:
        pdf.add_page()
        pdf.section_title('Top Antibody Score Profile')
        
        top_protein = viz_df.iloc[0]
        
        # Prepare radar chart data
        categories = [col.replace('_Score', '') for col in raw_score_cols]
        values = [top_protein[col] for col in raw_score_cols]
        
        # Create radar chart
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, polar=True)
        
        # Number of categories
        N = len(categories)
        
        # Angle of each axis
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop
        
        # Values for the chart
        values += values[:1]  # Close the loop
        
        # Draw the polygon
        ax.plot(angles, values, linewidth=1, linestyle='solid')
        ax.fill(angles, values, alpha=0.25)
        
        # Add category labels
        plt.xticks(angles[:-1], categories, size=10)
        
        # Add value grid lines
        ax.set_rlabel_position(0)
        plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0], ["0.2", "0.4", "0.6", "0.8", "1.0"], 
                  color="grey", size=8)
        plt.ylim(0, 1)
        
        plt.title('Score Profile of Top Ranked Antibody', y=1.1, fontsize=14)
        
        # Save and add to PDF
        temp_plot = 'temp_radar_plot.png'
        plt.savefig(temp_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        pdf.image(temp_plot, x=40, y=None, w=200)
        
        # Clean up
        if os.path.exists(temp_plot):
            os.remove(temp_plot)

def generate_ranking_report(input_csv: str, output_pdf: str, top_n: int = None):
    """Generate a comprehensive PDF report from the protein ranking results CSV."""
    # Convert paths to Path objects
    input_csv = Path(input_csv)
    output_pdf = Path(output_pdf)
    
    # Read the data
    df = pd.read_csv(input_csv)
    
    # Create PDF
    pdf = ProteinReportPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    
    # Executive Summary
    pdf.add_page()
    pdf.section_title('Executive Summary')
    pdf.set_font('SpaceGrotesk', '', 10)
    pdf.multi_cell(0, 7, "This comprehensive report presents a detailed analysis of protein rankings based on multiple evaluation criteria including binding affinity, immunogenicity, stability, and developability metrics. The rankings are based on both raw scores and weighted scores that reflect the relative importance of each category.")
    pdf.ln(5)

    # Key Findings
    pdf.section_title('Key Findings')
    findings = analyze_data(df)
    for finding in findings:
        pdf.set_font('SpaceGrotesk', '', 10)
        pdf.cell(0, 7, "• " + finding, 0, 1)
    pdf.ln(5)

    # Category Scores and Rankings
    pdf.add_page()
    pdf.section_title('Antibody Rankings')
        
    # Define column widths for score table
    score_widths = [15, 80]  # Rank, Antigen
    
    # Find all raw score columns
    raw_score_cols = [col for col in df.columns if col.endswith('_Score') and not col.startswith('Weighted')]
    
    # Find all weighted score columns
    weighted_score_cols = [col for col in df.columns if col.startswith('Weighted_')]
    
    # Add equal width for all score columns
    raw_col_width = 25
    score_widths.extend([raw_col_width] * (len(raw_score_cols) + 1))  # +1 for total score
        
    # Create header row for Raw Scores
    raw_headers = ['Rank', 'Antigen'] + [col.replace('_Score', '') for col in raw_score_cols] + ['Total']
    pdf.create_table_header(raw_headers, score_widths)
        
    # Add data rows for Raw Scores
    if top_n is not None and top_n < len(df):
        display_df = df.head(top_n)
    else:
        display_df = df
    
    # Sort by total score if available
    if 'total_score' in display_df.columns:
        display_df = display_df.sort_values('total_score', ascending=False)
        
    for idx, row in display_df.iterrows():
        data = [
            idx + 1,
            row['antigen']
        ]
        
        # Add raw scores
        for col in raw_score_cols:
            data.append(format_numeric(row.get(col, 0.0)))
        
        # Add total score
        data.append(format_numeric(row.get('total_score', 0.0)))
        
        pdf.create_table_row(data, score_widths)
    
    # Weighted Scores Table
    pdf.add_page()
    pdf.section_title('Weighted Category Scores')
    
    # Define column widths for weighted score table
    weighted_score_widths = [15, 80]  # Rank, Antigen
    weighted_col_width = 25
    weighted_score_widths.extend([weighted_col_width] * len(weighted_score_cols))
    
    # Create header row for Weighted Scores
    weighted_headers = ['Rank', 'Antigen'] + [col.replace('Weighted_', '') for col in weighted_score_cols]
    pdf.create_table_header(weighted_headers, weighted_score_widths)
    
    # Add data rows for Weighted Scores
    for idx, row in display_df.iterrows():
        data = [
            idx + 1,
            row['antigen']
        ]
        
        # Add weighted scores
        for col in weighted_score_cols:
            data.append(format_numeric(row.get(col, 0.0)))
        
        pdf.create_table_row(data, weighted_score_widths)
    
    # Data Visualizations
    create_visualizations(pdf, df, top_n)

    # Statistical Analysis
    pdf.add_page()
    pdf.section_title('Statistical Analysis by Category')
    stats_headers = ['Category', 'Mean', 'Std Dev', 'Min', 'Max']
    stats_widths = [80, 30, 30, 30, 30]
    pdf.create_table_header(stats_headers, stats_widths)
    
    # Statistics for raw scores
    for col in raw_score_cols:
        stats = df[col].describe()
        category = col.replace('_Score', '')
        data = [
            category,
            format_numeric(stats['mean'], 3),
            format_numeric(stats['std'], 3),
            format_numeric(stats['min'], 3),
            format_numeric(stats['max'], 3)
        ]
        pdf.create_table_row(data, stats_widths, is_alternating=True)
    
    # Statistics for total score
    if 'total_score' in df.columns:
        stats = df['total_score'].describe()
        data = [
            'Total Score',
            format_numeric(stats['mean'], 3),
            format_numeric(stats['std'], 3),
            format_numeric(stats['min'], 3),
            format_numeric(stats['max'], 3)
        ]
        pdf.create_table_row(data, stats_widths, is_alternating=True)

    # Save the PDF
    pdf.output(str(output_pdf))

def generate_report(antigen: str, mol_formula: str, output_dir: Path, timestamp: str):
    """
    Generate a PDF report for an protein
    """
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filenames using molecular formula
    base_filename = f"{mol_formula}_{timestamp}"
    pdf_path = output_dir / f"{base_filename}.pdf"
    csv_path = output_dir / f"{base_filename}.csv"

def generate_basic_ranking_report(input_csv: str, output_pdf: str):
    """Generate a simple PDF report from CSV when the enhanced report fails."""
    # Convert paths to Path objects
    input_csv = Path(input_csv)
    output_pdf = Path(output_pdf)
    
    # Read the data
    df = pd.read_csv(input_csv)
    
    # Create PDF
    pdf = ProteinReportPDF()
    
    # Add title page
    pdf.add_page()
    pdf.set_font('SpaceGrotesk', 'B', 16)
    pdf.cell(0, 10, 'Antibody Ranking Report', 0, 1, 'C')
    pdf.set_font('SpaceGrotesk', '', 10)
    pdf.cell(0, 10, f'Generated on: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', 0, 1, 'C')
    
    # Basic table of rankings
    pdf.add_page()
    pdf.section_title('Antibody Rankings')
    
    # Simplify the DataFrame for display
    display_cols = ['antigen']
    score_cols = [col for col in df.columns if col.endswith('_Score')]
    display_cols.extend(score_cols)
    
    if 'total_score' in df.columns:
        df = df.sort_values('total_score', ascending=False)
    
    # Create a simplified table
    col_widths = [15]  # Rank
    col_widths.append(80)  # Antigen
    col_widths.extend([25] * len(score_cols))  # Score columns
    
    headers = ['Rank', 'Antigen'] + [col.replace('_Score', '') for col in score_cols]
    pdf.create_table_header(headers, col_widths)
    
    for idx, row in df.iterrows():
        data = [idx + 1, row['antigen']]
        for col in score_cols:
            data.append(format_numeric(row[col]))
        pdf.create_table_row(data, col_widths)
    
    # Save the PDF
    pdf.output(str(output_pdf))

if __name__ == "__main__":
    input_csv = "protein_rankings_20250318_100659.csv"
    output_pdf = "enhanced_protein_ranking_report.pdf"
    generate_ranking_report(input_csv, output_pdf, top_n=10) 