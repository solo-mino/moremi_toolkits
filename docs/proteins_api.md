# Moremi BioKit - Proteins Subpackage API (`moremi_biokit.proteins`)

This document provides detailed API documentation for the `moremi_biokit.proteins` subpackage, designed for developers and users working with protein sequence analysis, antibody validation, and ranking within the Moremi BioKit framework.

## 1. Overview

The `proteins` subpackage offers a comprehensive workflow for processing protein sequences, particularly antibodies:

1. **Validation & Metric Calculation:** Validate protein sequences and calculate a wide range of physicochemical, structural, binding affinity, and functional properties using various analysis tools.
2. **Ranking:** Score and rank proteins based on calculated metrics according to a defined configuration with weighted scoring.
3. **Reporting:** Generate detailed reports (CSV, PDF) for individual proteins and overall rankings.
4. **Real-time Backup:** Provides data resilience through incremental backup during long-running processes.

## 2. Core Components

These are the main classes orchestrating the workflow.

### 2.1. `ProteinRanker`

**Purpose:** Manages the ranking of proteins based on calculated metrics, providing comprehensive scoring and report generation capabilities.

**Import:**

```python
from moremi_biokit.proteins import ProteinRanker, ScoringConfig
# Or via the main package
from moremi_biokit import ProteinRanker 
```

**Initialization:**

```python
def __init__(self, 
             config: Optional[ScoringConfig] = None, 
             generate_pdf: bool = False, 
             generate_csv: bool = True):
```

- `config` (Optional[ScoringConfig]): Custom scoring configuration defining weights and scoring functions. If `None`, uses default configuration with predefined category weights.
- `generate_pdf` (Optional[bool]): If `True`, generates PDF reports for individual proteins during ranking. Defaults to `False`.
- `generate_csv` (Optional[bool]): If `True` (default), generates CSV reports for individual proteins during ranking.

**Key Methods:**

- `set_output_directory(output_dir: str)`: Sets the directory where ranking reports and individual protein reports will be saved.
- `rank_proteins(metrics_input: Union[List[ProteinMetrics], ProteinMetrics, str, Path]) -> pd.DataFrame`: Main ranking method that accepts protein metrics and returns a ranked DataFrame. Can accept a list of `ProteinMetrics` objects, a single object, or a path to a JSON file containing serialized metrics.
- `calculate_overall_score(metrics: ProteinMetrics) -> Dict[str, Any]`: Calculates comprehensive scores for a single protein, returning overall score, category scores, and detailed metric scores.
- `save_rankings(output_dir: Path, timestamp: str) -> Optional[Path]`: Saves ranking results to CSV with proper column ordering.
- `generate_ranking_report(csv_file: str, output_dir: str, timestamp: str) -> Optional[Path]`: Generates comprehensive PDF ranking reports.

**Key Data Structures:**

- **`ScoringConfig`**: A dataclass holding the ranking configuration.
  - `category_weights` (Dict[MetricCategory, float]): Weights applied to each metric category for overall score calculation.
  - `property_configs` (Dict[str, Dict]): Defines scoring functions for individual metrics with optimal ranges and custom scoring functions.

**Example Usage:**

```python
from moremi_biokit.proteins import ProteinRanker, ScoringConfig, MetricCategory

# Create custom configuration
custom_config = ScoringConfig()
custom_config.category_weights[MetricCategory.BINDING_AFFINITY] = 0.4
custom_config.category_weights[MetricCategory.STABILITY] = 0.3

# Initialize ranker
ranker = ProteinRanker(config=custom_config, generate_pdf=True, generate_csv=True)
ranker.set_output_directory("ranking_outputs")

# Rank proteins (assuming metrics_list is available)
ranked_df = ranker.rank_proteins(metrics_list)

print("Top 5 Ranked Proteins:")
print(ranked_df[['sequence', 'antigen', 'total_score']].head())
```

### 2.2. `ProteinValidatorV2` (from `protein_validator_v2.py`)

**Purpose:** Validates protein sequences and orchestrates the calculation of comprehensive protein metrics across multiple analysis categories.

**Import:**

```python
from moremi_biokit.proteins import ProteinValidatorV2, MetricCategory
```

**Initialization:**

```python
def __init__(self, 
             pdb_files_path: str = "./pdb_files", 
             metrics_to_run: Optional[List[MetricCategory]] = None):
```

- `pdb_files_path` (str): Directory path where PDB files will be stored during structure analysis.
- `metrics_to_run` (Optional[List[MetricCategory]]): Specific metrics to calculate. If `None`, runs all available metrics.

**Key Methods:**

- `validate_protein_list(input_source: Union[str, List[str]], output_csv_path: str, realtime_csv_backup_path: Optional[str] = None, realtime_json_backup_path: Optional[str] = None) -> List[ProteinMetrics]`: Processes multiple proteins with real-time backup support.
- `process_protein(sequence: str) -> ProteinMetrics`: Validates and analyzes a single protein sequence.
- `set_antigen_context(target_antigen_sequence: Optional[str], target_antigen_pdb_file_path: Optional[str], target_antigen_pdb_chain_id: Optional[str], target_pdb_id: Optional[str] = None, antigen_pdb_download_dir: str = "./antigen_pdbs") -> bool`: Sets antigen context for binding affinity calculations.

**Key Data Structures:**

- **`ProteinMetrics`**: A dataclass holding all calculated protein metrics.
  - `sequence` (str): Protein sequence
  - `antigen` (Optional[str]): Antigen sequence if provided
  - `molecular_formula` (str): Molecular formula
  - `molecular_weight` (float): Molecular weight
  - `protparam` (Dict): ProtParam analysis results
  - `binding_affinity` (Dict): Binding affinity metrics
  - `structure` (Dict): Structural quality metrics
  - `glycosylation` (Dict): Glycosylation site analysis
  - `aggregation` (Dict): Aggregation propensity analysis
  - `stability` (Dict): Thermal stability predictions
  - `immunogenicity` (Dict): Immunogenicity predictions
  - `conservancy` (Dict): Conservancy analysis
  - `epitope` (Dict): Epitope mapping results
  - `developability` (Dict): Developability assessment
  - `blast` (Dict): BLAST search results
  - `warnings` (List[str]): Validation warnings
  - `to_dict() -> Dict`: Method to serialize metrics
  - `from_dict(data: Dict) -> ProteinMetrics`: Class method to deserialize metrics

- **`MetricCategory`**: An Enum defining metric categories:
  - `PROTPARAM`: Basic physicochemical properties
  - `BINDING_AFFINITY`: Antigen-antibody binding analysis
  - `STRUCTURE`: Structural quality assessment
  - `GLYCOSYLATION`: Glycosylation site analysis
  - `AGGREGATION`: Aggregation propensity
  - `STABILITY`: Thermal stability
  - `IMMUNOGENICITY`: Immunogenicity prediction
  - `CONSERVANCY`: Sequence conservancy
  - `EPITOPE`: Epitope mapping
  - `DEVELOPABILITY`: Developability assessment
  - `BLAST`: Database searches

**Example Usage:**

```python
from moremi_biokit.proteins import ProteinValidatorV2, MetricCategory

# Initialize validator for specific metrics
validator = ProteinValidatorV2(
    pdb_files_path="./analysis_pdbs",
    metrics_to_run=[
        MetricCategory.PROTPARAM,
        MetricCategory.BINDING_AFFINITY,
        MetricCategory.STABILITY
    ]
)

# Set antigen context for binding analysis
validator.set_antigen_context(
    target_antigen_pdb_chain_id="1A2Y_A",
    antigen_pdb_download_dir="./antigens"
)

# Process protein sequences
sequences = ["EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKG...", ...]
metrics_list = validator.validate_protein_list(
    input_source=sequences,
    output_csv_path="validation_results.csv"
)

print(f"Successfully validated {len(metrics_list)} proteins")
```

### 2.3. Main Processing Function

**`rank_proteins_from_metrics`**

**Purpose:** Main entry point function that orchestrates the complete protein analysis workflow from sequence input to final rankings.

**Import:**

```python
from moremi_biokit.proteins import rank_proteins_from_metrics
```

**Function Signature:**

```python
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
```

**Parameters:**

- `protein_sequences_input`: Input sequences as file path, list of sequences, or Path object
- `output_dir`: Directory for all output files and reports
- `config`: Optional custom scoring configuration
- `generate_pdf`: Whether to generate PDF reports for individual proteins
- `generate_csv`: Whether to generate CSV reports for individual proteins
- `metrics_to_run`: Specific metrics to calculate (if None, uses all metrics)
- `target_antigen_sequence`: Antigen sequence for binding analysis
- `target_antigen_pdb_file_path`: Path to antigen PDB file
- `target_antigen_pdb_chain_id`: PDB chain ID for antigen download
- `antigen_pdb_download_path`: Directory for downloaded antigen structures

**Returns:**

- `pd.DataFrame`: Ranked proteins with comprehensive scoring data

**Example Usage:**

```python
from moremi_biokit.proteins import rank_proteins_from_metrics, MetricCategory

# Complete analysis workflow
ranked_df = rank_proteins_from_metrics(
    protein_sequences_input="antibody_sequences.fasta",
    output_dir="protein_analysis_results",
    metrics_to_run=[
        MetricCategory.PROTPARAM,
        MetricCategory.BINDING_AFFINITY,
        MetricCategory.STABILITY,
        MetricCategory.IMMUNOGENICITY
    ],
    target_antigen_pdb_chain_id="6M0J_E",  # SARS-CoV-2 spike protein
    generate_pdf=True,
    generate_csv=True
)

print(f"Analysis complete. Ranked {len(ranked_df)} proteins.")
print("Top 3 antibodies:")
print(ranked_df[['sequence', 'antigen', 'total_score']].head(3))
```

## 3. Analysis Tools

The `analysis_tools` submodule contains specialized analysis components.

### 3.1. Key Analysis Modules

**Import:**

```python
from moremi_biokit.proteins.analysis_tools import (
    ProteinComparisonTool,  # from developability.py
    # Other analysis tools as needed
)
```

**Overview:**

- `developability.py`: `ProteinComparisonTool` for developability assessment using SAbDab database comparisons
- Other specialized analysis tools for various protein properties

**Example (Direct Usage):**

```python
from moremi_biokit.proteins.analysis_tools.developability import ProteinComparisonTool

# Initialize comparison tool
comp_tool = ProteinComparisonTool()
result = comp_tool.compare_sequence("EVQLVESGGGLVQPGGSLRLSCAAS...")
print(f"Developability score: {result.get('developability_score', 'N/A')}")
```

## 4. Reports Module

**Purpose:** Handles generation of comprehensive output reports.

**Import:**

```python
from moremi_biokit.proteins import reports
# Or import specific functions/classes
from moremi_biokit.proteins.reports import generate_enhanced_report, generate_ranking_report
```

**Key Components:**

- `protein_enhanced_report_generator.py`:
  - `EnhancedAntibodyReport` class: Creates detailed PDF reports with structure analysis, radar charts, and comprehensive property tables
  - `generate_enhanced_report(protein_data: Dict, output_dir: str, generate_pdf: bool, generate_csv: bool)`: Main function for individual protein reports
  
- `protein_report_generator.py`:
  - `ProteinReportPDF` class: Generates summary ranking reports
  - `generate_ranking_report(input_csv: str, output_pdf: str, top_n: int)`: Creates comprehensive ranking summary PDFs

**Example:**

```python
from moremi_biokit.proteins.reports import generate_enhanced_report

# Generate detailed report for a single protein
protein_data = {
    "sequence": "EVQLVESGGGLVQ...",
    "molecular_formula": "C627H1254N627O627",
    "protparam": {...},
    "binding_affinity": {...},
    # ... other metrics
}

report_paths = generate_enhanced_report(
    protein_data=protein_data,
    output_dir="protein_reports",
    generate_pdf=True,
    generate_csv=True
)
```

## 5. Real-time Backup and Data Resilience

The proteins subpackage provides robust data backup mechanisms for long-running analyses:

### 5.1. Backup Features

- **Real-time CSV Backup:** All validation attempts (successful and failed) are incrementally saved
- **Real-time JSON Backup:** Successfully processed protein metrics are saved for recovery
- **Resume Capability:** Rankings can be resumed from backup files if processes are interrupted

### 5.2. Backup File Locations

When using `rank_proteins_from_metrics`, backup files are automatically created:

```plaintext
output_dir/
├── ranking_run_YYYYMMDD_HHMMSS/
│   ├── realtime_validation_attempts_backup.csv
│   ├── realtime_successful_metrics_for_ranking.json
│   ├── validation_attempts_summary.csv
│   ├── protein_reports/
│   ├── rankings/
│   └── ranking_process.log
```

### 5.3. Recovery Usage

```python
from moremi_biokit.proteins import ProteinRanker

# Resume ranking from backup JSON
ranker = ProteinRanker()
ranker.set_output_directory("recovery_output")

# Load from backup and resume ranking
backup_json = "previous_run/realtime_successful_metrics_for_ranking.json"
ranked_df = ranker.rank_proteins(backup_json)
```

## 6. Configuration and Customization

### 6.1. Custom Scoring Configuration

```python
from moremi_biokit.proteins import ScoringConfig, MetricCategory

# Create custom configuration
config = ScoringConfig()

# Adjust category weights
config.category_weights = {
    MetricCategory.BINDING_AFFINITY: 0.35,
    MetricCategory.STABILITY: 0.25,
    MetricCategory.IMMUNOGENICITY: 0.20,
    MetricCategory.PROTPARAM: 0.10,
    MetricCategory.AGGREGATION: 0.10,
    # Other categories: 0.0 (disabled)
}

# Customize individual metric scoring
config.property_configs['dissociation_constant']['optimal_range'] = (1e-11, 1e-6)
config.property_configs['melting_temperature_celsius']['optimal_range'] = (70, 95)
```

### 6.2. Dynamic Weight Assignment

When `metrics_to_run` is specified, weights are automatically distributed equally among selected metrics:

```python
# Equal weighting among selected metrics
ranked_df = rank_proteins_from_metrics(
    protein_sequences_input="sequences.txt",
    output_dir="results",
    metrics_to_run=[
        MetricCategory.BINDING_AFFINITY,  # Will get 1/3 weight
        MetricCategory.STABILITY,         # Will get 1/3 weight  
        MetricCategory.IMMUNOGENICITY    # Will get 1/3 weight
    ]
)
```

## 7. Key Data Classes Summary

- **`ProteinMetrics`**: (`proteins.protein_validator_v2`) Comprehensive protein analysis results with serialization support
- **`ScoringConfig`**: (`proteins.protein_ranker`) Ranking configuration with category weights and scoring functions
- **`MetricCategory`**: (`proteins.protein_validator_v2`) Enum defining all available analysis categories
- **`ProteinRanker`**: Main ranking engine with scoring and report generation
- **`ProteinValidatorV2`**: Comprehensive protein validation and metric calculation

## 8. Integration with Other Subpackages

The proteins subpackage integrates seamlessly with other Moremi BioKit components:

```python
# Combined analysis workflow
from moremi_biokit import BatchMoleculeProcessor, rank_proteins_from_metrics

# Process small molecules
molecule_processor = BatchMoleculeProcessor("compounds.smi", "molecule_results")
molecule_processor.process_batch()

# Process proteins
protein_results = rank_proteins_from_metrics(
    "antibodies.fasta", 
    "protein_results",
    target_antigen_pdb_chain_id="6M0J_E"
)

print("Integrated analysis complete!")
```

This documentation provides a comprehensive guide to the `moremi_biokit.proteins` subpackage. For implementation details and advanced usage, refer to the specific module and class docstrings in the source code.
