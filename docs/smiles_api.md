# Moremi BioKit - SMILES Subpackage API (`moremi_biokit.smiles`)

This document provides detailed API documentation for the `moremi_biokit.smiles` subpackage, designed for developers and users working with SMILES-based small molecule data within the Moremi BioKit framework.

## 1. Overview

The `smiles` subpackage offers a comprehensive workflow for processing small molecules:

1. **Batch Processing:** Read and manage collections of SMILES strings from input files.
2. **Validation & Metric Calculation:** Validate SMILES strings and calculate a wide range of physicochemical and ADMET properties using various models and calculators.
3. **Ranking:** Score and rank molecules based on calculated metrics according to a defined configuration.
4. **Reporting:** Generate detailed reports (CSV, PDF) for individual molecules and overall rankings.

## 2. Core Components

These are the main classes orchestrating the workflow.

### 2.1. `BatchMoleculeProcessor`

**Purpose:** Manages the processing of multiple SMILES strings from an input file, coordinating validation, ranking, and report generation.

**Import:**

```python
from moremi_biokit.smiles import BatchMoleculeProcessor
# Or via the main package
from moremi_biokit import BatchMoleculeProcessor 
```

**Initialization:**

```python
def __init__(self, 
             input_file: str, 
             output_dir: Optional[str] = None, 
             generate_pdf: Optional[bool] = True, 
             generate_csv: Optional[bool] = True):
```

- `input_file` (str): Path to the input file containing SMILES strings (one per line).
- `output_dir` (Optional[str]): Path to the directory where all output files (logs, reports, rankings) will be saved. If `None`, a timestamped directory named `molecule_analysis_results/YYYYMMDD_HHMMSS` is created.
- `generate_pdf` (Optional[bool]): If `True` (default), generates PDF reports for each successfully processed molecule.
- `generate_csv` (Optional[bool]): If `True` (default), generates CSV reports for each successfully processed molecule.

**Key Methods:**

- `process_batch()`: Reads the input file, processes each SMILES string through the validator and ranker, generates reports based on the flags set during initialization, and saves ranking results. Handles logging and reports progress.
- `get_ranked_results_as_dict() -> List[Dict[str, Any]]`: Returns the final ranked results (after `process_batch()` has run) as a list of dictionaries, suitable for JSON serialization or other downstream processing. Returns an empty list if ranking failed or produced no results.

**Example Usage:**

```python
from moremi_biokit.smiles import BatchMoleculeProcessor

input_smiles = "path/to/your/molecules.smi"
output = "analysis_output"

# Process batch, generate only CSV reports
processor = BatchMoleculeProcessor(input_smiles, output, generate_pdf=False, generate_csv=True)
processor.process_batch()

# Get results as dictionary list
ranked_data = processor.get_ranked_results_as_dict()
if ranked_data:
    print(f"Processed and ranked {len(ranked_data)} molecules.")
```

**CLI Usage:**

The package provides a command-line script `process-smiles-batch` (defined in `pyproject.toml`).

```bash
process-smiles-batch <input_file> [-o <output_dir>] [--pdf] [--csv] [--no-pdf] [--no-csv]
```

- `<input_file>`: Path to the SMILES file.
- `-o`, `--output_dir`: Optional output directory.
- `--pdf`: Generate PDF reports.
- `--csv`: Generate CSV reports.
- `--no-pdf`: Explicitly disable PDF reports.
- `--no-csv`: Explicitly disable CSV reports.

*Default Behavior (CLI):* If no report flags (`--pdf`, `--csv`, `--no-pdf`, `--no-csv`) are provided, both PDF and CSV reports are generated.

### 2.2. `SmallMoleculeValidator` (from `small_molecule_validator_v3.py`)

**Purpose:** Validates individual SMILES strings and orchestrates the calculation of a wide range of molecular metrics using the `property_calculators`.

**Import:**

```python
from moremi_biokit.smiles import SmallMoleculeValidator
```

**Initialization:**

```python
def __init__(self, ranges: Optional[MetricRanges] = None):
```

- `ranges` (Optional[`MetricRanges`]): Optionally provide an instance of `MetricRanges` (though currently, ranges are primarily used internally for reference within the validator/ranker and not for pass/fail criteria).

**Key Methods:**

- `process_molecule(smiles: str) -> ProcessingResult`: Takes a single SMILES string, validates it, calculates all defined metrics using the various property calculators, and returns a `ProcessingResult` object.
- `process_molecules(input_source: Union[str, List[str]], output_dir: str) -> List[ProcessingResult]`: Processes multiple molecules from a file path or a list of SMILES strings. Returns a list of `ProcessingResult` objects.
- `get_successful_metrics(results: List[ProcessingResult]) -> List[MoleculeMetrics]`: A helper function to filter a list of `ProcessingResult` objects and return only the `MoleculeMetrics` from successful validations.

**Key Data Structures:**

- **`ProcessingResult`**: A dataclass returned by `process_molecule` and `process_molecules`.
  - `smiles` (str): The input SMILES string.
  - `metrics` (Optional[`MoleculeMetrics`]): Contains the calculated metrics if processing was successful.
  - `error` (Optional[str]): An error message if processing failed.
  - `success` (bool): Indicates whether processing was successful.
  - **`MoleculeMetrics`**: A dataclass holding the calculated metrics.
  - `smiles` (str): Canonical SMILES.
  - `molecular_formula` (str)
  - `molecular_weight` (float)
  - `physicochemical` (Dict[str, float])
  - `medicinal_chemistry` (Dict[str, float])
  - `lipophilicity` (Dict[str, float])
  - `druglikeness` (Dict[str, float])
  - `absorption` (Dict[str, float])
  - `distribution` (Dict[str, float])
  - `metabolism` (Dict[str, float])
  - `excretion` (Dict[str, float])
  - `toxicity` (Dict[str, float])
  - `warnings` (List[str])
  - `to_dict() -> Dict`: Method to convert the metrics object into a dictionary.
  - **`MetricCategory`**: An Enum defining the categories used for organizing metrics (e.g., `MetricCategory.PHYSICOCHEMICAL`, `MetricCategory.ABSORPTION`).

**Example Usage:**

```python
from moremi_biokit.smiles import SmallMoleculeValidator, MoleculeMetrics, ProcessingResult

validator = SmallMoleculeValidator()
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O" # Aspirin
result: ProcessingResult = validator.process_molecule(smiles)

if result.success and isinstance(result.metrics, MoleculeMetrics):
    print("Validation Successful!")
    metrics_dict = result.metrics.to_dict()
    print(f"Formula: {metrics_dict.get('molecular_formula')}")
    print(f"Absorption Metrics: {metrics_dict.get('metrics', {}).get('absorption')}")
else:
    print(f"Validation Failed: {result.error}")
```

### 2.3. `SmallMoleculeRankerV4` (from `small_molecule_ranker_v4.py`)

**Purpose:** Scores and ranks molecules based on the metrics collected by the `SmallMoleculeValidator`. Uses a `ScoringConfig` to define how properties are scored and weighted.

**Import:**

```python
from moremi_biokit.smiles import SmallMoleculeRankerV4, ScoringConfig
```

**Initialization:**

```python
def __init__(self, 
             config: Optional[ScoringConfig] = None, 
             generate_pdf: Optional[bool] = False, 
             generate_csv: Optional[bool] = True):
```

- `config` (Optional[`ScoringConfig`]): An instance of `ScoringConfig` defining weights and scoring functions. If `None`, uses the default configuration (`ScoringConfig` defined in the module, which uses equal weighting).
- `generate_pdf` (Optional[bool]): If `True` (default is `False`), instructs the ranker to generate individual PDF reports when `rank_molecules` is called (assuming `generate_enhanced_report` is used internally).
- `generate_csv` (Optional[bool]): If `True` (default), instructs the ranker to generate individual CSV reports when `rank_molecules` is called.

**Key Methods:**

- `set_output_directory(output_dir: str)`: Sets the directory where ranking reports and individual molecule reports (if generated by the ranker) will be saved.
- `rank_molecules(metrics_list: List[MoleculeMetrics]) -> pd.DataFrame`: Takes a list of `MoleculeMetrics` objects, calculates scores (overall, category, individual metrics) based on the configuration, ranks the molecules (higher score is better), optionally generates individual molecule reports, saves the final ranking report (CSV and PDF), and returns a pandas DataFrame containing the ranked molecules and their scores.
- `get_ranked_data_as_dict() -> List[Dict[str, Any]]`: Returns the ranked results (after `rank_molecules()` has run) as a list of dictionaries.
- `save_rankings(output_dir: Path, timestamp: str)`: Saves the internal DataFrame of ranked molecules to a timestamped CSV file in the specified output directory's `rankings` subfolder.
- `generate_report(output_dir: str, timestamp: str)`: Generates the final summary ranking report (CSV and PDF) in the specified output directory's `rankings` subfolder.

**Key Data Structures:**

- **`ScoringConfig`**: A dataclass (defined within `small_molecule_ranker_v4.py`) holding the ranking configuration.
- `category_weights` (Dict[`MetricCategory`, float]): Weights applied to each metric category score when calculating the overall score.
- `property_configs` (Dict[str, Dict]): Defines how individual metrics are scored. Each entry maps a metric name (e.g., 'tpsa', 'logp') to a dictionary containing an `optimal_range` (optional) and a `score_func` (lambda function taking the value and returning a score 0-1).

**Example Usage:**

```python
from moremi_biokit.smiles import SmallMoleculeRankerV4, MoleculeMetrics

# Assume metrics_list is a list of MoleculeMetrics objects obtained from a validator
metrics_list: List[MoleculeMetrics] = [...] 

# Use default ranking configuration, generate CSV reports only
ranker = SmallMoleculeRankerV4(generate_pdf=False, generate_csv=True)
ranker.set_output_directory("ranking_output")

ranked_df = ranker.rank_molecules(metrics_list)

print("Top 5 Ranked Molecules:")
print(ranked_df.head()['SMILES', 'Molecular Formula', 'Overall Score'])

# Get results as dictionary
ranked_dict = ranker.get_ranked_data_as_dict()
```

## 3. Sub-modules

These modules contain the underlying calculators and utilities.

### 3.1. `property_calculators`

**Purpose:** Provides individual classes for calculating specific groups of molecular properties.

**Import:**

```python
from moremi_biokit.smiles import property_calculators
# Or import a specific calculator
from moremi_biokit.smiles.property_calculators import ADMETPredictor
```

**Overview:**
This subpackage contains modules like:

- `admet_predictor.py` (`ADMETPredictor` class)
- `druglikeness.py` (`DruglikenessProperties` class)
- `lipophilicity.py` (`LipophilicityProperties` class)
- `medicinal_chemistry.py` (`MedicinalChemistryProperties` class)
- `pharmacokinetics.py` (`PharmacokineticsProperties` class)
- `physicochemical.py` (`PhysiochemicalProperties` class)
- `solubility.py` (`SolubilityProperties` class)

Each class is typically initialized with a SMILES string and provides methods (e.g., `calculate_admet_properties()`, `calculate_druglikeness()`) to compute the relevant metrics.

**Note:** While these calculators can be used directly, their primary role is to be invoked by the `SmallMoleculeValidator`.

**Example (Direct Usage):**

```python
from moremi_biokit.smiles.property_calculators import SolubilityProperties

solub_calc = SolubilityProperties("CCO") # Ethanol
results = solub_calc.calculate_solubility()
print(f"Ethanol ESOL LogS: {results.get('ESOL', {}).get('log_s')}")
```

### 3.2. `reports`

**Purpose:** Handles the generation of output reports.

**Import:**

```python
from moremi_biokit.smiles import reports
# Or import specific functions/classes
from moremi_biokit.smiles.reports import generate_enhanced_report, EnhancedMoleculeReport
```

**Key Components:**

- `enhanced_report_generator.py`:
- `EnhancedMoleculeReport` class: Inherits from `fpdf.FPDF` to create detailed PDF reports with custom formatting, including structure images and radar charts.
- `generate_enhanced_report(molecule_data: Dict, output_dir: str, generate_pdf: bool, generate_csv: bool)`: Function to generate both PDF (optional) and CSV (optional) reports for a single molecule's data. Requires a dictionary (`molecule_data`) containing all metrics and scores.
- `generate_csv_report(data: Dict, output_path: str)`: Generates a CSV from a (potentially nested) dictionary of molecule data.
- `report_generator.py`: (Assumed to exist based on imports in ranker)
- `generate_ranking_report(csv_file: str, pdf_file: str, top_n: int)`: Generates a summary PDF report from a ranking CSV file.

**Note:** These are primarily used internally by the `BatchMoleculeProcessor` and `SmallMoleculeRankerV4` but can be used independently if you have the required data structures.

### 3.3. `utils`

**Purpose:** Contains utility functions supporting the SMILES workflow.

**Import:**

```python
from moremi_biokit.smiles import utils
# Or import specific components
from moremi_biokit.smiles.utils import MoleculeConverter
```

**Key Components:**

- `convert_smiles.py`:
- `MoleculeConverter` class: Converts SMILES strings to 3D structures (PDB format) and generates interactive 3D visualizations (HTML using py3Dmol). Has methods to process single or multiple SMILES.
- Provides a `main()` function for CLI usage via the `convert-smiles` script.

**CLI Usage (`convert-smiles`):**

```bash
convert-smiles --smiles <input_smiles_file> --output <output_directory>
```

- `--smiles`: Input file containing SMILES strings.
- `--output`: Directory to save generated PDB and HTML files.

## 4. Key Data Classes Summary

- **`MoleculeMetrics`**: (`smiles.small_molecule_validator_v3`) Holds all calculated properties for a molecule, categorized. Includes a `to_dict()` method.
- **`ProcessingResult`**: (`smiles.small_molecule_validator_v3`) Wraps the result of `process_molecule`, indicating success/failure and holding either `MoleculeMetrics` or an error message.
- **`ScoringConfig`**: (`smiles.small_molecule_ranker_v4`) Defines weights and functions for scoring metrics during ranking.
- **`MetricCategory`**: (`smiles.small_molecule_validator_v3`) Enum for metric categories.
- **`MetricRanges`**: (`smiles.small_molecule_validator_v3`) Holds reference ranges (mostly informational).

This documentation provides a guide to the main components and usage patterns of the `moremi_biokit.smiles` subpackage. For implementation details, refer to the specific module and class docstrings in the source code.
