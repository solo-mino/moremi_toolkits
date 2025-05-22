# Moremi Biokit (`moremi-biokit`) Documentation

DRAFT V 0.1.0

## 1. Introduction

**`moremi-biokit`** is a foundational Python toolkit designed to provide robust processing, validation, analysis, and reporting capabilities for common biological entities. It is a core component used by the Moremi Bio Autonomous Agent and other internal projects requiring these functionalities.

Currently, `moremi-biokit` includes two main sub-packages:

- **`moremi_biokit.smiles`**: For handling small molecules represented by SMILES strings.
- **`moremi_biokit.proteins`**: For processing and analyzing protein data.

This document provides an overview of the package structure, installation, usage, and development guidelines.

**Important Note:** `moremi-biokit` is a library providing tools and functions. It is *not* an autonomous system by itself. The Moremi Bio agent (a separate project) utilizes this biokit to perform its tasks.

## 2. Package Structure

The `moremi_biokit` package is meticulously organized to provide a clear and maintainable codebase. It contains distinct sub-packages for `smiles` and `proteins` processing, each with its own internal structure reflecting its specific functionalities.

```plaintext
moremi_toolkits/               # Monorepo Root
├── .git/
├── .gitignore                 # Monorepo .gitignore (can include common Python patterns)
├── LICENSE                    # Monorepo LICENSE file (for the whole monorepo)
├── README.md                  # Monorepo root README (overview of the entire 'moremi_toolkits' project)
├── pyproject.toml             # Optional: Monorepo root pyproject.toml (for dev tools like pre-commit, linters applied to whole repo)
├── .pre-commit-config.yaml    # Optional: For pre-commit hooks across the monorepo
│
└── components/                # Directory housing all installable components/packages
    └── moremi_biokit/         # === Root directory for the 'moremi-biokit' PACKAGE ===
        ├── moremi_biokit/     # --- The actual Python package source code ---
        │   ├── __init__.py    # Main __init__ for the 'moremi_biokit' package
        │   │
        │   ├── smiles/        # Sub-package for Small Molecule (SMILES) related tools
        │   │   ├── __init__.py
        │   │   ├── batch_processor.py  
        │   │   ├── validator.py
        │   │   ├── ranker.py
        │   │   ├── property_calculators/
        │   │   │   ├── __init__.py
        │   │   │   ├── physicochemical.py
        │   │   │   ├── lipophilicity.py
        │   │   │   ├── solubility.py
        │   │   │   ├── druglikeness.py
        │   │   │   ├── medicinal_chemistry.py
        │   │   │   ├── pharmacokinetics.py
        │   │   │   └── admet_predictor.py
        │   │   └── reports/
        │   │       ├── __init__.py
        │   │       ├── basic_generator.py
        │   │       └── enhanced_generator.py
        │   │
        │   └── proteins/    # Sub-package for Antibody related tools
        │       ├── __init__.py
        │       ├── batch_processor.py
        │       ├── validator.py
        │       ├── ranker.py
        │       ├── analysis_tools/
        │       │   ├── __init__.py
        │       │   ├── blast_analyzer.py
        │       │   ├── protparam_analyzer.py
        │       │   ├── immunogenicity_predictor.py
        │       │   ├── stability_predictor.py
        │       │   ├── aggregation_predictor.py
        │       │   ├── glycosylation_analyzer.py
        │       │   ├── structure_predictor.py
        │       │   ├── affinity_predictor.py
        │       │   ├── epitope_predictor.py
        │       │   ├── conservancy_analyzer.py
        │       │   └── developability_assessor.py
        │       ├── reports/
        │       │   ├── __init__.py
        │       │   ├── basic_generator.py
        │       │   └── enhanced_generator.py
        │       └── utils/
        │           ├── __init__.py
        │           └── antigens.py
        │
        ├── pyproject.toml     # === BUILD CONFIGURATION for the 'moremi-biokit' package ===
        │                      # Defines dependencies, version, package name ("moremi-biokit"), etc.
        │
        ├── README.md          # === README specific to the 'moremi-biokit' package ===
        │                      # Details on how to use and develop this specific package.
        │
        └── tests/             # === UNIT TESTS for the 'moremi-biokit' package ===
            ├── __init__.py    # So 'tests' can be a package itself
            ├── conftest.py    # Pytest fixtures shared among tests for 'moremi-biokit'
            │
            ├── smiles/        # Tests specifically for the 'moremi_biokit.smiles' sub-package
            │   ├── __init__.py
            │   ├── test_batch_processor.py
            │   ├── test_validator.py
            │   ├── test_ranker.py
            │   └── property_calculators/
            │       ├── __init__.py
            │       └── test_physicochemical.py
            │       └── ... (other property calculator tests)
            │
            └── proteins/    # Tests specifically for the 'moremi_biokit.proteins' sub-package
                ├── __init__.py
                ├── test_batch_processor.py
                ├── test_validator.py
                └── analysis_tools/
                    ├── __init__.py
                    └── test_immunogenicity_predictor.py
                    └── ... (other analysis tool tests)
```

### **Explanation of Key Structural Choices:**

- **`moremi_biokit` (Root Package):** The single installable package. Its `__init__.py` could potentially re-export the main classes from `smiles` and `proteins` for easier top-level imports if desired (e.g., `from moremi_biokit import SMILESValidator`).
- **`smiles/` and `proteins/` (Sub-packages):** Clear separation of concerns. Each acts as a self-contained toolkit for its respective entity type.
  - **`__init__.py` in sub-packages:** These are crucial. They should import and expose the main public classes and functions from the modules within them. For example, `moremi_biokit/smiles/__init__.py` might contain:

        ```python
        from .batch_processor import SMILESBatchProcessor
        from .validator import SMILESValidator # Assuming this is the new name
        from .ranker import SMILESRanker     # Assuming this is the new name
        # ... and potentially key data classes like MoleculeMetrics
        ```

        This allows users to do `from moremi_biokit.smiles import SMILESValidator` instead of delving deeper.
  - **`batch_processor.py`, `validator.py`, `ranker.py`:** These core workflow components are present in both `smiles` and `proteins`, tailored to their specific data types. The class names should be distinct (e.g., `SMILESValidator`, `ProteinValidator`).
  - **`property_calculators/` (in `smiles`):** Groups all modules responsible for calculating specific properties of small molecules.
  - **`analysis_tools/` (in `proteins`):** Analogous to `property_calculators`, this groups modules for various protein-specific analyses and predictions.
  - **`reports/` (in both):** Centralizes report generation logic for each sub-package.
  - **`utils/` (in `proteins`):** For utilities or data files specifically related to protein processing (like the `antigens.py` example). A similar `utils/` could exist in `smiles/` if needed.

### **Impact on Usage (from Section 4):**

The usage examples would reflect this more granular structure if accessing deeper modules directly:

```python
# Accessing a specific property calculator within smiles
from moremi_biokit.smiles.property_calculators import physicochemical
# mw = physicochemical.calculate_molecular_weight("CCO") # If function exposed

# Or if Validator exposes metrics in a structured way:
# metrics = smiles_validator.process_molecule("CCO").metrics
# mw = metrics.physicochemical.molecular_weight # (Depends on MoleculeMetrics structure)

# Accessing an protein analysis tool
from moremi_biokit.proteins.analysis_tools import immunogenicity_predictor
# immunogenicity_score = immunogenicity_predictor.predict_t_cell_epitopes(sequence_vh, sequence_vl)
```

However, the primary interaction for users of the `biokit` would ideally be through the main `Validator`, `BatchProcessor`, and `Ranker` classes, which internally orchestrate calls to these more granular modules. The `__init__.py` files play a key role in defining this public API of each sub-package.

This detailed structure should give your development team a clear blueprint for organizing the code from your existing pipelines into the new `moremi-biokit` package.

## 3. Installation

`moremi-biokit` is designed to be installed as a Python package dependency.

### 3.1. From Internal Git Repository (Recommended)

Add the following to your project's `requirements.txt` (or equivalent for Poetry, PDM, etc.), replacing `<your_internal_git_repo_url>` with the actual URL of the `moremi_toolkits` monorepo and `<tag_or_branch>` with the desired version/branch:

```bash
moremi-biokit @ git@github.com:solo-mino/moremi_toolkits.git
```

```bash
moremi-biokit @ git@github.com:solo-mino/moremi_toolkits.git@<tag_or_branch>#subdirectory=components/moremi_biokit
```

**Examples:**

- **Using a specific tag (e.g., `v0.1.0` - recommended for stability):**

    ```bash
    moremi-biokit @ git@github.com:solo-mino/moremi_toolkits.git@v0.1.0#subdirectory=components/moremi_biokit
    ```

- **Using the `main` branch (for development builds):**

    ```bash
    moremi-biokit @ git@github.com:solo-mino/moremi_toolkits.git@main#subdirectory=components/moremi_biokit
    ```

Then, install using pip:

```bash
pip install -r requirements.txt
```

### 3.2. For Local Development of `moremi-biokit`

If you are actively developing `moremi-biokit` itself:

1. Ensure you have cloned the `moremi_toolkits` monorepo.
2. Navigate to the root of the monorepo.
3. Install in "editable" mode into your active virtual environment:

    ```bash
    pip install -e ./components/moremi_biokit/
    ```

    Changes made to the source code in `./components/moremi_biokit/moremi_biokit/` will be immediately reflected.

## 4. Usage Examples

### 4.1. Using `moremi_biokit.smiles`

Refer to your existing `smiles_pipeline.md` for detailed component functionalities. The import paths will now be prefixed with `moremi_biokit.smiles`.

```python
# Assuming class names from your smiles_pipeline.md, adapted for the new structure
from moremi_biokit.smiles.validator import SmallMoleculeValidator # Or your actual class name
from moremi_biokit.smiles.batch_processor import BatchMoleculeProcessor
from moremi_biokit.smiles.property_calculators import physicochemical, lipophilicity

# --- Example: Validating a single molecule ---
smiles_validator = SmallMoleculeValidator()
molecule_smiles = "CCO" # Ethanol
processing_result = smiles_validator.process_molecule(molecule_smiles)

if processing_result.success:
    print(f"Successfully processed SMILES: {molecule_smiles}")
    # Access metrics, e.g.:
    # print(f"Molecular Weight: {processing_result.metrics.physicochemical_props.get('MolecularWeight')}")
    # print(f"Consensus LogP: {processing_result.metrics.lipophilicity_props.get('ConsensusLogP')}")
else:
    print(f"Failed to process {molecule_smiles}: {processing_result.error_message}")


# --- Example: Conceptual Batch Processing ---
# Ensure your BatchMoleculeProcessor is adapted to be initialized and called appropriately
# smiles_file = "path/to/input_smiles.smi"
# output_directory = "path/to/output_reports"
#
# batch_proc = BatchMoleculeProcessor(input_file_path=smiles_file, output_dir=output_directory)
# batch_proc.initialize_components(...) # If needed
# batch_results = batch_proc.process_batch()
#
# print(f"Successfully processed: {len(batch_results.get('successful_metrics', []))}")
# print(f"Failed: {len(batch_results.get('failed_molecules', []))}")

# --- Example: Accessing a specific property calculator ---
# (This depends on how you expose functions/classes in property_calculators/__init__.py)
# mol_weight = physicochemical.calculate_molecular_weight("CN(C)C(=O)c1ccccc1N")
# print(f"RDKit Mol Weight: {mol_weight}")
```

### 4.2. Using `moremi_biokit.proteins` (Hypothetical)

This section should be filled out as the protein tools are developed.

```python
from moremi_biokit.proteins.sequence_analyzer import AntibodySequenceAnalyzer # Example class
# from moremi_biokit.proteins.processor import AntibodyDataProcessor

# analyzer = AntibodySequenceAnalyzer(sequence="QVQLQESGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYWRSDTKDYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR...")
# print(f"Calculated pI: {analyzer.calculate_isoelectric_point()}")
# print(f"Heavy Chain CDR3: {analyzer.get_cdr3_heavy()}")
```

## 5. Development and Contribution

This section is for developers working directly on `moremi-biokit`.

### 5.1. Project Setup (Recap for `moremi-biokit` context)

1. Clone the `moremi_toolkits` monorepo.
2. Create and activate a Python virtual environment.
3. Install `moremi-biokit` in editable mode:
    `pip install -e ./components/moremi_biokit/`
4. Install development dependencies (for tests, linters, etc.):
    `pip install -e "./components/moremi_biokit/[dev]"` (assuming a `[dev]` extra is defined in `pyproject.toml`)

### 5.2. `pyproject.toml`

The main build configuration and dependency list for `moremi-biokit` is in `components/moremi_biokit/pyproject.toml`.
Ensure it includes:

- `name = "moremi-biokit"`
- Correct `version` (follow SemVer)
- List of `dependencies` (e.g., `rdkit-pypi`, `pandas`, `biopython`)
- Optional `[project.optional-dependencies]` for `dev` tools (e.g., `pytest`, `black`, `flake8`).

Example snippet for `components/moremi_biokit/pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "moremi-biokit"
version = "0.1.0" # Increment as features/fixes are added
description = "Core biological processing toolkit for the Moremi ecosystem."
readme = "README.md"
requires-python = ">=3.8"
license = { file = "../../../LICENSE" } # Points to the monorepo root LICENSE
authors = [
  { name = "Moremi Development Team", email = "dev-contact@yourcompany.com" },
]
# Add classifiers if you plan to ever publish to an internal PyPI
# classifiers = [
#     "Development Status :: 3 - Alpha",
#     "Intended Audience :: Developers",
#     "Programming Language :: Python :: 3",
#     # ...
# ]

dependencies = [
    "rdkit-pypi >= 2022.09.5", # For SMILES processing
    "pandas >= 1.3.0",         # General data handling
    "biopython >= 1.79",       # For protein sequence analysis (example)
    # ... other core dependencies for SMILES or proteins
]

[project.optional-dependencies]
dev = [
    "pytest >= 7.0",
    "pytest-cov",
    "black",
    "flake8",
    "mypy", # Optional: for static type checking
    "pre-commit", # For git hooks
]

# If you want to expose any command-line scripts directly from moremi-biokit
# [project.scripts]
# analyze-smiles-file = "moremi_biokit.smiles.cli_tool:main"
```

### 5.3. Testing

- Unit tests are located in `components/moremi_biokit/tests/`.
- Organize tests by module (e.g., `test_smiles_validator.py`, `test_protein_sequencer.py`).
- Aim for high test coverage.
- Use `pytest` (or your chosen test runner).

    ```bash
    # From the monorepo root:
    pytest ./components/moremi_biokit/tests/
    # With coverage (if pytest-cov is installed):
    pytest --cov=components/moremi_biokit/moremi_biokit ./components/moremi_biokit/tests/
    ```

### 5.4. Code Style and Quality

- **Formatter:** Use `black` for consistent code formatting.
- **Linter:** Use `flake8` for style guide enforcement.
- **Type Hinting:** Use Python type hints and consider `mypy` for static analysis.
- **Pre-commit Hooks:** Set up `pre-commit` (via a `.pre-commit-config.yaml` in the monorepo root) to automatically run linters/formatters before commits. This helps maintain code quality across all contributions.

### 5.5. Versioning and Releasing

- Follow Semantic Versioning (Major.Minor.Patch - e.g., `0.1.0`, `1.0.0`).
- Increment the `version` string in `components/moremi_biokit/pyproject.toml`.
- After changes are merged to the main branch and a new version is decided:
    1. Update the version in `pyproject.toml`.
    2. Commit the version change (e.g., `chore: Bump version to v0.2.0`).
    3. Create a Git tag for the release (e.g., `git tag v0.2.0`).
    4. Push the commit and the tag to the remote repository (e.g., `git push && git push --tags`).
- Consuming projects can then pin to this new tag.

## 6. Future Development

(Outline potential future additions or areas of improvement for `moremi-biokit`)

- Integration of new property predictors.
- Support for other biological entity types.
- Advanced reporting features.
- Performance optimizations.

## 7. Troubleshooting

(Add common issues and solutions as they arise)

- Dependency conflicts (especially RDKit with other packages).
- Environment setup issues.

---

**Key considerations for your internal team:**

1. **Clarity on `Moremi Bio` vs. `moremi-biokit`:** Reinforce this distinction in team communications, onboarding, and any architectural diagrams.
2. **Internal Git Repo URL:** Make sure all placeholders like `<your_internal_git_repo_url>` are updated with the actual URL.
3. **LICENSE file:** Create an appropriate internal `LICENSE` file in the root of the `moremi_toolkits` monorepo.
4. **`pyproject.toml` details:**
    - Fill in actual author names/emails.
    - Thoroughly list all dependencies for `moremi-biokit`.
    - Define the `dev` dependencies your team uses.
5. **Contribution Guidelines in root README:** Flesh this out according to your team's established practices.
6. **Pre-commit hooks:** Strongly recommend setting up a `.pre-commit-config.yaml` in the monorepo root to enforce code style and quality automatically. Example:

    ```yaml
    # .pre-commit-config.yaml (in monorepo root)
    repos:
    -   repo: https://github.com/pre-commit/pre-commit-hooks
        rev: v4.5.0
        hooks:
        -   id: trailing-whitespace
        -   id: end-of-file-fixer
        -   id: check-yaml
        -   id: check-added-large-files
    -   repo: https://github.com/psf/black
        rev: 23.12.1 # Or your team's preferred version
        hooks:
        -   id: black
            args: ["./components/moremi_biokit/"] # Target specific component
    -   repo: https://github.com/PyCQA/flake8
        rev: 7.0.0 # Or your team's preferred version
        hooks:
        -   id: flake8
            args: ["./components/moremi_biokit/"] # Target specific component
    # Add mypy if you use it
    # -   repo: https://github.com/pre-commit/mirrors-mypy
    #     rev: 'v1.8.0'
    #     hooks:
    #     -   id: mypy
    #         args: ["./components/moremi_biokit/"]
    #         additional_dependencies: [ "pandas-stubs", ...] # for typed libraries
    ```

    Developers would then run `pre-commit install` once after cloning.
