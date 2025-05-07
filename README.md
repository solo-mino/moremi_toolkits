
# Moremi Toolkits Monorepo

## Overview

Welcome to the `moremi_toolkits` monorepo! This repository serves as the central location for foundational Python packages used within the Moremi ecosystem, primarily supporting the **Moremi Bio Autonomous Agent**.

The primary package housed here is **`moremi-biokit`**, a comprehensive toolkit for processing and analyzing biological entities, currently including:
-   Small molecule analysis (via SMILES)
-   Antibody analysis

This monorepo structure allows for unified versioning, issue tracking, and development of these core components.

**Key Distinction:**
-   **Moremi Bio (Autonomous Agent):** The higher-level AI-driven system that orchestrates research tasks (lives in its own separate project/repository).
-   **`moremi-biokit` (This Toolkit):** An installable Python library providing the underlying processing functions, classes, and utilities that the Moremi Bio agent calls upon. It is *not* autonomous on its own.

## Repository Structure

```
moremi_toolkits/
├── .git/
├── components/
│   └── moremi_biokit/               # Source and build configuration for the 'moremi-biokit' package
│       ├── moremi_biokit/           # The actual Python package source code
│       │   ├── __init__.py
│       │   ├── smiles/              # Sub-package for SMILES tools
│       │   │   ├── __init__.py
│       │   │   ├── batch_processor.py
│       │   │   ├── validator.py
│       │   │   ├── ranker.py
│       │   │   └── property_calculators/
│       │   │   └── reports/
│       │   └── antibodies/          # Sub-package for Antibody tools
│       │       ├── __init__.py
│       │       └── ...
│       ├── pyproject.toml           # Build config for moremi-biokit
│       ├── README.md                # Detailed README for the moremi-biokit package
│       └── tests/                   # Unit tests for moremi_biokit
│
├── .gitignore
├── LICENSE                        # (Ensure you have an appropriate internal license)
└── README.md                      # This file (main monorepo README)
```

## Getting Started

### Prerequisites

-   Python (version specified in `components/moremi_biokit/pyproject.toml`, e.g., >=3.8)
-   `pip` and `build` (for building/installing the package)
-   Git

### For Developers of `moremi-biokit`

1.  **Clone the repository:**
    ```bash
    git clone <your_internal_git_repo_url>/moremi_toolkits.git
    cd moremi_toolkits
    ```

2.  **Set up a virtual environment:**
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```

3.  **Install in editable mode (for development):**
    This allows you to make changes to the `moremi_biokit` code and have them immediately available without reinstalling.
    ```bash
    pip install -e ./components/moremi_biokit/
    ```

4.  **Install development dependencies (e.g., for testing, linting):**
    (Assuming these are defined in `pyproject.toml` under `[project.optional-dependencies]`)
    ```bash
    pip install -e "./components/moremi_biokit/[dev]"
    ```

5.  **Run tests:**
    Navigate to the `components/moremi_biokit/` directory or configure your test runner.
    ```bash
    # Example if using pytest (install it first if not in dev dependencies)
    pytest ./components/moremi_biokit/tests/
    ```

### For Consumers of `moremi-biokit` (e.g., the Moremi Bio Agent project)

The `moremi-biokit` package should be installed as a dependency.

1.  **Add to `requirements.txt` (or equivalent for your dependency manager):**

    *   **From a specific tag (recommended for stable builds):**
        ```
        moremi-biokit @ git+<your_internal_git_repo_url>/moremi_toolkits.git@v0.1.0#subdirectory=components/moremi_biokit
        ```
        (Replace `v0.1.0` with the desired tag and `<your_internal_git_repo_url>` with the actual URL)

    *   **From a specific branch (e.g., `main` - for development or cutting-edge builds):**
        ```
        moremi-biokit @ git+<your_internal_git_repo_url>/moremi_toolkits.git@main#subdirectory=components/moremi_biokit
        ```

2.  **Install dependencies in the consuming project:**
    ```bash
    pip install -r requirements.txt
    ```

## Contribution Guidelines

(Detail your internal contribution process here)
-   Branching strategy (e.g., `feature/`, `bugfix/`, `chore/`)
-   Commit message conventions
-   Code review process
-   Testing requirements
-   Code style (e.g., Black, Flake8 - consider pre-commit hooks)

## License

This project is proprietary and licensed under the company's internal licensing terms. Please refer to the `LICENSE` file for details. (Ensure you create a `LICENSE` file appropriate for internal use).

## Contact

-   For questions about `moremi-biokit` development: [Lead Developer/Team Contact]
-   For issues, please use the issue tracker in this Git repository.
