# Moremi Toolkits Monorepo 🧬

![Version](https://img.shields.io/badge/Version-0.1.0-blue)
![Status](https://img.shields.io/badge/Status-DRAFT-orange)

## 1. Overview 📋

Welcome to the `moremi_toolkits` monorepo! This repository serves as the central location for foundational Python packages used within the Moremi ecosystem, primarily supporting the **Moremi Bio Autonomous Agent**.

The primary package housed here is **`moremi-biokit`**, a comprehensive toolkit for processing and analyzing biological entities, currently including:

- 🧪 Small molecule analysis (via SMILES)
- 🔬 Antibody analysis

This monorepo structure allows for unified versioning, issue tracking, and development of these core components.

### 1.1. Key Distinction

- **Moremi Bio (Autonomous Agent):** The higher-level AI-driven system that orchestrates research tasks (lives in its own separate project/repository).
- **`moremi-biokit` (This Toolkit):** An installable Python library providing the underlying processing functions, classes, and utilities that the Moremi Bio agent calls upon. It is *not* autonomous on its own.

## 2. Important Requirements ⚠️

- **Linux Environment Required:** Due to the nature of some core dependencies (e.g., Open Babel), `moremi-biokit` is intended to be run in a Linux environment. For Windows users, this means using Windows Subsystem for Linux (WSL).
- **Conda for Core Dependencies:** Packages like `openbabel` have complex binary dependencies. It is **strongly recommended** to install them using Conda *before* installing `moremi-biokit` with pip. This will avoid common build failures.

## 3. Repository Structure 📁

```plaintext
moremi_toolkits/
├── .git/                        # Git version control directory
├── .gitignore                   # Git ignore rules
├── README.md                    # Main project documentation (this file)
├── requirements.txt             # Top-level requirements (may be used for dev or meta-deps)
├── docs/                        # Project documentation
│   ├── main_docs.md
│   ├── smiles_api.md
│   └── coding_guide.md
├── notebooks/                   # Example and exploratory Jupyter notebooks
│   ├── smiles_usage_example.ipynb
│   └── proteins_usage_example.ipynb
├── moremi-biokit/               # Main Python package source and build config
│   ├── pyproject.toml           # Build config for 'moremi-biokit'
│   ├── tests/                   # Unit tests for 'moremi_biokit'
│   │   ├── __init__.py
│   │   ├── conftest.py
│   │   ├── proteins/
│   │   └── smiles/
│   ├── moremi_biokit/           # The actual Python package
│   │   ├── __init__.py
│   │   ├── pdb_fetcher.py
│   │   ├── assets/
│   │   │   ├── minologo.png
│   │   │   └── fonts/
│   │   ├── connectors/
│   │   │   ├── __init__.py
│   │   │   ├── rcsb.py
│   │   │   ├── local_pdb_fasta_parser.py
│   │   │   ├── uniprot.py
│   │   │   ├── ncbi.py
│   │   │   └── _utils.py
│   │   ├── proteins/
│   │   │   ├── __init__.py
│   │   │   ├── batch_protein_processor.py
│   │   │   ├── protein_ranker.py
│   │   │   ├── protein_validator.py
│   │   │   ├── utils/
│   │   │   ├── reports/
│   │   │   └── analysis_tools/
│   │   ├── smiles/
│   │   │   ├── __init__.py
│   │   │   ├── batch_molecule_processor.py
│   │   │   ├── small_molecule_validator_v3.py
│   │   │   ├── small_molecule_ranker_v4.py
│   │   │   ├── utils/
│   │   │   ├── reports/
│   │   │   └── property_calculators/
│   │   └── test_pdb_targets/
│   │       ├── __init__.py
│   │       ├── proteins/
│   │       └── smiles/
```

## 4. Getting Started 🚀

### 4.1. Prerequisites

- 🐧 A Linux environment (e.g., native Linux or WSL on Windows)
- 🐍 Miniconda or Anaconda (for managing environments and complex dependencies)
- 📦 Python (version specified in `moremi-biokit/pyproject.toml`, e.g., >=3.8 - Conda will handle this)
- 📂 Git

### 4.2. Setting up SSH Keys for GitHub (for private repositories) 🔑

If you need to clone or install from this private repository, you'll need to authenticate with GitHub, typically using SSH keys.

1. **Check for existing SSH keys:**

   ```bash
   ls -al ~/.ssh
   ```

   Look for files like `id_rsa.pub` or `id_ed25519.pub`. If they exist, you might already have a key.

2. **Generate a new SSH key if you don't have one:**

   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   ```

   Follow the prompts. It's generally fine to accept the default file location and skip setting a passphrase for local development convenience, or add one for more security.

3. **Add your SSH private key to the ssh-agent:**

   ```bash
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519  # Or id_rsa if you generated an RSA key
   ```

4. **Add your SSH public key to your GitHub account:**
   Copy the contents of your public key file:

   ```bash
   cat ~/.ssh/id_ed25519.pub  # Or id_rsa.pub
   ```

   Go to your GitHub account settings -> SSH and GPG keys -> New SSH key. Paste your key there and save it.

5. **Test your SSH connection:**

   ```bash
   ssh -T git@github.com
   ```

   You should see a message like: "Hi username! You've successfully authenticated, but GitHub does not provide shell access."

## 5. Installation & Development 💻

### 5.1. For Developers of `moremi-biokit`

1. **Clone the repository (using SSH):**

   ```bash
   git clone git@github.com:solo-mino/moremi_toolkits.git
   cd moremi_toolkits
   ```

2. **Create and activate a Conda environment:**
   It's recommended to install `openbabel` and `rdkit` (if used by your project, it's listed in `dependencies`) via Conda first.

   ```bash
   conda create -n moremi_env python=3.10  # Or your desired Python version (e.g., 3.8, 3.9)
   conda activate moremi_env
   conda install -c conda-forge openbabel rdkit --yes # Add other Conda-specific deps here if any
   ```

   *Note: Check the `moremi-biokit/pyproject.toml` for the specific Python version required.*

3. **Install `moremi-biokit` in editable mode:**
   This allows you to make changes to the code and have them immediately available.

   ```bash
   pip install -e ./moremi-biokit/
   ```

4. **Install development dependencies:**
   (Assuming these are defined in `moremi-biokit/pyproject.toml` under `[project.optional-dependencies]`)

   ```bash
   pip install -e "./moremi-biokit/[dev]"
   ```

5. **Run tests:**

   ```bash
   pytest ./moremi-biokit/tests/
   ```

### 5.2. For Consumers of `moremi-biokit` (e.g., the Moremi Bio Agent project)

1. **Set up a Conda environment with necessary dependencies:**
   It is **highly recommended** to install `openbabel` and `rdkit` via Conda into the target environment *before* installing `moremi-biokit`.

   ```bash
   conda create -n my_agent_env python=3.10 # Or your desired Python version
   conda activate my_agent_env
   conda install -c conda-forge openbabel
   ```

2. **Pip install:**

   - **From a specific branch (e.g., `main`, using SSH URL for private repo):**

    ```bash
    pip install git+ssh://git@github.com/solo-mino/moremi_toolkits.git@main#subdirectory=moremi-biokit
    ```

3. **Install dependencies in the consuming project (after activating the Conda environment):**

   ```bash
   pip install -r requirements.txt
   ```

   Pip should now find `openbabel` and `rdkit` already installed via Conda and will proceed to install `moremi-biokit` and its other pure Python dependencies.

   **NB:** This step is *optional* as the package in the `step 2` will install all the needed dependencies.

## 6. Contribution Guidelines 🤝

(Detail your internal contribution process here)

- Branching strategy (e.g., `feature/`, `bugfix/`, `chore/`)
- Commit message conventions
- Code review process
- Testing requirements
- Code style (e.g., Black, Flake8 - consider pre-commit hooks)

## 7. License 📄

This project is proprietary and licensed under the company's internal licensing terms. Please refer to the `LICENSE` file for details. (Ensure you create a `LICENSE` file appropriate for internal use).

## 8. Contact 📧

- For questions about `moremi-biokit` development: [Lead Developer/Team Contact]
- For issues, please use the issue tracker in this Git repository.
