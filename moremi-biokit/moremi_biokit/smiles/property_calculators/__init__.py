"""
Subpackage for SMILES property calculators.

This subpackage provides a collection of modules, each dedicated to calculating
a specific set of molecular properties from SMILES strings. These calculators
are used by the `SMILESValidator` to gather comprehensive data for 
molecular analysis.

Available Calculators:
    - ADMETPredictor: Predicts Absorption, Distribution, Metabolism, Excretion, 
      and Toxicity properties using deep learning models and RDKit.
    - DruglikenessProperties: Calculates various drug-likeness rules (Lipinski, 
      Ghose, Veber, Egan, Muegge) and bioavailability scores.
    - LipophilicityProperties: Computes different LogP values (iLOGP, XLOGP3, 
      WLOGP, MLOGP, SILICOS-IT) and a consensus LogP.
    - MedicinalChemistryProperties: Assesses medicinal chemistry aspects like PAINS/Brenk
      alerts, lead-likeness, and synthetic accessibility.
    - PharmacokineticsProperties: Estimates pharmacokinetic parameters such as GI 
      absorption, BBB permeability, P-gp substrate status, CYP inhibition, and LogKp.
    - PhysiochemicalProperties: Determines fundamental physicochemical properties like
      molecular weight, TPSA, H-bond donors/acceptors, rotatable bonds, etc.
    - SolubilityProperties: Calculates water solubility using methods like ESOL, Ali, 
      and SILICOS-IT.

Example:
    To use a specific calculator, you can import it directly after importing 
    the `property_calculators` module or its parent `smiles` module.

    >>> from moremi_biokit.smiles.property_calculators import PhysiochemicalProperties
    >>> aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    >>> phys_chem_calc = PhysiochemicalProperties(aspirin_smiles)
    >>> properties = phys_chem_calc.calculate_physiochemical_properties()
    >>> print(f"Aspirin MW: {properties.get('MW')}")

    Alternatively, these calculators are typically orchestrated by the 
    `moremi_biokit.smiles.validator.SMILESValidator`.
"""

from .admet_predictor import ADMETPredictor
from .druglikeness import DruglikenessProperties
from .lipophilicity import LipophilicityProperties
from .medicinal_chemistry import MedicinalChemistryProperties
from .pharmacokinetics import PharmacokineticsProperties
from .physicochemical import PhysiochemicalProperties
from .solubility import SolubilityProperties

__all__ = [
    "ADMETPredictor",
    "DruglikenessProperties",
    "LipophilicityProperties",
    "MedicinalChemistryProperties",
    "PharmacokineticsProperties",
    "PhysiochemicalProperties",
    "SolubilityProperties",
]
