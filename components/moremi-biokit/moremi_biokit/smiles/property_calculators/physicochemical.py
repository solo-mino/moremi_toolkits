"""
Molecular Physicochemical Properties

References:
    - TPSA: Ertl P, Rohde B, Selzer P. Fast calculation of molecular polar surface area 
            as a sum of fragment-based contributions and its application to the prediction 
            of drug transport properties. J Med Chem. 2000;43(20):3714-3717.
    - Molar Refractivity: Ghose AK, Crippen GM. Atomic physicochemical parameters for 
                         three-dimensional structure-directed quantitative structure-activity 
                         relationships. J Comput Chem. 1987;8(6):1109-1126.
    - Fraction CSP3: Lovering F, Bikker J, Humblet C. Escape from flatland: increasing 
                    saturation as an approach to improving clinical success. 
                    J Med Chem. 2009;52(21):6752-6756.
"""

from typing import Dict, Optional, Union
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


class PhysiochemicalProperties:
    """
    A class to calculate and store physiochemical properties of molecules using SMILES notation.
    
    This class provides methods to calculate various molecular properties that are commonly
    used in drug discovery and molecular analysis, matching SwissADME calculations.
    
    Attributes:
        smile (str): The SMILES string representation of the molecule
        molecule (rdkit.Chem.rdchem.Mol): The RDKit molecule object
        _properties_cache (Dict): Cache for storing calculated properties
    
    Example:
        >>> mol = PhysiochemicalProperties("CC(=O)OC1=CC=CC=C1C(=O)O")
        >>> properties = mol.calculate_physiochemical_properties()
        >>> report = mol.generate_report()
        >>> print(report)
    """
    
    def __init__(self, smile: str):
        """
        Initialize the PhysiochemicalProperties object.
        
        Args:
            smile (str): SMILES string representation of the molecule
            
        Raises:
            ValueError: If the SMILES string is invalid or cannot be parsed
        """
        self.smile = smile
        self._properties_cache = {}
        self._initialize_molecule()
    
    def _initialize_molecule(self) -> None:
        """
        Initialize the RDKit molecule object from SMILES string.
        
        Raises:
            ValueError: If the SMILES string is invalid or cannot be parsed
        """
        self.molecule = Chem.MolFromSmiles(self.smile)
        if self.molecule is None:
            raise ValueError(f"Invalid SMILES string: {self.smile}") 
            
    def _get_cached_properties(self) -> Optional[Dict]:
        """
        Retrieve cached properties if available.
        
        Returns:
            Dict or None: Cached properties if available, None otherwise
        """
        return self._properties_cache.get('basic_properties')
    
    def calculate_physiochemical_properties(self) -> Dict[str, Union[float, int, str]]:
        """
        Calculate basic molecular properties matching SwissADME calculations.
        
        This method calculates various molecular properties including:
        - Molecular Formula
        - Molecular Weight (MW)
        - Number of Heavy Atoms
        - Number of Aromatic Heavy Atoms
        - Fraction of sp3 Carbons
        - Number of Rotatable Bonds
        - Number of H-bond Donors (HBD)
        - Number of H-bond Acceptors (HBA)
        - Molar Refractivity
        - Topological Polar Surface Area (TPSA)
        
        Returns:
            Dict[str, Union[float, int, str]]: Dictionary containing calculated properties
            
        Note:
            Results are cached for subsequent retrievals
        """
        # Check cache first
        cached_properties = self._get_cached_properties()
        if cached_properties:
            return cached_properties
            
        try:
            properties = {}
            
            # Molecular Formula
            properties['formula'] = rdMolDescriptors.CalcMolFormula(self.molecule)
            
            # Molecular Weight - using MonoisotopicMolWt for better precision
            properties['MW'] = round(Descriptors.ExactMolWt(self.molecule), 2)
            
            # Heavy atoms
            properties['num_heavy_atoms'] = self.molecule.GetNumHeavyAtoms()
            
            # Aromatic heavy atoms
            properties['num_aromatic_heavy_atoms'] = len([atom for atom in self.molecule.GetAtoms() 
                                                        if atom.GetIsAromatic()])
            
            # Fraction of sp3 hybridized carbons
            properties['fraction_csp3'] = round(Descriptors.FractionCSP3(self.molecule), 2)
            
            # Rotatable bonds
            properties['num_rotatable_bonds'] = rdMolDescriptors.CalcNumRotatableBonds(self.molecule)
            
            # H-bond donors and acceptors (matching SwissADME counting)
            properties['num_hbd'] = rdMolDescriptors.CalcNumHBD(self.molecule)
            properties['num_hba'] = rdMolDescriptors.CalcNumHBA(self.molecule)
            
            # Molar Refractivity
            properties['molar_refractivity'] = round(Crippen.MolMR(self.molecule), 2)
            
            # TPSA - using only polar surface area from N and O
            properties['TPSA'] = round(Descriptors.TPSA(self.molecule, includeSandP=False), 2)
            
            # Cache the results
            self._properties_cache['basic_properties'] = properties
            
            return properties
            
        except Exception as e:
            raise RuntimeError(f"Error calculating properties: {str(e)}")
    
    
    def get_smile(self) -> str:
        """
        Get the SMILES string representation of the molecule.
        
        Returns:
            str: SMILES string
        """
        return self.smile
    
    def get_molecule(self) -> Chem.rdchem.Mol:
        """
        Get the RDKit molecule object.
        
        Returns:
            rdkit.Chem.rdchem.Mol: RDKit molecule object
        """
        return self.molecule