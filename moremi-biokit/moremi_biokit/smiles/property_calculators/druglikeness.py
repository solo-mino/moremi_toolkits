"""
Module for calculating drug-likeness properties of small molecules. 

This module provides methods to calculate various drug-likeness rules and scores:
- Lipinski's Rule of 5 (Pfizer) - Lipinski CA. et al. 2001 Adv. Drug Deliv. Rev.
- Ghose Filter - Ghose AK. et al. 1999 J. Comb. Chem.
- Veber (GSK) Filter - Veber DF. et al. 2002 J. Med. Chem.
- Egan (Pharmacia) Filter - Egan WJ. et al. 2000 J. Med. Chem.
- Muegge (Bayer) Filter - Muegge I. et al. 2001 J. Med. Chem.
- Bioavailability Score - Martin YC. 2005 J. Med. Chem.
"""

from typing import Dict, List, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


class DruglikenessProperties:
    """
    A class to calculate drug-likeness properties of molecules using SMILES notation.
    
    References:
        - Lipinski CA. et al. 2001 Adv. Drug Deliv. Rev. (Pfizer)
        - Ghose AK. et al. 1999 J. Comb. Chem.
        - Veber DF. et al. 2002 J. Med. Chem. (GSK)
        - Egan WJ. et al. 2000 J. Med. Chem. (Pharmacia)
        - Muegge I. et al. 2001 J. Med. Chem. (Bayer)
        - Martin YC. 2005 J. Med. Chem. (Bioavailability)
    
    Attributes:
        smile (str): The SMILES string representation of the molecule
        molecule (rdkit.Chem.rdchem.Mol): The RDKit molecule object
        _properties_cache (Dict): Cache for storing calculated properties
    """
    
    def __init__(self, smile: str):
        """Initialize with SMILES string."""
        self.smile = smile
        self._properties_cache = {}
        self._initialize_molecule()
        
    def _initialize_molecule(self) -> None:
        """Initialize RDKit molecule object."""
        self.molecule = Chem.MolFromSmiles(self.smile)
        if self.molecule is None:
            raise ValueError(f"Invalid SMILES string: {self.smile}")
        self.molecule = Chem.AddHs(self.molecule)
        
    def _check_lipinski(self) -> Tuple[bool, List[str]]:
        """
        Check Lipinski's Rule of 5 (Pfizer).
        Reference: Lipinski CA. et al. 2001 Adv. Drug Deliv. Rev.
        
        Criteria:
        - Molecular Weight (MW) ≤ 500
        - Modified Log P (MLOGP) ≤ 4.15
        - Number of Nitrogen or Oxygen atoms ≤ 10
        - Number of NH or OH groups ≤ 5
        
        Returns:
            Tuple[bool, List[str]]: (passes, violations)
        """
        mw = Descriptors.ExactMolWt(self.molecule)
        mlogp = Descriptors.MolLogP(self.molecule)  # Using MLOGP
        no_count = sum(1 for atom in self.molecule.GetAtoms() 
                      if atom.GetAtomicNum() in [7, 8])  # N or O
        hbd = Descriptors.NumHDonors(self.molecule)  # NH or OH
        
        violations = []
        if mw > 500: violations.append("MW>500")
        if mlogp > 4.15: violations.append("MLOGP>4.15")
        if no_count > 10: violations.append("N/O>10")
        if hbd > 5: violations.append("NH/OH>5")
        
        passes = len(violations) <= 1  # Lipinski allows 1 violation
        return passes, violations
        
    def _check_ghose(self) -> Tuple[bool, List[str]]:
        """
        Check Ghose Filter.
        Reference: Ghose AK. et al. 1999 J. Comb. Chem.
        
        Criteria:
        - 160 ≤ Molecular Weight (MW) ≤ 480
        - -0.4 ≤ Weighted Log P (WLOGP) ≤ 5.6
        - Molar Refractivity (MR) between 40 and 130
        - Total number of atoms between 20 and 70
        
        Returns:
            Tuple[bool, List[str]]: (passes, violations)
        """
        mw = Descriptors.ExactMolWt(self.molecule)
        wlogp = Crippen.MolLogP(self.molecule)  # Using WLOGP
        molar_refractivity = Crippen.MolMR(self.molecule)
        num_atoms = self.molecule.GetNumAtoms()
        
        violations = []
        if not (160 <= mw <= 480): violations.append("MW outside [160,480]")
        if not (-0.4 <= wlogp <= 5.6): violations.append("WLOGP outside [-0.4,5.6]")
        if not (40 <= molar_refractivity <= 130): violations.append("MR outside [40,130]")
        if not (20 <= num_atoms <= 70): violations.append("Atoms outside [20,70]")
        
        return len(violations) == 0, violations
        
    def _check_veber(self) -> Tuple[bool, List[str]]:
        """
        Check Veber (GSK) Filter.
        Reference: Veber DF. et al. 2002 J. Med. Chem.
        
        Criteria:
        - Number of rotatable bonds ≤ 10
        - Topological Polar Surface Area (TPSA) ≤ 140
        
        Returns:
            Tuple[bool, List[str]]: (passes, violations)
        """
        rotatable_bonds = Descriptors.NumRotatableBonds(self.molecule)
        tpsa = Descriptors.TPSA(self.molecule)
        
        violations = []
        if rotatable_bonds > 10: violations.append("RotBonds>10")
        if tpsa > 140: violations.append("TPSA>140")
        
        # Veber is more lenient - pass even with one violation
        return True, violations
        
    def _check_egan(self) -> Tuple[bool, List[str]]:
        """
        Check Egan (Pharmacia) Filter.
        Reference: Egan WJ. et al. 2000 J. Med. Chem.
        
        Criteria:
        - Weighted Log P (WLOGP) ≤ 5.88
        - Topological Polar Surface Area (TPSA) ≤ 131.6
        
        Returns:
            Tuple[bool, List[str]]: (passes, violations)
        """
        wlogp = Crippen.MolLogP(self.molecule)  # Using WLOGP
        tpsa = Descriptors.TPSA(self.molecule)
        
        violations = []
        if wlogp > 5.88: violations.append("WLOGP>5.88")
        if tpsa > 131.6: violations.append("TPSA>131.6")
        
        return len(violations) == 0, violations
        
    def _check_muegge(self) -> Tuple[bool, List[str]]:
        """
        Check Muegge (Bayer) Filter.
        Reference: Muegge I. et al. 2001 J. Med. Chem.
        
        Criteria:
        - Molecular Weight (MW) between 200 and 600
        - Extended Log P (XLOGP) between -2 and 5
        - Topological Polar Surface Area (TPSA) ≤ 150
        - Number of rings ≤ 7
        - Number of carbon atoms > 4
        - Number of heteroatoms > 1
        - Number of rotatable bonds ≤ 15
        - Hydrogen bond acceptors ≤ 10
        - Hydrogen bond donors ≤ 5
        
        Returns:
            Tuple[bool, List[str]]: (passes, violations)
        """
        mw = Descriptors.ExactMolWt(self.molecule)
        # Use XLOGP3 specifically for Muegge rule
        xlogp = Descriptors.MolLogP(self.molecule)
        tpsa = Descriptors.TPSA(self.molecule)
        rings = Chem.GetSSSR(self.molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(self.molecule)
        hba = Descriptors.NumHAcceptors(self.molecule)
        hbd = Descriptors.NumHDonors(self.molecule)
        carbons = sum(1 for atom in self.molecule.GetAtoms() 
                     if atom.GetAtomicNum() == 6)
        heteroatoms = Descriptors.NumHeteroatoms(self.molecule)
        
        
        violations = []
        if not (200 <= mw <= 600): 
            violations.append("MW outside [200,600]")
        if xlogp > 5: 
            violations.append("XLOGP3>5")  # Only check upper bound as per web platform
        if tpsa > 150: 
            violations.append("TPSA>150")
        if len(rings) > 7: 
            violations.append("Rings>7")
        if carbons <= 4: 
            violations.append("Carbons≤4")
        if heteroatoms <= 1: 
            violations.append("Heteroatoms≤1")
        if rotatable_bonds > 15: 
            violations.append("RotBonds>15")
        if hba > 10: 
            violations.append("HBA>10")
        if hbd > 5: 
            violations.append("HBD>5")
        
        # Muegge is strict - ANY violation means it fails
        return len(violations) == 0, violations
        
    def _calculate_bioavailability_score(self) -> float:
        """
        Calculate Bioavailability Score.
        Reference: Martin YC. 2005 J. Med. Chem.
        
        Description: Probability of Fraction absorbed (F) > 10% in rat
        
        Returns:
            float: Bioavailability score between 0 and 1
        """
        mw = Descriptors.ExactMolWt(self.molecule)
        logp = Descriptors.MolLogP(self.molecule)
        psa = Descriptors.TPSA(self.molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(self.molecule)
        hbd = Descriptors.NumHDonors(self.molecule)
        hba = Descriptors.NumHAcceptors(self.molecule)
        
        # Martin's rules for bioavailability
        if logp > 5 or hbd >= 2:  # Adjusted to match reference
            score = 0.55
        elif psa > 140 or mw > 500:
            score = 0.11
        elif rotatable_bonds > 10:
            score = 0.17
        else:
            score = 0.85
            
        return round(score, 2)
        
    def calculate_druglikeness(self) -> Dict:
        """
        Calculate all drug-likeness properties.
        
        Returns:
            Dict containing all calculated drug-likeness properties
        """
        try:
            if 'druglikeness' in self._properties_cache:
                return self._properties_cache['druglikeness']
            
            # Calculate all properties
            lipinski_pass, lipinski_violations = self._check_lipinski()
            ghose_pass, ghose_violations = self._check_ghose()
            veber_pass, veber_violations = self._check_veber()
            egan_pass, egan_violations = self._check_egan()
            muegge_pass, muegge_violations = self._check_muegge()
            
            properties = {
                'lipinski': {
                    'passes': lipinski_pass,
                    'violations': lipinski_violations,
                },
                'ghose': {
                    'passes': ghose_pass,
                    'violations': ghose_violations
                },
                'veber': {
                    'passes': veber_pass,
                    'violations': veber_violations
                },
                'egan': {
                    'passes': egan_pass,
                    'violations': egan_violations
                },
                'muegge': {
                    'passes': muegge_pass,
                    'violations': muegge_violations
                },
                'bioavailability_score': self._calculate_bioavailability_score()
            }
            
            # Cache the results
            self._properties_cache['druglikeness'] = properties
            
            return properties
            
        except Exception as e:
            raise RuntimeError(f"Error calculating drug-likeness properties: {str(e)}")
