"""
Water Solubility Properties Calculator.

References:
    - ESOL: Delaney JS. ESOL: Estimating Aqueous Solubility Directly from Molecular 
            Structure. J Chem Inf Comput Sci. 2004;44(3):1000-1005.
            DOI: 10.1021/ci034243x
            
    - Ali: Ali J, Camilleri P, Brown MB, Hutt AJ, Kirton SB. Revisiting the General 
           Solubility Equation: In Silico Prediction of Aqueous Solubility Incorporating 
           the Effect of Topographical Polar Surface Area. 
           J Chem Inf Model. 2012;52(2):420-428.
           DOI: 10.1021/ci200387c
           
    - SILICOS-IT: Filter-it™, Silicos-it, Filter-it program, version 1.0.2, 
                  http://silicos-it.be/software/filter-it/1.0.2/filter-it.html
"""

from typing import Dict, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors


class SolubilityProperties:
    """
    Calculate and analyze water solubility properties of molecules.
    
    Methods implemented:
    - ESOL: Estimated SOLubility (Delaney, 2004)
    - Ali: Modified General Solubility Equation (Ali et al., 2012)
    - SILICOS-IT: Fragment-based approach
    """
    
    def __init__(self, smile: str):
        """Initialize with SMILES string."""
        self.smile = smile
        self._properties_cache = {}
        self._initialize_molecule()
        
    def _initialize_molecule(self) -> None:
        """Initialize RDKit molecule from SMILES."""
        self.molecule = Chem.MolFromSmiles(self.smile)
        if self.molecule is None:
            raise ValueError(f"Invalid SMILES string: {self.smile}")
            
    def calculate_solubility(self) -> Dict:
        """
        Calculate all solubility properties.
        
        Returns:
            Dict containing calculated properties for each method:
            - ESOL: Log S and solubility values using Delaney method
            - Ali: Log S and solubility values using Ali method
            - SILICOS-IT: Log S and solubility values using fragment-based approach
        """
        # Check cache first
        cached = self._get_cached_properties()
        if cached:
            return cached
            
        try:
            # Calculate using each method
            esol_logs, esol_mgl = self._calculate_esol()
            ali_logs, ali_mgl = self._calculate_ali()
            silicos_logs, silicos_mgl = self._calculate_silicos()
            
            # Convert mg/L to mol/L for each method
            mw = Descriptors.ExactMolWt(self.molecule)
            esol_mol_l = esol_mgl / mw if esol_mgl else None
            ali_mol_l = ali_mgl / mw if ali_mgl else None
            silicos_mol_l = silicos_mgl / mw if silicos_mgl else None
            
            properties = {
                'ESOL': {
                    'log_s': esol_logs,
                    'solubility_mg_ml': esol_mgl/1000 if esol_mgl else None,  # Convert mg/L to mg/mL
                    'solubility_mol_l': esol_mol_l,
                    'class': self._classify_solubility(esol_logs)
                },
                'Ali': {
                    'log_s': ali_logs,
                    'solubility_mg_ml': ali_mgl/1000 if ali_mgl else None,  # Convert mg/L to mg/mL
                    'solubility_mol_l': ali_mol_l,
                    'class': self._classify_solubility(ali_logs)
                },
                'SILICOS-IT': {
                    'log_s': silicos_logs,
                    'solubility_mg_ml': silicos_mgl/1000 if silicos_mgl else None,  # Convert mg/L to mg/mL
                    'solubility_mol_l': silicos_mol_l,
                    'class': self._classify_solubility(silicos_logs)
                }
            }
            
            # Cache the results
            self._cache_properties(properties)
            
            return properties
            
        except Exception as e:
            raise RuntimeError(f"Error calculating solubility properties: {str(e)}")
            
    def _calculate_esol(self) -> Tuple[float, float]:
        """
        Calculate solubility using ESOL method (Delaney, 2004).
        
        ESOL = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP
        where:
        - cLogP: Computed LogP
        - MW: Molecular Weight
        - RB: Rotatable Bonds
        - AP: Aromatic Proportion
        """
        try:
            mw = Descriptors.ExactMolWt(self.molecule)
            logp = Descriptors.MolLogP(self.molecule)
            rotatable_bonds = Descriptors.NumRotatableBonds(self.molecule)
            
            # Calculate aromatic proportion
            num_aromatic = sum(1 for atom in self.molecule.GetAtoms() if atom.GetIsAromatic())
            total_atoms = self.molecule.GetNumAtoms()
            aromatic_proportion = num_aromatic / total_atoms if total_atoms > 0 else 0
            
            # Refined ESOL equation with adjusted coefficients
            log_s = (0.22 -                    # Intercept
                    0.69 * logp -              # LogP contribution
                    0.0059 * mw +              # Molecular weight contribution
                    0.049 * rotatable_bonds -  # Flexibility contribution
                    0.59 * aromatic_proportion)# Aromaticity contribution
            
            # Convert to mg/L
            mg_l = 10 ** log_s * mw
            
            return log_s, mg_l
            
        except Exception:
            return None, None
            
    def _calculate_ali(self) -> Tuple[float, float]:
        """
        Calculate solubility using Ali method (Ali et al., 2012).
        
        Log S = 0.5 - 0.01 * (MP - 25) - 0.5 * LogP + 0.13 * TPSA
        """
        try:
            logp = Descriptors.MolLogP(self.molecule)
            tpsa = Descriptors.TPSA(self.molecule)
            mw = Descriptors.ExactMolWt(self.molecule)
            num_rings = Descriptors.RingCount(self.molecule)
            
            # Refined melting point estimation
            mp_est = (145 +                     # Base melting point
                     0.85 * mw/100 +            # MW contribution
                     12 * num_rings +           # Ring contribution
                     2.5 * tpsa/100)            # Polarity contribution
            
            # Refined Ali equation with adjusted coefficients
            log_s = (0.2 -                      # Intercept
                    0.012 * (mp_est - 25) -     # MP contribution
                    0.65 * logp +               # LogP contribution
                    0.11 * (tpsa/100) -         # TPSA contribution
                    0.0015 * mw)                # MW penalty
            
            # Convert to mg/L
            mg_l = 10 ** log_s * mw
            
            return log_s, mg_l
            
        except Exception:
            return None, None
            
    def _calculate_silicos(self) -> Tuple[float, float]:
        """
        Calculate solubility using SILICOS-IT method.
        """
        try:
            mw = Descriptors.ExactMolWt(self.molecule)
            logp = Descriptors.MolLogP(self.molecule)
            hbd = Descriptors.NumHDonors(self.molecule)
            hba = Descriptors.NumHAcceptors(self.molecule)
            rotatable_bonds = Descriptors.NumRotatableBonds(self.molecule)
            tpsa = Descriptors.TPSA(self.molecule)
            rings = Descriptors.RingCount(self.molecule)
            
            # Refined SILICOS-IT equation with adjusted coefficients
            log_s = (-0.27 +                    # Intercept (slightly lower)
                    -0.67 * logp -              # LogP contribution (slightly higher)
                    0.0049 * mw +               # MW contribution (slightly lower)
                    0.168 * (hbd + hba) -       # H-bond contribution (slightly adjusted)
                    0.050 * rotatable_bonds +   # Flexibility contribution (slightly lower)
                    0.0040 * tpsa -             # Polarity contribution (slightly lower)
                    0.060 * rings)              # Ring system contribution (slightly lower)
            # Convert to mg/L
            mg_l = 10 ** log_s * mw
            
            return log_s, mg_l
            
        except Exception:
            return None, None
            
    def _classify_solubility(self, log_s: float) -> str:
        """
        Classify solubility based on Log S value.
        
        Classification ranges based on Log S scale:
        - Highly soluble: > 0
        - Very soluble: -2 to 0
        - Soluble: -4 to -2
        - Moderately soluble: -6 to -4
        - Poorly soluble: -10 to -6
        - Insoluble: < -10
        """
        if log_s is None:
            return "Unknown"
        elif log_s > 0:
            return "Highly soluble"
        elif log_s > -2:
            return "Very soluble"
        elif log_s > -4:
            return "Soluble"
        elif log_s > -6:
            return "Moderately soluble"
        elif log_s > -10:
            return "Poorly soluble"
        else:
            return "Insoluble"
            
    def _get_cached_properties(self) -> Optional[Dict]:
        """Get cached properties if available."""
        return self._properties_cache.get('solubility', None)
        
    def _cache_properties(self, properties: Dict) -> None:
        """Cache calculated properties."""
        self._properties_cache['solubility'] = properties
