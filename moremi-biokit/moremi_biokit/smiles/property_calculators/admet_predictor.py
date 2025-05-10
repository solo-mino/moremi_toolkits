"""
ADMET Property Predictor

This module implements a comprehensive ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) 
prediction system using deep learning models (ADMET_AI) and molecular property calculations (RDKit).

Key Features:
1. Molecular property calculations using RDKit
2. ADMET predictions using state-of-the-art deep learning models
3. Comprehensive property coverage including:
   - Physicochemical properties
   - ADME properties
   - Toxicity profiles
   - Nuclear receptor interactions
   - Stress response characteristics
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from typing import Dict
import math
from openbabel import openbabel as ob
from admet_ai import ADMETModel


# Initialize the ADMET Models
model = ADMETModel()


class ADMETPredictor:
    """
    Comprehensive ADMET property predictor using deep learning models and molecular descriptors.
    Provides predictions for absorption, distribution, metabolism, excretion, and toxicity properties.
    """
    
    def __init__(self, smiles: str):
        """
        Initialize predictor with SMILES string.
        
        Args:
            smiles (str): SMILES representation of the molecule
        
        Raises:
            ValueError: If SMILES string is invalid
        """
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(self.smiles)
        if not self.mol:
            raise ValueError("Invalid SMILES string")
            
        # Calculate molecular properties
        self.properties = {}
        self._calculate_molecular_properties()
        
    def calculate_pka_openbabel(self):
        """
        Calculate pKa values using OpenBabel.
        
        Returns:
            tuple: (acidic_pka, basic_pka) where each may be None if not found
        """
        # Convert RDKit mol to OpenBabel
        obmol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("smi", "smi")
        conv.ReadString(obmol, Chem.MolToSmiles(self.mol))
        
        # Calculate pKa
        charges = ob.OBChargeModel.FindType("gasteiger")
        charges.ComputeCharges(obmol)
        
        acidic_pkas = []
        basic_pkas = []
        
        for atom in ob.OBMolAtomIter(obmol):
            if atom.GetFormalCharge() < 0:
                acidic_pkas.append(atom.GetData("pKa").GetDoubleValue())
            elif atom.GetFormalCharge() > 0:
                basic_pkas.append(atom.GetData("pKa").GetDoubleValue())
                
        return (
            min(acidic_pkas) if acidic_pkas else None,
            max(basic_pkas) if basic_pkas else None
        )
        
    def _calculate_molecular_properties(self):
        """
        Calculate fundamental molecular properties using RDKit.
        Includes 3D structure generation, volume, density, ring analysis, and flexibility metrics.
        """
        try:
            # Basic properties that don't require 3D structure
            mass = Descriptors.ExactMolWt(self.mol)
            
            # Initialize with default values for 3D-dependent properties
            volume = None
            density = None
            
            # Try 3D Structure Generation with error handling
            mol_3d = Chem.AddHs(self.mol)
            
            # Use more robust embedding parameters
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = True
            params.maxIterations = 2000  # Increase iterations for better convergence
            
            # Attempt embedding and handle failure gracefully
            conf_id = AllChem.EmbedMolecule(mol_3d, params)
            
            if conf_id >= 0:  # Successful embedding
                try:
                    # Try force field optimization
                    AllChem.MMFFOptimizeMolecule(mol_3d, confId=0, maxIters=2000)
                    
                    # Calculate volume after successful 3D generation
                    volume = AllChem.ComputeMolVolume(mol_3d)
                    density = mass / volume if volume and volume > 0 else None
                except Exception as e:
                    # If optimization fails, log it but continue with other calculations
                    print(f"⚠️ Warning: 3D optimization failed: {str(e)}")
            else:
                # If embedding fails, log it but continue
                print("⚠️ Warning: 3D structure generation failed")

            # Ring Analysis - doesn't require 3D structure
            ring_info = self.mol.GetRingInfo()
            max_ring_size = max([len(ring) for ring in ring_info.AtomRings()]) if ring_info.AtomRings() else 0
            
            # Rigidity and Flexibility - doesn't require 3D structure
            n_rotatable = Descriptors.NumRotatableBonds(self.mol)
            total_bonds = self.mol.GetNumBonds()
            flexibility = n_rotatable / total_bonds if total_bonds > 0 else 0
            rigidity = 1 - flexibility
            
            # Add pKa values
            try:
                pka_acid, pka_base = self.calculate_pka_openbabel()
            except Exception as e:
                print(f"⚠️ Warning: pKa calculation failed: {str(e)}")
                pka_acid, pka_base = None, None
            
            self.properties = {
                'molecular_weight': mass,                    # [g/mol] Exact molecular weight
                'molecular_volume': volume,                  # [Å³] Van der Waals volume
                'molecular_density': density,                # [g/cm³] Molecular density
                'hydrogen_bond_acceptors': Descriptors.NumHAcceptors(self.mol),  # [count] Number of H-bond acceptors
                'hydrogen_bond_donors': Descriptors.NumHDonors(self.mol),        # [count] Number of H-bond donors
                'rotatable_bonds': n_rotatable,             # [count] Number of rotatable bonds
                'ring_count': Descriptors.RingCount(self.mol),  # [count] Total number of rings
                'max_ring_size': max_ring_size,             # [count] Size of largest ring
                'heteroatom_count': len([atom for atom in self.mol.GetAtoms()   # [count] Number of non-C/H atoms
                                if atom.GetSymbol() not in ['C', 'H']]),
                'formal_charge': Chem.GetFormalCharge(self.mol),  # [e] Net molecular charge
                'molecular_rigidity': rigidity,             # [0-1] Molecular rigidity score
                'molecular_flexibility': flexibility,        # [0-1] Molecular flexibility score
                'topological_polar_surface_area': Descriptors.TPSA(self.mol, includeSandP=False),  # [Å²] TPSA
                'logP': Descriptors.MolLogP(self.mol)       # [log] Octanol-water partition coefficient
            }
            
            # Add pKa values separately to avoid dict update issues
            if pka_acid is not None:
                self.properties['pKa_acid'] = pka_acid  # Acid dissociation constant
            if pka_base is not None:
                self.properties['pKa_base'] = pka_base  # Base dissociation constant
            
        except Exception as e:
            raise ValueError(f"Error calculating molecular properties: {str(e)}")
            
    def calculate_admet_properties(self) -> Dict[str, float]:
        """
        Calculate comprehensive ADMET properties using deep learning models.
        
        Returns:
            Dict[str, float]: Dictionary containing all predicted ADMET properties with units
        """
        pred = model.predict(self.smiles)
        mdck_permeability = self.predict_mdck_permeability()
        
        admet_properties = {
            # Drug-likeness Properties
            'qed': pred['QED'],                    # [0-1] Quantitative Estimate of Drug-likeness
            'stereogenic_centers': pred['stereo_centers'],         # [count] Number of chiral centers
            
            # Absorption Properties
            'bbb': pred['BBB_Martins'],      # [0-1] Probability of BBB penetration
            'oral_bioavailability': pred['Bioavailability_Ma'],          # [0-1] Probability of >20% oral bioavailability
            'hia': pred['HIA_Hou'],              # [0-1] Fraction absorbed from GI tract
            'pampa_permeability': pred['PAMPA_NCATS'],                  # [0-1] Artificial membrane permeability
            'caco2_permeability': pred['Caco2_Wang'],                   # log(10-⁶ cm/s] Caco-2 cell permeability
            'mdck_permeability': mdck_permeability,                     # [log cm/s] MDCK cell permeability
            'pgp_substrate': pred['Pgp_Broccatelli'],        # [0-1] Probability of P-gp substrate
            'pgp_inhibitor': self.predict_pgp_inhibition(),  # [0-1] Probability of P-gp inhibition
            
            # Distribution Properties
            'vdss': pred['VDss_Lombardo'],            # [L/kg] Steady-state volume of distribution
            'ppb': min(100, max(0, pred['PPBR_AZ'])),                  # [%] Percent bound to plasma proteins (capped at 0-100%)
            'fraction_unbound': min(100, max(0, 100 - pred['PPBR_AZ'])),     # [%] Free fraction in plasma (capped at 0-100%)
            'aqueous_solubility': pred['Solubility_AqSolDB'],          # [log mol/L] Water solubility
            'lipophilicity': pred['Lipophilicity_AstraZeneca'],        # [log] Lipophilicity measure
            'hydration_free_energy': pred['HydrationFreeEnergy_FreeSolv'],  # [kcal/mol] Solvation energy
            
            # Clearance and Half-life
            'hepatic_clearance': pred['Clearance_Hepatocyte_AZ'],      # [μL/min/mg] Hepatocyte clearance
            'microsomal_clearance': pred['Clearance_Microsome_AZ'],    # [μL/min/mg] Microsomal clearance
            'plasma_clearance': self.predict_plasma_clearance(
                vdss=pred['VDss_Lombardo'],
                hepatic_cl=pred['Clearance_Hepatocyte_AZ'],
                fu=(100 - pred['PPBR_AZ'])
            ),# [mL/min/kg] Total plasma clearance
            'half_life': pred['Half_Life_Obach'],                      # [h] Elimination half-life
            
            # CYP450 Interactions
            'cyp1a2_inhibition': pred['CYP1A2_Veith'],                # [0-1] CYP1A2 inhibition probability
            'cyp2c19_inhibition': pred['CYP2C19_Veith'],              # [0-1] CYP2C19 inhibition probability
            'cyp2c9_inhibition': pred['CYP2C9_Veith'],                # [0-1] CYP2C9 inhibition probability
            'cyp2d6_inhibition': pred['CYP2D6_Veith'],                # [0-1] CYP2D6 inhibition probability
            'cyp3a4_inhibition': pred['CYP3A4_Veith'],                # [0-1] CYP3A4 inhibition probability
            'cyp2c9_substrate': pred['CYP2C9_Substrate_CarbonMangels'],  # [0-1] CYP2C9 substrate probability
            'cyp2d6_substrate': pred['CYP2D6_Substrate_CarbonMangels'],  # [0-1] CYP2D6 substrate probability
            'cyp3a4_substrate': pred['CYP3A4_Substrate_CarbonMangels'],  # [0-1] CYP3A4 substrate probability
            
            # Toxicity Properties
            'ames_mutagenicity': pred['AMES'],                         # [0-1] Bacterial mutagenicity
            'clinical_toxicity': pred['ClinTox'],                      # [0-1] Clinical trial toxicity risk
            'drug_induced_liver_injury': pred['DILI'],                 # [0-1] Hepatotoxicity risk
            'herg_inhibition': pred['hERG'],                          # [0-1] hERG channel blockade risk
            'lethal_dose_50': pred['LD50_Zhu'],                       # [log mol/kg] Median lethal dose
            'skin_sensitization': pred['Skin_Reaction'],              # [0-1] Skin reaction probability
            
            # Nuclear Receptor Interactions
            'androgen_receptor_binding': pred['NR-AR'],               # [0-1] AR binding probability
            'androgen_receptor_lbd': pred['NR-AR-LBD'],              # [0-1] AR-LBD binding probability
            'aryl_hydrocarbon_receptor': pred['NR-AhR'],             # [0-1] AhR activation probability
            'aromatase_inhibition': pred['NR-Aromatase'],            # [0-1] Aromatase inhibition
            'estrogen_receptor_binding': pred['NR-ER'],              # [0-1] ER binding probability
            'estrogen_receptor_lbd': pred['NR-ER-LBD'],             # [0-1] ER-LBD binding probability
            'ppar_gamma': pred['NR-PPAR-gamma'],                     # [0-1] PPAR-γ activation probability
            
            # Stress Response Pathways
            'oxidative_stress_response': pred['SR-ARE'],             # [0-1] ARE pathway activation
            'dna_damage_response': pred['SR-ATAD5'],                 # [0-1] DNA damage response
            'heat_shock_response': pred['SR-HSE'],                   # [0-1] Heat shock response
            'mitochondrial_toxicity': pred['SR-MMP'],               # [0-1] Mitochondrial membrane potential
            'p53_pathway_activation': pred['SR-p53']                # [0-1] p53 pathway activation
        }
        
        self.properties.update(admet_properties)
        # 
        # print(f"hepatic clearnce {self.properties['hepatic_clearance']}")
        return self.properties
        
    def predict_mdck_permeability(self) -> float:
        """
        Predict MDCK cell permeability using molecular descriptors.
        
        Returns:
            float: Predicted log Papp [cm/s] for MDCK cell permeability
        """
        desc = self.properties
        
        # Base permeability value
        log_papp = -4.5  # Starting point based on literature
        
        # 1. Size and Shape Effects
        mw_factor = -0.004 * (desc['molecular_weight'] - 350)/100  # Optimal MW around 350
        log_papp += mw_factor
        
        # 2. Lipophilicity Contribution
        logp = desc.get('logP', 0)
        if 1 <= logp <= 3:  # Optimal range
            log_papp += 0.3
        elif logp > 3:
            log_papp += 0.3 - 0.1 * (logp - 3)  # Penalty for high lipophilicity
        
        # 3. Hydrogen Bonding
        hbd = desc.get('hydrogen_bond_donors', 0)
        hba = desc.get('hydrogen_bond_acceptors', 0)
        h_bond_factor = -0.1 * (hbd + hba/2)/5  # Weighted penalty
        log_papp += h_bond_factor
        
        # 4. Flexibility Impact
        flex = desc.get('molecular_flexibility', 0)
        if flex > 0.5:  # Penalty for high flexibility
            log_papp -= 0.2 * (flex - 0.5)
            
        return max(-6.0, min(-3.0, log_papp))  # Constrain to realistic range

    def predict_pgp_inhibition(self) -> float:
        """
        Predict P-glycoprotein inhibition probability using established pharmacophore model.
        Based on the following references:
        1. Wang et al. (2011) J Chem Inf Model. DOI: 10.1021/ci200254h
        2. Broccatelli et al. (2011) J Med Chem. DOI: 10.1021/jm101421d
        3. Chen et al. (2012) Mol Pharm. DOI: 10.1021/mp300202c
        
        The model considers:
        - Molecular weight (optimal range: 400-800 Da)
        - LogP (optimal range: 2.5-5)
        - Number of aromatic rings (optimal: 3-4)
        - Number of rotatable bonds (optimal: 6-12)
        - Number of H-bond acceptors (optimal: 2-4)
        - TPSA (optimal range: 70-150 Å²)
        
        Returns:
            float: Probability [0-1] of the molecule being a P-gp inhibitor
        """
        desc = self.properties
        
        # Initialize base score
        score = 0.5  # Start at middle probability
        
        # 1. Molecular Weight contribution (Wang et al., 2011)
        mw = desc['molecular_weight']
        if 400 <= mw <= 800:
            score += 0.1
        elif mw > 800:
            score -= 0.1
            
        # 2. LogP contribution (Broccatelli et al., 2011)
        logp = desc['logP']
        if 2.5 <= logp <= 5:
            score += 0.15
        elif logp > 5:
            score -= 0.1
            
        # 3. Aromatic Rings (Chen et al., 2012)
        aromatic_rings = Descriptors.NumAromaticRings(self.mol)
        if 3 <= aromatic_rings <= 4:
            score += 0.1
        elif aromatic_rings > 4:
            score -= 0.05
            
        # 4. Rotatable Bonds
        rot_bonds = desc['rotatable_bonds']
        if 6 <= rot_bonds <= 12:
            score += 0.05
        elif rot_bonds > 12:
            score -= 0.05
            
        # 5. H-bond Acceptors
        hba = desc['hydrogen_bond_acceptors']
        if 2 <= hba <= 4:
            score += 0.05
        elif hba > 4:
            score -= 0.05
            
        # 6. TPSA contribution
        tpsa = desc['topological_polar_surface_area']
        if 70 <= tpsa <= 150:
            score += 0.05
        elif tpsa > 150:
            score -= 0.05
            
        # Normalize score to [0-1] range
        return max(0, min(1, score))

    def predict_plasma_clearance(self, vdss: float, hepatic_cl: float, fu: float) -> float:
        """
        Predict total plasma clearance using physiological parameters and molecular properties.
        
        The model considers:
        - Hepatic clearance
        - Renal clearance (based on molecular properties)
        - Plasma protein binding
        - Volume of distribution
        
        Args:
            vdss (float): Volume of distribution [L/kg]
            hepatic_cl (float): Hepatic clearance [μL/min/mg]
            fu (float): Fraction unbound (as a percentage 0-100)
        
        Returns:
            float: Predicted total plasma clearance [mL/min/kg]
            
        References:
        1. Gibaldi & Perrier (1982) Pharmacokinetics, 2nd Ed. Marcel Dekker, NY
            - Basic principles of clearance calculation
            - Relationship between protein binding and clearance
        
        2. Smith et al. (2015) J Med Chem. DOI: 10.1021/jm501400z
            - Modern approaches to clearance prediction
            - Molecular descriptors influence on clearance
        
        3. Rowland & Tozer (2011) Clinical Pharmacokinetics and Pharmacodynamics
            - Volume of distribution influence on clearance
            - Protein binding corrections
        
        4. Lin & Lu (1997) Clin Pharmacokinet 33:184-209
            - Hepatic clearance models
            - Well-stirred model implementation
        """
        try:
            # Get required properties
            ppb = self.properties.get('ppb', 0)  # % bound to plasma proteins
            mw = self.properties.get('molecular_weight', 0)  # Da
            logp = self.properties.get('logP', 0)
            
            # Handle extreme values for fraction unbound
            if fu < 0:
                # If fu is negative (because PPBR_AZ > 100), set to a small positive value
                fraction_unbound = 0.01  # 1% unbound as minimum
            elif fu > 100:
                # If fu is > 100% (because PPBR_AZ < 0), cap at 100%
                fraction_unbound = 1.0
            else:
                # Normal case: convert percentage to fraction
                fraction_unbound = fu / 100
            
            # Handle extreme values for hepatic clearance
            if hepatic_cl < 0:
                hepatic_cl = 0  # Set negative clearance to zero
            
            # Fix ppb if it's out of range
            if ppb > 100:
                ppb = 100
            elif ppb < 0:
                ppb = 0
            
            # 1. Convert hepatic clearance to mL/min/kg using well-stirred model
            # Reference: Lin & Lu (1997)
            Q_h = 21  # Hepatic blood flow [mL/min/kg]
            
            # Calculate hepatic clearance contribution
            if hepatic_cl > 0:
                E_h = hepatic_cl / (Q_h + hepatic_cl)  # Hepatic extraction ratio
                hepatic_cl_total = Q_h * E_h * (fraction_unbound / (1 - E_h + fraction_unbound * E_h))
            else:
                hepatic_cl_total = 0  # If hepatic_cl is 0, set the total to 0
            
            # 2. Calculate renal clearance
            renal_cl = self._calculate_renal_clearance(mw, logp, fraction_unbound)
            
            # 3. Apply Volume of Distribution correction
            # Reference: Rowland & Tozer (2011)
            vd_factor = 1.0
            if vdss > 0:
                # Correction factor based on Vdss
                vd_factor = math.sqrt(vdss / 0.7)  # 0.7 L/kg is typical Vdss
                vd_factor = max(0.5, min(2.0, vd_factor))  # Limit correction
            
            # 4. Calculate total clearance with Vdss correction
            total_cl = (hepatic_cl_total + renal_cl) * vd_factor
            
            # 5. Apply protein binding adjustment
            # Reference: Gibaldi & Perrier (1982)
            if ppb > 0:
                fu_plasma = (100 - ppb) / 100
                # Ensure fu_plasma is reasonable
                fu_plasma = max(0.01, min(1.0, fu_plasma))
                
                # Avoid division by zero and ensure meaningful calculation
                total_cl *= (fraction_unbound / fu_plasma)
            
            # 6. Apply physiological constraints
            max_cl = 100  # Maximum physiological clearance ~100 mL/min/kg
            total_cl = min(max_cl, total_cl)
            
            return round(total_cl, 3)
            
        except Exception as e:
            raise ValueError(f"Error calculating plasma clearance: {str(e)}")
    
    def _calculate_renal_clearance(self, mw: float, logp: float, fu: float) -> float:
        """
        Calculate renal clearance based on molecular properties.
        Reference: Smith et al. (2015)
        
        Args:
            mw (float): Molecular weight [Da]
            logp (float): LogP value
            fu (float): Fraction unbound
            
        Returns:
            float: Renal clearance [mL/min/kg]
        """
        renal_cl = 0
        gfr = 120  # Standard glomerular filtration rate [mL/min]
        fu /= 100 # Convert to fraction
        
        # Glomerular filtration contribution
        if mw < 500:  # Small molecules can be filtered
            renal_cl += gfr * fu * 0.01  # Scale to per kg basis
        
        # Active secretion contribution
        if -1 <= logp <= 3:  # Optimal range for renal secretion
            secretion_factor = math.exp(-abs(logp - 1))  # Peak at logP = 1
            renal_cl += 30 * secretion_factor * fu * 0.01  # Max 30 mL/min/kg
        
        return renal_cl
