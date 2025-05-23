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
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors, Crippen
from typing import Dict, Optional
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
    
    def __init__(self, smiles: Optional[str] = None):
        """
        Initialize predictor. SMILES string can be provided at instantiation or later
        to the `calculate_admet_properties` method.
        
        Args:
            smiles (Optional[str]): SMILES representation of the molecule. Defaults to None.
        
        Raises:
            ValueError: If an initial SMILES string is provided and is invalid.
        """
        self.smiles: Optional[str] = None
        self.mol: Optional[Chem.Mol] = None
        self.properties: Dict[str, any] = {} # Initialize properties dictionary

        if smiles is not None:
            self._set_molecule_and_base_properties(smiles)
            
    def _set_molecule_and_base_properties(self, smiles_string: str):
        """
        Sets the molecule from a SMILES string and calculates RDKit-based physicochemical properties.

        Args:
            smiles_string (str): The SMILES string for the molecule.

        Raises:
            ValueError: If the SMILES string is invalid.
        """
        if not smiles_string:
            raise ValueError("SMILES string cannot be empty.")
            
        current_mol = Chem.MolFromSmiles(smiles_string)
        if not current_mol:
            raise ValueError(f"Invalid SMILES string provided: {smiles_string}")
        
        self.smiles = smiles_string
        self.mol = current_mol
        self.properties = {} # Reset properties for the new molecule
        self._calculate_molecular_properties() # Populate RDKit descriptors

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
        Also includes properties previously in PhysiochemicalProperties class.
        """
        try:
            # Basic properties that don't require 3D structure
            mass = Descriptors.ExactMolWt(self.mol)
            mol_formula = rdMolDescriptors.CalcMolFormula(self.mol)
            num_heavy_atoms_val = self.mol.GetNumHeavyAtoms()
            num_aromatic_heavy_atoms_val = len([atom for atom in self.mol.GetAtoms() if atom.GetIsAromatic()])
            fraction_csp3_val = round(Descriptors.FractionCSP3(self.mol), 2) if self.mol.GetNumHeavyAtoms() > 0 else 0.0 # Avoid error on empty mol
            molar_refractivity_val = round(Crippen.MolMR(self.mol), 2) if self.mol.GetNumHeavyAtoms() > 0 else 0.0 # Avoid error on empty mol

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
                'molecular_formula': mol_formula,           # Molecular formula string
                'num_heavy_atoms': num_heavy_atoms_val,     # [count] Number of heavy atoms
                'num_aromatic_heavy_atoms': num_aromatic_heavy_atoms_val, # [count] Number of aromatic heavy atoms
                'fraction_csp3': fraction_csp3_val,         # [0-1] Fraction of sp3 hybridized carbons
                'molar_refractivity': molar_refractivity_val, # [Å³] Molar refractivity
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
                'logP': Descriptors.MolLogP(self.mol),       # [log] Octanol-water partition coefficient
                'pKa_acid': pka_acid,                        # Acid dissociation constant
                'pKa_base': pka_base,                        # Base dissociation constant
            }
            
        except Exception as e:
            raise ValueError(f"Error calculating molecular properties: {str(e)}")
            
    def calculate_admet_properties(self, smiles: Optional[str] = None) -> Dict[str, float]:
        """
        Calculate comprehensive ADMET properties using deep learning models.
        SMILES string can be provided here if not set during initialization.

        Args:
            smiles (Optional[str]): SMILES representation of the molecule. 
                                    If provided, it overrides any SMILES set at initialization
                                    for this calculation run. Defaults to None.
        
        Returns:
            Dict[str, float]: Dictionary containing all predicted ADMET properties with units.
        
        Raises:
            ValueError: If no SMILES string is available (neither at init nor as an argument)
                        or if the provided SMILES is invalid.
        """
        target_smiles_str = smiles if smiles is not None else self.smiles

        if target_smiles_str is None:
            raise ValueError("SMILES string must be provided either during initialization or to this method.")

        # If a new SMILES is provided to this method, or if self.mol is not set (e.g. SMILES wasn't provided at init)
        # then (re)initialize the molecule and its base RDKit properties.
        if smiles is not None or self.mol is None:
            self._set_molecule_and_base_properties(target_smiles_str)
        
        # At this point, self.smiles, self.mol, and self.properties (with RDKit descriptors) are set.
        # Ensure self.smiles is used for ADMET-AI prediction, as it's confirmed valid by _set_molecule_and_base_properties
        pred = model.predict(self.smiles) # Use the validated self.smiles
        
        # These methods rely on self.properties which are populated by _calculate_molecular_properties
        mdck_permeability = self.predict_mdck_permeability()
        pgp_inhibitor_value = self.predict_pgp_inhibition()
        
        # Prepare values for plasma_clearance calculation
        # These also rely on self.properties for fallbacks if ADMET-AI doesn't provide them
        vdss_for_clearance = pred.get('VDss_Lombardo', self.properties.get('vdss', 0.0))
        hepatic_cl_for_clearance = pred.get('Clearance_Hepatocyte_AZ', self.properties.get('hepatic_clearance_val', 0.0)) # Assuming 'hepatic_clearance_val' if stored from another source
        ppbr_az_for_clearance = pred.get('PPBR_AZ', self.properties.get('ppb_val', 0.0)) # Assuming 'ppb_val' if stored from another source
        
        # Ensure ppbr_az_for_clearance is a number before subtraction
        if not isinstance(ppbr_az_for_clearance, (int, float)):
            ppbr_az_for_clearance = 0.0 # Default to 0 if not a number
            
        fu_for_clearance = 100.0 - ppbr_az_for_clearance
        
        plasma_clearance_value = self.predict_plasma_clearance(
            vdss=vdss_for_clearance,
            hepatic_cl=hepatic_cl_for_clearance,
            fu=fu_for_clearance
        )

        admet_properties = {
            # Physicochemical Property Percentiles (for RDKit calculated values)
            'molecular_weight_drugbank_approved_percentile': pred.get('molecular_weight_drugbank_approved_percentile'), # [%] Molecular Weight - DrugBank Approved Percentile
            'logP_drugbank_approved_percentile': pred.get('logP_drugbank_approved_percentile'),                         # [%] logP - DrugBank Approved Percentile
            'hydrogen_bond_acceptors_drugbank_approved_percentile': pred.get('hydrogen_bond_acceptors_drugbank_approved_percentile'), # [%] H-Bond Acceptors - DrugBank Approved Percentile
            'hydrogen_bond_donors_drugbank_approved_percentile': pred.get('hydrogen_bond_donors_drugbank_approved_percentile'),    # [%] H-Bond Donors - DrugBank Approved Percentile
            'tpsa_drugbank_approved_percentile': pred.get('tpsa_drugbank_approved_percentile'),                         # [%] Topological Polar Surface Area (TPSA) - DrugBank Approved Percentile

            # Drug-likeness Properties
            'qed': pred.get('QED'),                                         # [0-1] Quantitative Estimate of Drug-likeness
            'qed_drugbank_approved_percentile': pred.get('QED_drugbank_approved_percentile'), # [%] QED - DrugBank Approved Percentile
            'stereogenic_centers': pred.get('stereo_centers'),              # [count] Number of chiral centers
            'stereogenic_centers_drugbank_approved_percentile': pred.get('stereo_centers_drugbank_approved_percentile'), # [%] Stereogenic Centers - DrugBank Approved Percentile
            'lipinski_drugbank_approved_percentile': pred.get('Lipinski_drugbank_approved_percentile'), # [%] Lipinski Rule Compliance - DrugBank Approved Percentile
            
            # Absorption Properties
            'bbb': pred.get('BBB_Martins'),                                 # [0-1] Probability of BBB penetration
            'bbb_drugbank_approved_percentile': pred.get('BBB_Martins_drugbank_approved_percentile'), # [%] BBB Martins - DrugBank Approved Percentile
            'oral_bioavailability': pred.get('Bioavailability_Ma'),         # [0-1] Probability of >20% oral bioavailability
            'oral_bioavailability_drugbank_approved_percentile': pred.get('Bioavailability_Ma_drugbank_approved_percentile'), # [%] Oral Bioavailability (Ma) - DrugBank Approved Percentile
            'hia': pred.get('HIA_Hou'),                                     # [0-1] Fraction absorbed from GI tract
            'hia_drugbank_approved_percentile': pred.get('HIA_Hou_drugbank_approved_percentile'), # [%] HIA (Hou) - DrugBank Approved Percentile
            'pampa_permeability': pred.get('PAMPA_NCATS'),                  # [0-1] Artificial membrane permeability
            'pampa_permeability_drugbank_approved_percentile': pred.get('PAMPA_NCATS_drugbank_approved_percentile'), # [%] PAMPA (NCATS) - DrugBank Approved Percentile
            'caco2_permeability': pred.get('Caco2_Wang'),                   # log(10-⁶ cm/s] Caco-2 cell permeability
            'caco2_permeability_drugbank_approved_percentile': pred.get('Caco2_Wang_drugbank_approved_percentile'), # [%] Caco2 Wang - DrugBank Approved Percentile
            'mdck_permeability': mdck_permeability,                         # [log cm/s] MDCK cell permeability
            'pgp_substrate': pred.get('Pgp_Broccatelli'),                   # [0-1] Probability of P-gp substrate
            'pgp_substrate_drugbank_approved_percentile': pred.get('Pgp_Broccatelli_drugbank_approved_percentile'), # [%] Pgp Substrate (Broccatelli) - DrugBank Approved Percentile
            'pgp_inhibitor': pgp_inhibitor_value,                           # [0-1] Probability of P-gp inhibition
            
            # Distribution Properties
            'vdss': pred.get('VDss_Lombardo'),                              # [L/kg] Steady-state volume of distribution
            'vdss_drugbank_approved_percentile': pred.get('VDss_Lombardo_drugbank_approved_percentile'), # [%] VDss (Lombardo) - DrugBank Approved Percentile
            'ppb': min(100, max(0, pred.get('PPBR_AZ', 0))),               # [%] Percent bound to plasma proteins (capped at 0-100%)
            'ppb_drugbank_approved_percentile': pred.get('PPBR_AZ_drugbank_approved_percentile'), # [%] PPBR_AZ (underlying PPB) - DrugBank Approved Percentile
            'fraction_unbound': min(100, max(0, 100 - pred.get('PPBR_AZ', 0))), # [%] Free fraction in plasma (capped at 0-100%)
            'aqueous_solubility': pred.get('Solubility_AqSolDB'),           # [log mol/L] Water solubility
            'aqueous_solubility_drugbank_approved_percentile': pred.get('Solubility_AqSolDB_drugbank_approved_percentile'), # [%] AqSolDB Solubility - DrugBank Approved Percentile
            'lipophilicity': pred.get('Lipophilicity_AstraZeneca'),         # [log] Lipophilicity measure
            'lipophilicity_drugbank_approved_percentile': pred.get('Lipophilicity_AstraZeneca_drugbank_approved_percentile'), # [%] Lipophilicity (AstraZeneca) - DrugBank Approved Percentile
            'hydration_free_energy': pred.get('HydrationFreeEnergy_FreeSolv'), # [kcal/mol] Solvation energy
            'hydration_free_energy_drugbank_approved_percentile': pred.get('HydrationFreeEnergy_FreeSolv_drugbank_approved_percentile'), # [%] Hydration Free Energy (FreeSolv) - DrugBank Approved Percentile
            
            # Clearance and Half-life
            'hepatic_clearance': pred.get('Clearance_Hepatocyte_AZ'),       # [μL/min/mg] Hepatocyte clearance
            'hepatic_clearance_drugbank_approved_percentile': pred.get('Clearance_Hepatocyte_AZ_drugbank_approved_percentile'), # [%] Hepatic Clearance (AZ) - DrugBank Approved Percentile
            'microsomal_clearance': pred.get('Clearance_Microsome_AZ'),     # [μL/min/mg] Microsomal clearance
            'microsomal_clearance_drugbank_approved_percentile': pred.get('Clearance_Microsome_AZ_drugbank_approved_percentile'), # [%] Microsomal Clearance (AZ) - DrugBank Approved Percentile
            'plasma_clearance': plasma_clearance_value,                     # [mL/min/kg] Total plasma clearance
            'half_life': pred.get('Half_Life_Obach'),                       # [h] Elimination half-life
            'half_life_drugbank_approved_percentile': pred.get('Half_Life_Obach_drugbank_approved_percentile'), # [%] Half-Life (Obach) - DrugBank Approved Percentile
            
            # CYP450 Interactions
            'cyp1a2_inhibition': pred.get('CYP1A2_Veith'),                  # [0-1] CYP1A2 inhibition probability
            'cyp1a2_inhibition_drugbank_approved_percentile': pred.get('CYP1A2_Veith_drugbank_approved_percentile'), # [%] CYP1A2 Inhibition (Veith) - DrugBank Approved Percentile
            'cyp2c19_inhibition': pred.get('CYP2C19_Veith'),                # [0-1] CYP2C19 inhibition probability
            'cyp2c19_inhibition_drugbank_approved_percentile': pred.get('CYP2C19_Veith_drugbank_approved_percentile'), # [%] CYP2C19 Inhibition (Veith) - DrugBank Approved Percentile
            'cyp2c9_inhibition': pred.get('CYP2C9_Veith'),                  # [0-1] CYP2C9 inhibition probability
            'cyp2c9_inhibition_drugbank_approved_percentile': pred.get('CYP2C9_Veith_drugbank_approved_percentile'), # [%] CYP2C9 Inhibition (Veith) - DrugBank Approved Percentile
            'cyp2d6_inhibition': pred.get('CYP2D6_Veith'),                  # [0-1] CYP2D6 inhibition probability
            'cyp2d6_inhibition_drugbank_approved_percentile': pred.get('CYP2D6_Veith_drugbank_approved_percentile'), # [%] CYP2D6 Inhibition (Veith) - DrugBank Approved Percentile
            'cyp3a4_inhibition': pred.get('CYP3A4_Veith'),                  # [0-1] CYP3A4 inhibition probability
            'cyp3a4_inhibition_drugbank_approved_percentile': pred.get('CYP3A4_Veith_drugbank_approved_percentile'), # [%] CYP3A4 Inhibition (Veith) - DrugBank Approved Percentile
            'cyp2c9_substrate': pred.get('CYP2C9_Substrate_CarbonMangels'), # [0-1] CYP2C9 substrate probability
            'cyp2c9_substrate_drugbank_approved_percentile': pred.get('CYP2C9_Substrate_CarbonMangels_drugbank_approved_percentile'), # [%] CYP2C9 Substrate (CarbonMangels) - DrugBank Approved Percentile
            'cyp2d6_substrate': pred.get('CYP2D6_Substrate_CarbonMangels'), # [0-1] CYP2D6 substrate probability
            'cyp2d6_substrate_drugbank_approved_percentile': pred.get('CYP2D6_Substrate_CarbonMangels_drugbank_approved_percentile'), # [%] CYP2D6 Substrate (CarbonMangels) - DrugBank Approved Percentile
            'cyp3a4_substrate': pred.get('CYP3A4_Substrate_CarbonMangels'), # [0-1] CYP3A4 substrate probability
            'cyp3a4_substrate_drugbank_approved_percentile': pred.get('CYP3A4_Substrate_CarbonMangels_drugbank_approved_percentile'), # [%] CYP3A4 Substrate (CarbonMangels) - DrugBank Approved Percentile
            
            # Toxicity Properties
            'ames_mutagenicity': pred.get('AMES'),                          # [0-1] Bacterial mutagenicity
            'ames_mutagenicity_drugbank_approved_percentile': pred.get('AMES_drugbank_approved_percentile'), # [%] AMES Mutagenicity - DrugBank Approved Percentile
            'clinical_toxicity': pred.get('ClinTox'),                       # [0-1] Clinical trial toxicity risk
            'clinical_toxicity_drugbank_approved_percentile': pred.get('ClinTox_drugbank_approved_percentile'), # [%] Clinical Toxicity (ClinTox) - DrugBank Approved Percentile
            'drug_induced_liver_injury': pred.get('DILI'),                  # [0-1] Hepatotoxicity risk
            'drug_induced_liver_injury_drugbank_approved_percentile': pred.get('DILI_drugbank_approved_percentile'), # [%] DILI - DrugBank Approved Percentile
            'herg_inhibition': pred.get('hERG'),                            # [0-1] hERG channel blockade risk
            'herg_inhibition_drugbank_approved_percentile': pred.get('hERG_drugbank_approved_percentile'), # [%] hERG Inhibition - DrugBank Approved Percentile
            'lethal_dose_50': pred.get('LD50_Zhu'),                         # [log mol/kg] Median lethal dose
            'lethal_dose_50_drugbank_approved_percentile': pred.get('LD50_Zhu_drugbank_approved_percentile'), # [%] LD50 (Zhu) - DrugBank Approved Percentile
            'skin_sensitization': pred.get('Skin_Reaction'),                # [0-1] Skin reaction probability
            'skin_sensitization_drugbank_approved_percentile': pred.get('Skin_Reaction_drugbank_approved_percentile'), # [%] Skin Reaction - DrugBank Approved Percentile
            'carcinogens_lagunin_drugbank_approved_percentile': pred.get('Carcinogens_Lagunin_drugbank_approved_percentile'), # [%] Carcinogenicity (Lagunin) - DrugBank Approved Percentile

            # Nuclear Receptor Interactions
            'androgen_receptor_binding': pred.get('NR-AR'),                 # [0-1] AR binding probability
            'androgen_receptor_binding_drugbank_approved_percentile': pred.get('NR-AR_drugbank_approved_percentile'), # [%] NR-AR Binding - DrugBank Approved Percentile
            'androgen_receptor_lbd': pred.get('NR-AR-LBD'),                 # [0-1] AR-LBD binding probability
            'androgen_receptor_lbd_drugbank_approved_percentile': pred.get('NR-AR-LBD_drugbank_approved_percentile'), # [%] NR-AR-LBD Binding - DrugBank Approved Percentile
            'aryl_hydrocarbon_receptor': pred.get('NR-AhR'),                # [0-1] AhR activation probability
            'aryl_hydrocarbon_receptor_drugbank_approved_percentile': pred.get('NR-AhR_drugbank_approved_percentile'), # [%] NR-AhR Activation - DrugBank Approved Percentile
            'aromatase_inhibition': pred.get('NR-Aromatase'),               # [0-1] Aromatase inhibition
            'aromatase_inhibition_drugbank_approved_percentile': pred.get('NR-Aromatase_drugbank_approved_percentile'), # [%] NR-Aromatase Inhibition - DrugBank Approved Percentile
            'estrogen_receptor_binding': pred.get('NR-ER'),                 # [0-1] ER binding probability
            'estrogen_receptor_binding_drugbank_approved_percentile': pred.get('NR-ER_drugbank_approved_percentile'), # [%] NR-ER Binding - DrugBank Approved Percentile
            'estrogen_receptor_lbd': pred.get('NR-ER-LBD'),                 # [0-1] ER-LBD binding probability
            'estrogen_receptor_lbd_drugbank_approved_percentile': pred.get('NR-ER-LBD_drugbank_approved_percentile'), # [%] NR-ER-LBD Binding - DrugBank Approved Percentile
            'ppar_gamma': pred.get('NR-PPAR-gamma'),                        # [0-1] PPAR-γ activation probability
            'ppar_gamma_drugbank_approved_percentile': pred.get('NR-PPAR-gamma_drugbank_approved_percentile'), # [%] NR-PPAR-gamma Activation - DrugBank Approved Percentile
            
            # Stress Response Pathways
            'oxidative_stress_response': pred.get('SR-ARE'),                # [0-1] ARE pathway activation
            'oxidative_stress_response_drugbank_approved_percentile': pred.get('SR-ARE_drugbank_approved_percentile'), # [%] SR-ARE Activation - DrugBank Approved Percentile
            'dna_damage_response': pred.get('SR-ATAD5'),                    # [0-1] DNA damage response
            'dna_damage_response_drugbank_approved_percentile': pred.get('SR-ATAD5_drugbank_approved_percentile'), # [%] SR-ATAD5 Activation - DrugBank Approved Percentile
            'heat_shock_response': pred.get('SR-HSE'),                      # [0-1] Heat shock response
            'heat_shock_response_drugbank_approved_percentile': pred.get('SR-HSE_drugbank_approved_percentile'), # [%] SR-HSE Activation - DrugBank Approved Percentile
            'mitochondrial_toxicity': pred.get('SR-MMP'),                   # [0-1] Mitochondrial membrane potential
            'mitochondrial_toxicity_drugbank_approved_percentile': pred.get('SR-MMP_drugbank_approved_percentile'), # [%] SR-MMP - DrugBank Approved Percentile
            'p53_pathway_activation': pred.get('SR-p53'),                    # [0-1] p53 pathway activation
            'p53_pathway_activation_drugbank_approved_percentile': pred.get('SR-p53_drugbank_approved_percentile') # [%] SR-p53 Activation - DrugBank Approved Percentile
        }
        
        # Ensure all keys from pred are accessed with .get() if they might be missing
        # For keys that are expected to be present (like 'QED', 'BBB_Martins', etc, based on admet_ai model's core output)
        # if they were missing, it might indicate an issue with the model or input SMILES.
        # Using .get() for all 'pred' accesses makes it more robust to variations in 'pred' dictionary content.

        # Remove entries where the value is None (i.e., percentile not available from pred)
        admet_properties_cleaned = {k: v for k, v in admet_properties.items() if v is not None}

        self.properties.update(admet_properties_cleaned)
        # 
        # print(f"hepatic clearnce {self.properties['hepatic_clearance']}")
        return self.properties
        
    def get_physicochemical_properties(self) -> Dict[str, float]:
        """Returns physicochemical properties and their DrugBank percentiles.

        Ensure `calculate_admet_properties()` has been called before this method,
        as percentiles are sourced from its full computation. If not, only base
        RDKit properties calculated during __init__ will be returned.

        Returns:
            Dict[str, float]: A dictionary of physicochemical properties.
        """
        # Keys for base physicochemical properties calculated in _calculate_molecular_properties
        base_physicochem_keys = [
            'molecular_weight', 'molecular_formula', 'num_heavy_atoms',
            'num_aromatic_heavy_atoms', 'fraction_csp3', 'molar_refractivity',
            'molecular_volume', 'molecular_density',
            'hydrogen_bond_acceptors', 'hydrogen_bond_donors', 'rotatable_bonds',
            'ring_count', 'max_ring_size', 'heteroatom_count', 'formal_charge',
            'molecular_rigidity', 'molecular_flexibility',
            'topological_polar_surface_area', 'logP', 'pKa_acid', 'pKa_base'
        ]
        # Keys for their DrugBank percentiles (added by calculate_admet_properties)
        percentile_keys = [
            'molecular_weight_drugbank_approved_percentile',
            'logP_drugbank_approved_percentile',
            'hydrogen_bond_acceptors_drugbank_approved_percentile',
            'hydrogen_bond_donors_drugbank_approved_percentile',
            'tpsa_drugbank_approved_percentile'
        ]

        if not self.properties or 'qed' not in self.properties: # 'qed' is a proxy for full ADMET calculation
            # Only return base properties if full ADMET calculation hasn't run
            return {key: self.properties.get(key) for key in base_physicochem_keys if self.properties.get(key) is not None}
        else:
            # Return base properties plus available percentiles
            all_phys_keys = base_physicochem_keys + percentile_keys
            return {key: self.properties.get(key) for key in all_phys_keys if self.properties.get(key) is not None}

    def get_drug_likeness_properties(self) -> Dict[str, float]:
        """Returns drug-likeness properties.

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of drug-likeness properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        drug_likeness_keys = [
            'qed', 
            'qed_drugbank_approved_percentile',
            'stereogenic_centers', 
            'stereogenic_centers_drugbank_approved_percentile',
            'lipinski_drugbank_approved_percentile'
        ]
        return {key: self.properties[key] for key in drug_likeness_keys if key in self.properties}

    def get_absorption_properties(self) -> Dict[str, float]:
        """Returns absorption properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of absorption properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        absorption_keys = [
            'bbb',
            'bbb_drugbank_approved_percentile',
            'oral_bioavailability',
            'oral_bioavailability_drugbank_approved_percentile',
            'hia',
            'hia_drugbank_approved_percentile',
            'pampa_permeability',
            'pampa_permeability_drugbank_approved_percentile',
            'caco2_permeability',
            'caco2_permeability_drugbank_approved_percentile',
            'mdck_permeability',
            'pgp_substrate', 'pgp_substrate_drugbank_approved_percentile',
            'pgp_inhibitor',
            
        ]
        return {key: self.properties[key] for key in absorption_keys if key in self.properties}

    def get_distribution_properties(self) -> Dict[str, float]:
        """Returns distribution properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of distribution properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        distribution_keys = [
            'vdss', 
            'vdss_drugbank_approved_percentile',
            'ppb',
            'ppb_drugbank_approved_percentile',
            'fraction_unbound',
            'aqueous_solubility',
            'aqueous_solubility_drugbank_approved_percentile',
            'lipophilicity',
            'lipophilicity_drugbank_approved_percentile',
            'hydration_free_energy', 
            'hydration_free_energy_drugbank_approved_percentile',
        ]
        return {key: self.properties[key] for key in distribution_keys if key in self.properties}

    def get_metabolism_properties(self) -> Dict[str, float]:
        """Returns metabolism/CYP450 interaction properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of metabolism properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        metabolism_keys = [
            'cyp1a2_inhibition', 
            'cyp1a2_inhibition_drugbank_approved_percentile',
            'cyp2c19_inhibition',
            'cyp2c19_inhibition_drugbank_approved_percentile',
            'cyp2c9_inhibition',
            'cyp2c9_inhibition_drugbank_approved_percentile',
            'cyp2d6_inhibition',
            'cyp2d6_inhibition_drugbank_approved_percentile',
            'cyp3a4_inhibition',
            'cyp3a4_inhibition_drugbank_approved_percentile',
            'cyp2c9_substrate',
            'cyp2c9_substrate_drugbank_approved_percentile',
            'cyp2d6_substrate',
            'cyp2d6_substrate_drugbank_approved_percentile',
            'cyp3a4_substrate',
            'cyp3a4_substrate_drugbank_approved_percentile'
        ]
        return {key: self.properties[key] for key in metabolism_keys if key in self.properties}

    def get_excretion_properties(self) -> Dict[str, float]:
        """Returns excretion (clearance and half-life) properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of excretion properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        excretion_keys = [
            'hepatic_clearance',
            'hepatic_clearance_drugbank_approved_percentile',
            'microsomal_clearance',
            'microsomal_clearance_drugbank_approved_percentile',
            'plasma_clearance',
            'half_life',
            'half_life_drugbank_approved_percentile'
        ]
        return {key: self.properties[key] for key in excretion_keys if key in self.properties}

    def get_toxicity_properties(self) -> Dict[str, float]:
        """Returns toxicity properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of toxicity properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        toxicity_keys = [
            'ames_mutagenicity',
            'ames_mutagenicity_drugbank_approved_percentile',
            'clinical_toxicity',
            'clinical_toxicity_drugbank_approved_percentile',
            'drug_induced_liver_injury',
            'drug_induced_liver_injury_drugbank_approved_percentile',
            'herg_inhibition',
            'herg_inhibition_drugbank_approved_percentile',
            'lethal_dose_50',
            'lethal_dose_50_drugbank_approved_percentile',
            'skin_sensitization',
            'skin_sensitization_drugbank_approved_percentile',
            'carcinogens_lagunin_drugbank_approved_percentile'
        ]
        return {key: self.properties[key] for key in toxicity_keys if key in self.properties}

    def get_nuclear_receptor_properties(self) -> Dict[str, float]:
        """Returns nuclear receptor interaction properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of nuclear receptor properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        nuclear_receptor_keys = [
            'androgen_receptor_binding', 
            'androgen_receptor_binding_drugbank_approved_percentile',
            'androgen_receptor_lbd',
            'androgen_receptor_lbd_drugbank_approved_percentile',
            'aryl_hydrocarbon_receptor',
            'aryl_hydrocarbon_receptor_drugbank_approved_percentile',
            'aromatase_inhibition',
            'aromatase_inhibition_drugbank_approved_percentile',
            'estrogen_receptor_binding',
            'estrogen_receptor_binding_drugbank_approved_percentile',
            'estrogen_receptor_lbd',
            'estrogen_receptor_lbd_drugbank_approved_percentile',
            'ppar_gamma',
            'ppar_gamma_drugbank_approved_percentile'
        ]
        return {key: self.properties[key] for key in nuclear_receptor_keys if key in self.properties}

    def get_stress_response_properties(self) -> Dict[str, float]:
        """Returns stress response pathway properties (ADMET).

        Ensure `calculate_admet_properties()` has been called before this method.

        Returns:
            Dict[str, float]: A dictionary of stress response properties.
        """
        if not self.properties or 'qed' not in self.properties:
            return {}
        stress_response_keys = [
            'oxidative_stress_response',
            'oxidative_stress_response_drugbank_approved_percentile',
            'dna_damage_response',
            'dna_damage_response_drugbank_approved_percentile',
            'heat_shock_response',
            'heat_shock_response_drugbank_approved_percentile',
            'mitochondrial_toxicity',
            'mitochondrial_toxicity_drugbank_approved_percentile',
            'p53_pathway_activation',
            'p53_pathway_activation_drugbank_approved_percentile'
        ]
        return {key: self.properties[key] for key in stress_response_keys if key in self.properties}

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
