import subprocess
import os
from datetime import datetime
import re

def combine_pdb_files(protein1_path, protein2_path, output_path):
    """Combine two PDB files into one, with the second protein's chains renamed"""
    print(f"Combining PDB files into {output_path}")
    
    with open(protein1_path, 'r') as f1, open(protein2_path, 'r') as f2, open(output_path, 'w') as out:
        # Write first protein (chain A)
        for line in f1:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                # Ensure chain ID is A
                new_line = line[:21] + 'A' + line[22:]
                out.write(new_line)
            elif not line.startswith('END'):
                out.write(line)
        out.write('TER\n')
        
        # Write second protein (chain B)
        for line in f2:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                # Ensure chain ID is B
                new_line = line[:21] + 'B' + line[22:]
                out.write(new_line)
            elif not line.startswith('END'):
                out.write(line)
        out.write('END\n')

def save_results(sequence_num, affinity, prodigy_output, protein_path, antigen_path, error=None, output_dir="hepatitis_binding_results"):
    """Save results to a text file with sequence number and timestamp"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"binding_affinity_seq{sequence_num}_{timestamp}.txt")
    
    with open(output_file, 'w') as f:
        f.write("Binding Affinity Prediction Summary\n")
        f.write("=================================\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Sequence Number: {sequence_num}\n")
        f.write(f"Protein: {os.path.basename(protein_path)}\n")
        f.write(f"Antigen: {os.path.basename(antigen_path)}\n\n")
        
        if error:
            f.write(f"Error: {error}\n")
        else:
            f.write(f"Predicted Binding Affinity: {affinity:.2f} kcal/mol\n")
            f.write("Note: More negative values indicate stronger binding\n\n")
            f.write("Method Information:\n")
            f.write("------------------\n")
            f.write("Tool: PRODIGY (PROtein binDIng enerGY prediction)\n")
            f.write("Version: 2.2.5\n")
            f.write("Source: https://github.com/haddocking/prodigy\n")
            f.write("Reference: Xue LC, et al. PRODIGY: a web server for predicting the binding affinity\n")
            f.write("           of protein-protein complexes. Bioinformatics (2016) 32:3676–3678.\n")
            f.write("           DOI: 10.1093/bioinformatics/btw514\n\n")
            f.write("\nRaw PRODIGY Output:\n")
            f.write("------------------\n")
            cleaned_output = []
            for line in prodigy_output.split('\n'):
                if '[+] Reading structure file:' in line:
                    line = '[+] Reading structure file: complex.pdb'
                elif 'Parsed structure file' in line:
                    line = line.replace('/Users/aideveloperminohealth/Desktop/useless/things/', '')
                cleaned_output.append(line)
            f.write('\n'.join(cleaned_output))
    
    print(f"Results saved to: {output_file}")

def predict_binding_affinity(antibody_path, antigen_path, system_type=None):
    """
    Predict binding affinity between two protein structures using PRODIGY
    
    Args:
        antibody_path (str): Path to first PDB file
        antigen_path (str): Path to second PDB file
        system_type (str): Operating system type ('windows_wsl', 'macos', 'linux', None)
                          If None, will auto-detect the system
    
    Returns:
        dict: Dictionary containing binding affinity prediction results including:
            - intermolecular_contacts: Total number of contacts
            - charged_charged_contacts: Number of charged-charged contacts
            - charged_polar_contacts: Number of charged-polar contacts
            - charged_apolar_contacts: Number of charged-apolar contacts
            - polar_polar_contacts: Number of polar-polar contacts
            - apolar_polar_contacts: Number of apolar-polar contacts
            - apolar_apolar_contacts: Number of apolar-apolar contacts
            - apolar_nis_percentage: Percentage of apolar NIS residues
            - charged_nis_percentage: Percentage of charged NIS residues
            - binding_affinity: Predicted binding affinity in kcal/mol
            - dissociation_constant: Predicted dissociation constant in M at 25°C
    """
    print(f"Processing files:\n- {antibody_path}\n- {antigen_path}")
    
    # Auto-detect system type if not provided
    if system_type is None:
        import platform
        system = platform.system().lower()
        if system == 'windows':
            # Check if running in WSL
            if os.path.exists('/proc/sys/fs/binfmt_misc/WSLInterop'):
                system_type = 'linux'
            else:
                system_type = 'windows_wsl'
        elif system == 'darwin':
            system_type = 'macos'
        else:
            system_type = 'linux'
    
    print(f"Using system type: {system_type}")
    
    # Combine PDB files
    combined_pdb = "complex.pdb"
    combine_pdb_files(antibody_path, antigen_path, combined_pdb)
    
    try:
        # Run PRODIGY prediction
        if system_type == 'windows_wsl':
            # Convert Windows paths to WSL paths
            wsl_protein1 = subprocess.getoutput(f'wsl wslpath "{antibody_path}"')
            wsl_protein2 = subprocess.getoutput(f'wsl wslpath "{antigen_path}"')
            
            # Create a temporary directory in WSL
            temp_dir = subprocess.getoutput('wsl mktemp -d')
            
            # Copy files to WSL temp directory
            subprocess.run(f'wsl cp "{wsl_protein1}" "{temp_dir}/protein1.pdb"', shell=True)
            subprocess.run(f'wsl cp "{wsl_protein2}" "{temp_dir}/protein2.pdb"', shell=True)
            
            try:
                # Change to temp directory and run prodigy in WSL
                cmd = f'wsl cd "{temp_dir}" && prodigy complex.pdb'
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            finally:
                # Cleanup
                subprocess.run(f'wsl rm -rf "{temp_dir}"', shell=True)
        else:
            # For macOS and Linux, run prodigy directly
            cmd = ['prodigy', combined_pdb]
            print(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
        
        # print("\nPRODIGY Output:")
        # print(result.stdout)
        
        if result.stderr:
            print("\nErrors:")
            print(result.stderr)
            
        if result.returncode != 0:
            raise Exception(result.stderr or "Command failed with no error message")
            
        # Parse PRODIGY output to extract all values
        output_dict = {}
        
        # Extract intermolecular contacts
        contacts_match = re.search(r'No\. of intermolecular contacts: (\d+)', result.stdout)
        if contacts_match:
            output_dict['intermolecular_contacts'] = int(contacts_match.group(1))
        
        # Extract charged-charged contacts
        charged_charged_match = re.search(r'No\. of charged-charged contacts: (\d+)', result.stdout)
        if charged_charged_match:
            output_dict['charged_charged_contacts'] = int(charged_charged_match.group(1))
        
        # Extract charged-polar contacts
        charged_polar_match = re.search(r'No\. of charged-polar contacts: (\d+)', result.stdout)
        if charged_polar_match:
            output_dict['charged_polar_contacts'] = int(charged_polar_match.group(1))
        
        # Extract charged-apolar contacts
        charged_apolar_match = re.search(r'No\. of charged-apolar contacts: (\d+)', result.stdout)
        if charged_apolar_match:
            output_dict['charged_apolar_contacts'] = int(charged_apolar_match.group(1))
        
        # Extract polar-polar contacts
        polar_polar_match = re.search(r'No\. of polar-polar contacts: (\d+)', result.stdout)
        if polar_polar_match:
            output_dict['polar_polar_contacts'] = int(polar_polar_match.group(1))
        
        # Extract apolar-polar contacts
        apolar_polar_match = re.search(r'No\. of apolar-polar contacts: (\d+)', result.stdout)
        if apolar_polar_match:
            output_dict['apolar_polar_contacts'] = int(apolar_polar_match.group(1))
        
        # Extract apolar-apolar contacts
        apolar_apolar_match = re.search(r'No\. of apolar-apolar contacts: (\d+)', result.stdout)
        if apolar_apolar_match:
            output_dict['apolar_apolar_contacts'] = int(apolar_apolar_match.group(1))
        
        # Extract apolar NIS percentage
        apolar_nis_match = re.search(r'Percentage of apolar NIS residues: (\d+\.\d+)', result.stdout)
        if apolar_nis_match:
            output_dict['apolar_nis_percentage'] = float(apolar_nis_match.group(1))
        
        # Extract charged NIS percentage
        charged_nis_match = re.search(r'Percentage of charged NIS residues: (\d+\.\d+)', result.stdout)
        if charged_nis_match:
            output_dict['charged_nis_percentage'] = float(charged_nis_match.group(1))
        
        # Extract binding affinity
        affinity_match = re.search(r'Predicted binding affinity \(kcal\.mol-1\):\s+(-?\d+\.?\d*)', result.stdout)
        if affinity_match:
            output_dict['binding_affinity'] = float(affinity_match.group(1))
        
        # Extract dissociation constant
        kd_match = re.search(r'Predicted dissociation constant \(M\) at 25\.0˚C:\s+(\d+\.?\d*e-\d+)', result.stdout)
        if kd_match:
            output_dict['dissociation_constant'] = kd_match.group(1)
        
        return output_dict
            
    except Exception as e:
        error_msg = str(e)
        print(f"Error during prediction: {error_msg}")
        return {"error": f"Error during prediction: {error_msg}"}
    
    finally:
        # Clean up combined PDB file
        if os.path.exists(combined_pdb):
            os.remove(combined_pdb)

def extract_sequence_number(filename):
    """Extract sequence number from filename"""
    match = re.search(r'sequence_(\d+)_structure', filename)
    return int(match.group(1)) if match else None

def analyze_sequence_gaps(pdb_files):
    """Analyze gaps in sequence numbers"""
    sequence_numbers = set()
    for pdb_file in pdb_files:
        seq_num = extract_sequence_number(pdb_file)
        if seq_num is not None:
            sequence_numbers.add(seq_num)
    
    # Find min and max sequence numbers
    min_seq = min(sequence_numbers)
    max_seq = max(sequence_numbers)
    
    # Find missing sequence numbers
    all_sequences = set(range(min_seq, max_seq + 1))
    missing_sequences = sorted(all_sequences - sequence_numbers)
    
    return min_seq, max_seq, missing_sequences

def save_summary(results, output_dir="hepatitis_binding_results"):
    """Save summary of all predictions"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_file = os.path.join(output_dir, f"binding_affinity_summary_{timestamp}.txt")
    
    successful_predictions = [r for r in results if r['affinity'] is not None]
    failed_predictions = [r for r in results if r['affinity'] is None]
    
    if successful_predictions:
        affinities = [r['affinity'] for r in successful_predictions]
        min_affinity = min(affinities)
        max_affinity = max(affinities)
        avg_affinity = sum(affinities) / len(affinities)
    
    with open(summary_file, 'w') as f:
        f.write("Binding Affinity Prediction Summary\n")
        f.write("=================================\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        f.write(f"Total predictions attempted: {len(results)}\n")
        f.write(f"Successful predictions: {len(successful_predictions)}\n")
        f.write(f"Failed predictions: {len(failed_predictions)}\n\n")
        
        if successful_predictions:
            f.write("Binding Affinity Statistics:\n")
            f.write("-------------------------\n")
            f.write(f"Minimum affinity: {min_affinity:.2f} kcal/mol\n")
            f.write(f"Maximum affinity: {max_affinity:.2f} kcal/mol\n")
            f.write(f"Average affinity: {avg_affinity:.2f} kcal/mol\n\n")
        
        if failed_predictions:
            f.write("Failed Predictions:\n")
            f.write("-----------------\n")
            for result in failed_predictions:
                f.write(f"Sequence {result['sequence']}: {result['error']}\n")
        
        f.write("\nTop 10 Strongest Binding Sequences:\n")
        f.write("--------------------------------\n")
        top_10 = sorted(successful_predictions, key=lambda x: x['affinity'])[:10]
        for result in top_10:
            f.write(f"Sequence {result['sequence']}: {result['affinity']:.2f} kcal/mol\n")
    
    print(f"\nSummary saved to: {summary_file}")

def main():
    # Paths
    malaria_dir = "hepatitis_3d"
    antigen_path = "antigen_hep/model_01_HCV_NS5B.pdb"
    
    # Get all PDB files
    pdb_files = [f for f in os.listdir(malaria_dir) if f.endswith('.pdb') and f.startswith('sequence_')]
    total_files = len(pdb_files)
    
    # Analyze sequence gaps
    min_seq, max_seq, missing_sequences = analyze_sequence_gaps(pdb_files)
    
    print(f"\nSequence Analysis:")
    print(f"Total PDB files found: {total_files}")
    print(f"Sequence range: {min_seq} to {max_seq} (expected {max_seq - min_seq + 1} files)")
    print(f"Missing sequences: {len(missing_sequences)}")
    if missing_sequences:
        print("Missing sequence numbers:")
        # Print missing sequences in a compact format
        for i in range(0, len(missing_sequences), 10):
            chunk = missing_sequences[i:i+10]
            print(f"  {', '.join(map(str, chunk))}")
    
    # Process each PDB file
    print(f"\nStarting binding affinity predictions...")
    results = []
    
    for i, pdb_file in enumerate(sorted(pdb_files, key=extract_sequence_number), 1):
        seq_num = extract_sequence_number(pdb_file)
        print(f"\nProcessing file {i}/{total_files}: {pdb_file} (Sequence {seq_num})")
        
        if seq_num is None:
            print(f"Could not extract sequence number from {pdb_file}, skipping...")
            continue
        
        protein_path = os.path.join(malaria_dir, pdb_file)
        
        # Predict binding affinity
        affinity, prodigy_output, error = predict_binding_affinity(protein_path, antigen_path)
        
        # Store result
        results.append({
            'sequence': seq_num,
            'affinity': affinity,
            'error': error
        })
        
        # Save individual result
        save_results(seq_num, affinity, prodigy_output, protein_path, antigen_path, error)
        
        if affinity is not None:
            print(f"\nPredicted binding affinity: {affinity:.2f} kcal/mol")
        else:
            print(f"\nFailed to predict binding affinity: {error}")
    
    # Save summary
    save_summary(results)

if __name__ == "__main__":
    main()
