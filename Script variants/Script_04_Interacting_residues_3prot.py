import os
import pymol
from pymol import cmd

# Initialize PyMOL in a script-friendly mode
pymol.finish_launching(["pymol", "-qc"])

def analyze_protein_cif(file_path, output_file, pidt_threshold, distance_cutoff):
    """
    Analyzes interactions between residues in three chains of a protein structure.
    Appends results to a single text file.
    Parameters:
        file_path (str): Path to the .cif file.
        output_file (str): The text file to which results will be appended.
        pidt_threshold (float): Minimum pIDDT value for residue selection.
        distance_cutoff (float): Distance cutoff in Å for residue interaction.
    """
    try:
        # Extract protein names dynamically from the filename (format: fold_PROTEINA_PROTEINB_PROTEINC_model_0.cif)
        structure_name = os.path.splitext(os.path.basename(file_path))[0]
        parts = structure_name.split("_")
        
        # Assuming the protein names are in positions 1 and 2 (e.g., fold_PROTEINA_PROTEINB_model_0)
        proteinA = parts[1]  # PROTEINA is the second element in the split filename
        proteinB = parts[2]  # PROTEINB is the third element in the split filename
        proteinC = parts[3]  # PROTEINC is the fourth element in the split filename

        # Reset PyMOL state before loading each new structure
        cmd.reinitialize()
        
        # Load the structure
        cmd.load(file_path, structure_name)

        # Step 1: Select residues with b > 49.999 for both chains
        cmd.select("chainA_b_factor", f"chain A and b > {pidt_threshold}")
        cmd.select("chainB_b_factor", f"chain B and b > {pidt_threshold}")
        cmd.select("chainC_b_factor", f"chain C and b > {pidt_threshold}")

        # Step 2: Find residues in proximity (within distance_cutoff)
        cmd.select("close_A_to_B", f"byres (chainA_b_factor within {distance_cutoff} of chainB_b_factor)")
        cmd.select("close_A_to_C", f"byres (chainA_b_factor within {distance_cutoff} of chainC_b_factor)")
        cmd.select("close_B_to_A", f"byres (chainB_b_factor within {distance_cutoff} of chainA_b_factor)")
        cmd.select("close_B_to_C", f"byres (chainB_b_factor within {distance_cutoff} of chainC_b_factor)")
        cmd.select("close_C_to_A", f"byres (chainC_b_factor within {distance_cutoff} of chainA_b_factor)")
        cmd.select("close_C_to_B", f"byres (chainC_b_factor within {distance_cutoff} of chainB_b_factor)")
        
        # Step 3: Get the length of each protein (number of unique residues in each chain)
        chainA_residues = set(atom.resi for atom in cmd.get_model("chain A").atom)
        chainB_residues = set(atom.resi for atom in cmd.get_model("chain B").atom)
        chainC_residues = set(atom.resi for atom in cmd.get_model("chain C").atom)
        
        # Calculate the number of unique residues in each chain
        chainA_length = len(chainA_residues)
        chainB_length = len(chainB_residues)
        chainC_length = len(chainC_residues)
        
        # Extract residue IDs for residues meeting the criteria
        close_A_to_B_IDs = sorted(set(int(atom.resi) for atom in cmd.get_model("close_A_to_B").atom))
        close_A_to_C_IDs = sorted(set(int(atom.resi) for atom in cmd.get_model("close_A_to_C").atom))
        close_B_to_A_IDs = sorted(set(int(atom.resi) for atom in cmd.get_model("close_B_to_A").atom))
        close_B_to_C_IDs = sorted(set(int(atom.resi) for atom in cmd.get_model("close_B_to_C").atom))
        close_C_to_A_IDs = sorted(set(int(atom.resi) for atom in cmd.get_model("close_C_to_A").atom))
        close_C_to_B_IDs = sorted(set(int(atom.resi) for atom in cmd.get_model("close_C_to_B").atom))
        
        # Write results to the output file
        with open(output_file, "a") as f:
            f.write(f"Residues in {proteinA} (close to {proteinB}):\n")
            if close_A_to_B_IDs:
                f.write(", ".join(map(str, close_A_to_B_IDs)) + "\n")
            else:
                f.write("\n")
                
            f.write(f"Residues in {proteinA} (close to {proteinC}):\n")
            if close_A_to_C_IDs:
                f.write(", ".join(map(str, close_A_to_C_IDs)) + "\n")
            else:
                f.write("\n")
                
            f.write(f"Residues in {proteinB} (close to {proteinA}):\n")
            if close_B_to_A_IDs:
                f.write(", ".join(map(str, close_B_to_A_IDs)) + "\n")
            else:
                f.write("\n")
            
            f.write(f"Residues in {proteinB} (close to {proteinC}):\n")
            if close_B_to_C_IDs:
                f.write(", ".join(map(str, close_B_to_C_IDs)) + "\n")
            else:
                f.write("\n")
                
            f.write(f"Residues in {proteinC} (close to {proteinA}):\n")
            if close_C_to_A_IDs:
                f.write(", ".join(map(str, close_C_to_A_IDs)) + "\n")
            else:
                f.write("\n")
                
            f.write(f"Residues in {proteinC} (close to {proteinB}):\n")
            if close_C_to_B_IDs:
                f.write(", ".join(map(str, close_C_to_B_IDs)) + "\n")
            else:
                f.write("\n")
            
            
        print(f"Analysis complete for {file_path}. Results appended to {output_file}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
    finally:
        # Clean up PyMOL objects
        cmd.delete("all")    

def analyze_folder(input_folder, output_file, pidt_threshold, distance_cutoff):
    """
    Analyzes all .cif files in a folder and appends results to a single text file.
    Parameters:
        input_folder (str): Path to the folder containing .cif files.
        output_file (str): The text file to which results will be appended.
        pidt_threshold (float): Minimum pIDDT value for residue selection.
        distance_cutoff (float): Distance cutoff in Å for residue interaction.
    """
    cif_files = [f for f in os.listdir(input_folder) if f.endswith(".cif")]
    for cif_file in cif_files:
        file_path = os.path.join(input_folder, cif_file)
        analyze_protein_cif(file_path, output_file, pidt_threshold, distance_cutoff)

if __name__ == "__main__":
    # User-specified parameters
    input_folder = "cif_3prot"  # Replace with your .cif files folder
    output_folder = "Interacting residues_short_3p_70.txt"  # Replace with your desired output file
    pidt_threshold = 70  # Minimum pIDDT score
    distance_cutoff = 4.5  # Distance cutoff in Å

    # Run batch analysis
    analyze_folder(input_folder, output_folder, pidt_threshold, distance_cutoff)

    # Shutdown PyMOL when done
    cmd.quit()
