import os
import pymol
from pymol import cmd

# Initialize PyMOL in a script-friendly mode
pymol.finish_launching(["pymol", "-qc"])

def analyze_protein_cif(file_path, output_folder, pidt_threshold, distance_cutoff):
    """
    Vizualize interacting residues between three proteins in a complex
    Save a .pse PyMOL session file for each structure.
    
    Parameters:
        file_path (str): Path to the .cif file.
        output_folder (str): Folder where .pse session files will be saved.
        pidt_threshold (float): Minimum pIDDT value for residue selection.
        distance_cutoff (float): Distance cutoff in Å for residue interaction.
    """
    try:
        # Extract protein names dynamically from the filename (format: fold_PROTEINA_PROTEINB_PROTEINC_model_0.cif)
        structure_name = os.path.splitext(os.path.basename(file_path))[0]
        parts = structure_name.split("_")
        
        # Assuming the protein names are in positions 2 and 3 and 4 (e.g., fold_PROTEINA_PROTEINB_PROTEINC_model_0)
        proteinA = parts[1]  # PROTEINA is the second element in the split filename
        proteinB = parts[2]  # PROTEINB is the third element in the split filename
        proteinC = parts[3]  # PROTEINC is the fourth element in the split filename

        # Reset PyMOL state before loading each new structure
        cmd.reinitialize()
        
        # Load the structure
        cmd.load(file_path, structure_name)

        # Step 1: Select residues with b > 49.999 for all chains
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
        
        # Step 3: Color chains and residues
        cmd.bg_color("white")
        cmd.color("hafnium", "chain A")
        cmd.color("actinium", "chain A and backbone")
        cmd.color("orange", "chain B")
        cmd.color("copper", "chain B and backbone")
        cmd.color("barium", "chain C")
        cmd.color("forest", "chain C and backbone")
        cmd.color("grey", f"b < {pidt_threshold}")
        
        # Step 4: Visualization
        cmd.select("chainA_b_factor", f"chainA_b_factor extend 30")
        cmd.select("chainB_b_factor", f"chainB_b_factor extend 30")
        cmd.select("chainC_b_factor", f"chainC_b_factor extend 30")
        
        cmd.hide("everything")
        cmd.show("cartoon", "chainA_b_factor chainB_b_factor chainC_b_factor")
        cmd.show("licorice", "close_A_to_B close_A_to_C close_B_to_A close_B_to_C close_C_to_A close_C_to_B")

        # Step 5: Save session as .pse file
        os.makedirs(output_folder, exist_ok=True)
        pse_filename = f"{structure_name}.pse"
        pse_path = os.path.join(output_folder, pse_filename)
        cmd.save(pse_path)
        print(f"Saved session: {pse_path}")
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

def analyze_folder(input_folder, output_folder, pidt_threshold, distance_cutoff):
    # Performs the vizualization for all .cif files in a folder
    cif_files = [f for f in os.listdir(input_folder) if f.endswith(".cif")]
    for cif_file in cif_files:
        file_path = os.path.join(input_folder, cif_file)
        analyze_protein_cif(file_path, output_folder, pidt_threshold, distance_cutoff)

if __name__ == "__main__":
    # User-specified parameters
    input_folder = "cif_3prot"  # Replace with your .cif files folder
    output_folder = "pymol_sessions"  # Replace with your desired output file
    pidt_threshold = 50  # Minimum pIDDT score
    distance_cutoff = 4.5  # Distance cutoff in Å

    # Run batch analysis
    analyze_folder(input_folder, output_folder, pidt_threshold, distance_cutoff)

    # Shutdown PyMOL when done
    cmd.quit()
