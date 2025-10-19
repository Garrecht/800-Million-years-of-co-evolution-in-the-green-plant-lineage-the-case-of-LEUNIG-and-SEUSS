import os
import pymol
from pymol import cmd

pymol.finish_launching(["pymol", "-qc"])

def analyze_protein_cif(file_path, output_file, pidt_threshold=70, distance_cutoff=5):
    
    try:
        # Extract protein names from the filename (format: fold_PROTEINA_PROTEINB_model_0.cif)
        structure_name = os.path.splitext(os.path.basename(file_path))[0]
        parts = structure_name.split("_")
        
        proteinA = parts[1]
        proteinB = parts[2]

        # Load the structure
        cmd.load(file_path, structure_name)

        # Step 1: Select residues with b > pidt threshold for both chains
        cmd.select("chainA_b_factor", f"chain A and b > {pidt_threshold}")
        cmd.select("chainB_b_factor", f"chain B and b > {pidt_threshold}")

        # Step 2: Find residues within distance_cutoff
        cmd.select("close_A", f"chainA_b_factor within {distance_cutoff} of chainB_b_factor")
        cmd.select("close_B", f"chainB_b_factor within {distance_cutoff} of chainA_b_factor")

        # Step 3: Get the length of each protein
        chainA_residues = set(atom.resi for atom in cmd.get_model("chain A").atom)
        chainB_residues = set(atom.resi for atom in cmd.get_model("chain B").atom)

        # Calculate the number of unique residues in each chain
        chainA_length = len(chainA_residues)
        chainB_length = len(chainB_residues)

        # Extract residue IDs for residues meeting the criteria above
        close_A_residues = sorted(set(int(atom.resi) for atom in cmd.get_model("close_A").atom))
        close_B_residues = sorted(set(int(atom.resi) for atom in cmd.get_model("close_B").atom))

        # Write results to the output file
        with open(output_file, "a") as f:
            f.write(f"Residues in {proteinA} (close to {proteinB}):\n")
            if close_A_residues:
                f.write(", ".join(map(str, close_A_residues)) + "\n")
            else:
                f.write("None\n")
            
            f.write(f"Residues in {proteinB} (close to {proteinA}):\n")
            if close_B_residues:
                f.write(", ".join(map(str, close_B_residues)) + "\n")
            else:
                f.write("None\n")


        print(f"Analysis complete for {file_path}. Results appended to {output_file}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
    finally:
        # Clean up PyMOL objects
        cmd.delete("all")

def analyze_folder(input_folder, output_file, pidt_threshold=50, distance_cutoff=5):
   
    cif_files = [f for f in os.listdir(input_folder) if f.endswith(".cif")]
    for cif_file in cif_files:
        file_path = os.path.join(input_folder, cif_file)
        analyze_protein_cif(file_path, output_file, pidt_threshold, distance_cutoff)

if __name__ == "__main__":
    # User-specified parameters
    input_folder = "cif"  # Path to the .cif file.
    output_file = "Interacting residues_short.txt"  # Output file
    pidt_threshold = 50  # Minimum pIDDT score
    distance_cutoff = 4.5  # Distance cutoff in Å for residue interaction

    # Run batch analysis
    analyze_folder(input_folder, output_file, pidt_threshold, distance_cutoff)

    # Shutdown PyMOL when done
    cmd.quit()
