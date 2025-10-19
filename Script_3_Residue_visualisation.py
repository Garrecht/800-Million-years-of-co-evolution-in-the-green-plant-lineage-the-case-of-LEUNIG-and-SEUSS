import os
import pymol
from pymol import cmd

# Set the folder containing the protein files and the output directory
input_folder = "E:/HomeOffice_current/Computational/Interacting residues/cif"
output_folder = "E:/HomeOffice_current/Computational/Interacting residues/PyMol visualisation"

os.makedirs(output_folder, exist_ok=True)

def visualize_protein(file_path, output_name):
    # Load the protein complex
    cmd.reinitialize()
    cmd.load(file_path)
    protein_name = os.path.splitext(os.path.basename(file_path))[0]

    # Selection of domains
    cmd.select("AtLUG_LUFS", f"/{protein_name}//A and resi 1-160")
    cmd.select("AtSEU_inter", f"/{protein_name}//B and b>50 within 10 of AtLUG_LUFS")

    # Selection of whole residues where at least one atom is within 4.5 Å and b>50 of the other chain
    cmd.select("AtLUG_AS", "byres AtLUG_LUFS and b>50 within 4.5 of AtSEU_inter and b>50")
    cmd.select("AtSEU_AS", "byres AtSEU_inter and b>50 within 4.5 of AtLUG_LUFS and b>50")

    # Grouping selections
    cmd.group("AtLUG_SEU", f"{protein_name} AtLUG_LUFS AtLUG_AS AtSEU_AS AtSEU_inter")

    # Coloring
    cmd.bg_color("white")
    cmd.color("hafnium", f"/{protein_name}//A")
    cmd.color("actinium", f"/{protein_name}//A and backbone")
    cmd.color("orange", f"/{protein_name}//B")
    cmd.color("copper", f"/{protein_name}//B and backbone")
    cmd.color("grey", f"/{protein_name}/ and b<50")  


    # Visualization
    cmd.select("AtSEU_inter", "AtSEU_inter extend 20")
    cmd.hide("everything")
    cmd.show("cartoon", "AtLUG_LUFS AtSEU_inter")
    cmd.show("licorice", "AtLUG_AS AtSEU_AS")

    # Adjust view to fill canvas
    cmd.orient("AtLUG_AS")
    cmd.zoom("AtLUG_AS", buffer=10)
    cmd.set("fog", 0)  # Disables the fog effect
    
    # Export with desired resolution and orientation
    output_file = os.path.join(output_folder, f"{output_name}.png")
    cmd.ray(4800, 3000)
    cmd.png(output_file, dpi=600)

# Loop through all files in the folder
for filename in os.listdir(input_folder):
    if filename.endswith(".pdb") or filename.endswith(".cif"):
        file_path = os.path.join(input_folder, filename)
        output_name = os.path.splitext(filename)[0]
        visualize_protein(file_path, output_name)

print("Visualization complete. Images saved to:", output_folder)
