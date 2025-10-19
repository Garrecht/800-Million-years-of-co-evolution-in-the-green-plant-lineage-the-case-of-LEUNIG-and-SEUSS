import os
import csv
import pymol
from pymol import cmd

def collapse_residue_ranges(residues):
    """Collapse a list of domain positions into residue ranges."""
    if not residues:
        return []
    
    residues = sorted(map(int, residues))  # Convert residues to integers
    collapsed = []
    
    start = residues[0]
    end = residues[0]
    
    for i in range(1, len(residues)):
        if residues[i] == end + 1:
            end = residues[i]
        else:
            if start == end:
                collapsed.append((start, end))  # Single residue
            else:
                collapsed.append((start, end))  # Range
            start = residues[i]
            end = residues[i]
    
    # Append the last range or single residue
    collapsed.append((start, end))
    return collapsed

def extract_sheet_positions(cif_file, threshold=70):
    """Extract beta-sheet positions for chains with pIDDT >= threshold."""
    cmd.load(cif_file)
    
    sheet_data = []
    for chain in cmd.get_chains():
        cmd.select(f'chain_{chain}', f'chain {chain}')
        
        atoms = cmd.get_model(f'chain_{chain}')
        selected_residues = [atom.resi for atom in atoms.atom if atom.b >= threshold]
        
        if selected_residues:
            residue_selection = '+'.join(map(str, selected_residues))
            cmd.select(f'sheet_{chain}', f'chain {chain} and ss S and resi {residue_selection}')
            
            sheet_residues = cmd.get_model(f'sheet_{chain}')
            if sheet_residues:
                sheet_residue_numbers = [int(atom.resi) for atom in sheet_residues.atom]
                sheet_data.extend(collapse_residue_ranges(sheet_residue_numbers))
    
    cmd.delete("all")
    return sheet_data

def process_directory(cif_dir, output_file, threshold=70):
    """Process all CIF files in the given directory and write results to a table."""
    results = []
    
    for cif_file in os.listdir(cif_dir):
        if cif_file.endswith('.cif'):
            full_cif_path = os.path.join(cif_dir, cif_file)
            print(f"Processing {cif_file}...")
            sheet_data = extract_sheet_positions(full_cif_path, threshold)
            
            # Prepare row data
            row = [os.path.basename(cif_file)]
            for start, end in sheet_data:
                row.append(start)
                row.append(end)
            
            results.append(row)
    
    # Determine the maximum number of sheet columns required
    max_sheet_count = max((len(row) - 1) // 2 for row in results)
    
    # Create header
    header = ["Name"]
    for i in range(1, max_sheet_count + 1):
        header.append(f"Sheet {i} Start")
        header.append(f"Sheet {i} End")
    
    # Write results to a CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in results:
            # Pad rows with empty values if they have fewer helices
            while len(row) < len(header):
                row.append("")
            writer.writerow(row)

# Main function to run the processing
def main():
    cif_dir = "cif"  # Directory containing CIF files
    output_file = "results_sheet_70.csv"  # Name of output file
    threshold = 70  # pIDDT threshold for selecting residues (B-factor >= 70)

    process_directory(cif_dir, output_file, threshold)

if __name__ == "__main__":
    main()
