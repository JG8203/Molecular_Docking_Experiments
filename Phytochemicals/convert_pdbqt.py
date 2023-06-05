
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def convert_sdf_to_pdbqt():
    # Get the current working directory
    input_directory = os.getcwd()

    # Create the output directory if it doesn't exist
    output_directory = os.path.join(input_directory, "pdbqt")
    os.makedirs(output_directory, exist_ok=True)

    # Get a list of all SDF files in the input directory
    sdf_files = [f for f in os.listdir(input_directory) if f.endswith(".sdf")]
    
    # Iterate over each SDF file
    for sdf_file in sdf_files:
        # Construct the input and output file paths
        input_file = os.path.join(input_directory, sdf_file)
        pdbqt_file = os.path.splitext(sdf_file)[0] + ".pdbqt"
        output_file = os.path.join(output_directory, pdbqt_file)
        
        # Load the SDF file
        supplier = Chem.SDMolSupplier(input_file)
        
        # Iterate over each molecule in the SDF file
        for molecule in supplier:
            # Generate 3D coordinates and assign charges
            AllChem.ComputeGasteigerCharges(molecule)
            AllChem.EmbedMolecule(molecule)
            AllChem.UFFOptimizeMolecule(molecule)
            
            # Write the molecule to a PDBQT file
            writer = Chem.rdmolfiles.PDBQTWriter(output_file)
            writer.write(molecule)
            writer.close()

# Usage
convert_sdf_to_pdbqt()
