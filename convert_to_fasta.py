
#!/usr/bin/env python
from rdkit import Chem
import sys

def convert_to_fasta(input_pdb, output_fasta):
    mol = Chem.MolFromPDBFile(input_pdb)
    if mol is None:
        print(f"Failed to read the PDB file: {input_pdb}")
        sys.exit(1)
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    sequence = Chem.MolToSequence(mol)
    fasta = ">Molecule\n" + sequence + "\n"
    with open(output_fasta, 'w') as file:
        file.write(fasta)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python convert_to_fasta.py input.pdb output.fasta")
        sys.exit(1)
    input_pdb = sys.argv[1]
    output_fasta = sys.argv[2]
    convert_to_fasta(input_pdb, output_fasta)

