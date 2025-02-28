'''
function that loads an openff molecule from a pdb file and write the smiles string into smiles.npy in the same folder
'''

from openff.toolkit.topology import Molecule
from pathlib import Path
import numpy as np

def write_smiles(pdb_folder):
    pdb_folder = Path(pdb_folder)
    pdb_file = pdb_folder/"pep.pdb"
    if not pdb_file.exists():
        raise FileNotFoundError(f"{pdb_file} does not exist.")
    mol = Molecule.from_polymer_pdb(pdb_file)
    smiles = mol.to_smiles(mapped=False)
    mapped_smiles = mol.to_smiles(mapped=True)
    np.save(str(pdb_folder/"smiles.npy"), smiles)
    np.save(str(pdb_folder/"mapped_smiles.npy"), mapped_smiles)


def write_smiles_in_folder(folder):
    folder = Path(folder)
    for pdb_folder in folder.iterdir():
        if pdb_folder.is_dir():
            write_smiles(pdb_folder)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Write the smiles string of the molecule in the pdb file to smiles.npy in the same folder.')
    parser.add_argument('folder', type=str, help='The folder containing the pdb files.')
    args = parser.parse_args()

    write_smiles_in_folder(args.folder)