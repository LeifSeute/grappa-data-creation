import numpy as np
from pathlib import Path

def add_nans(path):
    for subdir in path.iterdir():
        if subdir.is_dir():
            if (subdir/"psi4_energies.npy").exists():
                energies = np.load(subdir/"psi4_energies.npy")
                forces = np.load(subdir/"psi4_forces.npy")
                if energies.shape[0] < 50:
                    energies = np.concatenate((energies, np.full((50-energies.shape[0]), np.nan)), axis=0)
                    forces = np.concatenate((forces, np.full((50-forces.shape[0], forces.shape[1], forces.shape[2]), np.nan)), axis=0)
                    np.save(subdir/"psi4_energies.npy", energies)
                    np.save(subdir/"psi4_forces.npy", forces)
            else:
                # create the files
                energies = np.full((50), np.nan)
                forces = np.full((50, 1, 3), np.nan)
                np.save(subdir/"psi4_energies.npy", energies) 
                np.save(subdir/"psi4_forces.npy", forces)

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add nans to psi4_energies.npy and psi4_forces.npy')
    parser.add_argument('dir', type=str, help='Directory to search within.')
    args = parser.parse_args()

    # Convert the input argument to a Path object
    directory = Path(args.dir)
    if directory.exists() and directory.is_dir():
        add_nans(directory)
    else:
        print("The specified directory does not exist or is not a directory.")