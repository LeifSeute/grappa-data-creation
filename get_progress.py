import argparse
from pathlib import Path
import numpy as np

def sum_shapes_and_calculate_ratio(directory):
    # Paths to search for
    positions_pattern = '**/positions.npy'
    energies_pattern = '**/psi4_energies.npy'
    
    # Initialize counters
    total_positions = 0
    total_energies = 0
    
    # Sum up shapes[0] for positions.npy files
    for positions_file in directory.glob(positions_pattern):
        data = np.load(positions_file)
        total_positions += data.shape[0]
    
    # Sum up shapes[0] for psi4_energies.npy files
    for energies_file in directory.glob(energies_pattern):
        data = np.load(energies_file)
        total_energies += data.shape[0]
    
    # Calculate and report the results
    ratio = (total_energies / total_positions) * 100 if total_positions > 0 else 0
    print(f"Total shape[0] of positions.npy files: {total_positions}")
    print(f"Total shape[0] of psi4_energies.npy files: {total_energies}")
    print(f"Ratio of energies to positions (in percent): {ratio:.2f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sum shapes of .npy files and calculate ratio.')
    parser.add_argument('dir', type=str, help='Directory to search within.')
    args = parser.parse_args()

    # Convert the input argument to a Path object
    directory = Path(args.dir)
    if directory.exists() and directory.is_dir():
        sum_shapes_and_calculate_ratio(directory)
    else:
        print("The specified directory does not exist or is not a directory.")
