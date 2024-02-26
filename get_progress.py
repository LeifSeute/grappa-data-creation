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
    total_molecules = 0
    finished_molecules = 0
    unfinished_num_states_left = {}
    
    # Sum up shapes[0] for positions.npy files
    for positions_file in directory.glob(positions_pattern):
        data = np.load(positions_file)
        total_positions += data.shape[0]
        total_molecules += 1
    
    # Sum up shapes[0] for psi4_energies.npy files
    for energies_file in directory.glob(energies_pattern):
        data = np.load(energies_file)
        num_finite_energies = np.sum(np.isfinite(data))
        total_energies += num_finite_energies
        # HARD CODED THAT THE NUMBER OF STATES IS 50, change later
        if num_finite_energies == 50:
            finished_molecules += 1
        else:
            unfinished_num_states_left[energies_file.parent.stem] = 50-num_finite_energies
    
    # Calculate and report the results
    ratio = (total_energies / total_positions) * 100 if total_positions > 0 else 0
    print(f"Total shape[0] of positions.npy files: {total_positions}")
    print(f"Total shape[0] of psi4_energies.npy files: {total_energies}")
    print(f"Ratio of energies to positions (in percent): {ratio:.2f}%")
    print(f"Total molecules: {total_molecules}")
    print(f"Finished molecules: {finished_molecules}")
    print(f"States left: {unfinished_num_states_left}")


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
