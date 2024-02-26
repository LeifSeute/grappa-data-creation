from ase import Atoms
from ase.calculators.psi4 import Psi4
import numpy as np

from time import time



from pathlib import Path

from utils import Logger
###################
import sys
import os

###################


def calc_state(pdb_folder, memory=32, num_threads=4):
    """
    Calculates the energy and forces for one state defined by positions.npy, atomic_numbers.npy and charge.npy in the given folder. Several of these calls can run in parallel since the function signalizes that it is working on a state by writing inf to the psi4_energies.npy and psi4_forces.npy files.
    """

    log = Logger(Path(pdb_folder).parent, print_to_screen=True)

    METHOD = 'bmk'
    BASIS = '6-311+G(2df,p)'

    scratch = Path(__file__).parent/"psi_scratch"
    if not scratch.exists():
        os.makedirs(str(scratch), exist_ok=True)
    os.environ['PSI_SCRATCH'] = str(scratch)


    ACCURACY = 1e-2 # numerical threshold, in kcal/mol, we have 1kcal/mol ~50 meV

    if not memory is None and memory > 0:
        MEMORY = f'{int(memory)}GB'
    else:
        MEMORY = None
        
    if not num_threads is None and num_threads > 0:
        NUM_THREADS=num_threads
    else:
        NUM_THREADS = None


    pdb_folder = Path(pdb_folder)

    if not (pdb_folder/Path("positions.npy")).exists():
        return
    

    positions = np.load(str(pdb_folder/Path("positions.npy")))
    atomic_numbers = np.load(str(pdb_folder/Path("atomic_numbers.npy")))
    
    total_charge = np.load(str(pdb_folder/Path("charge.npy")))
    
    if not total_charge.shape == (1,):
        raise ValueError(f"total_charge.shape must be (1,), is: {total_charge.shape}")
    total_charge = int(total_charge[0])

    if not np.isclose(total_charge, round(total_charge,0), atol=1e-5):
        raise ValueError(f"total_charge is no integer: {total_charge}")
    

    multiplicity = 1

    if (pdb_folder/Path("multiplicity.npy")).exists():
        multiplicity = np.load(str(pdb_folder/Path("multiplicity.npy")))
        if not multiplicity.shape == (1,):
            raise ValueError(f"multiplicity.shape must be (1,), is: {multiplicity.shape}")
        multiplicity = int(multiplicity[0])



    # load if present:
    if (pdb_folder/Path("psi4_energies.npy")).exists():
        psi4_energies = np.load(str(pdb_folder/Path("psi4_energies.npy")))
    else:
        psi4_energies = np.zeros_like(positions[:,0,0])*np.nan

    if (pdb_folder/Path("psi4_forces.npy")).exists():       
        psi4_forces = np.load(str(pdb_folder/Path("psi4_forces.npy")))
    else:
        psi4_forces = np.zeros_like(positions)*np.nan

    if np.all(np.isfinite(psi4_energies)) and np.all(np.isfinite(psi4_forces)):
        log(f"all states have been calculated for {pdb_folder.stem}")
        return


    # now pick a state index for an uncalculated state, i.e. an index in the energies array where the energy is either nan or inf:
    state_index = np.where(np.isnan(psi4_energies) | np.isinf(psi4_energies))[0][0]

    # assert that all force entries are nan of inf too
    assert np.all(np.isnan(psi4_forces[state_index])) or np.all(np.isinf(psi4_forces[state_index]))

    # now store the arrays with inf at the state index (to signal that the state is being calculated)
    psi4_energies[state_index] = np.inf
    psi4_forces[state_index] = np.ones_like(psi4_forces[state_index])*np.inf

    np.save(str(pdb_folder/Path("psi4_energies.npy")), psi4_energies)
    np.save(str(pdb_folder/Path("psi4_forces.npy")), psi4_forces)


    start = time()


    log(f"calculating state using the config\n")
    log(f"\tMETHOD: {METHOD}")
    log(f"\tBASIS: {BASIS}")
    log(f"\tMEMORY: {MEMORY}")
    log(f"\tNUM_THREADS: {NUM_THREADS}")
    log(f"\ttotal_charge: {total_charge}")
    log(f"\tmultiplicity: {multiplicity}")


    
    msg = f"calculating state number {state_index}..."
    
    start = time()

    # Read the configuration
    atoms = Atoms(numbers=atomic_numbers, positions=positions[state_index])

    ###################
    # set up the calculator:
    kwargs = {"atoms":atoms, "method":METHOD, "basis":BASIS, "charge":total_charge, "multiplicity":1, "d_convergence":ACCURACY*23.06}

    if not MEMORY is None:
        kwargs["memory"] = MEMORY
    if not NUM_THREADS is None:
        kwargs["num_threads"] = NUM_THREADS

    atoms.set_calculator(Psi4(atoms=atoms, method=METHOD, memory=MEMORY, basis=BASIS, num_threads=NUM_THREADS, charge=total_charge, multiplicity=multiplicity))
    ###################

    energy = atoms.get_potential_energy(apply_constraint=False) # units: eV
    forces = atoms.get_forces(apply_constraint=False) # units: eV/Angstrom

    EV_IN_KCAL = 23.0609

    energy = energy * EV_IN_KCAL
    forces = forces * EV_IN_KCAL

    print(f"time elapsed: {round((time() - start)/60., 2)} min")


    # load the energies and forces again (another process might have written to them in the meantime)
    psi4_energies = np.load(str(pdb_folder/Path("psi4_energies.npy")))
    psi4_forces = np.load(str(pdb_folder/Path("psi4_forces.npy")))

    psi4_energies[state_index] = energy
    psi4_forces[state_index] = forces

    np.save(str(pdb_folder/Path("psi4_energies.npy")), psi4_energies)
    np.save(str(pdb_folder/Path("psi4_forces.npy")), psi4_forces)


def has_uncalculated_states(pdb_folder):
    """
    Returns True if there are any states with nan or inf in the psi4_energies.npy file or if the file does not exist.
    """
    pdb_folder = Path(pdb_folder)
    if not (pdb_folder/Path("psi4_energies.npy")).exists():
        return True
    psi4_energies = np.load(str(pdb_folder/Path("psi4_energies.npy")))
    return np.any(np.isnan(psi4_energies) | np.isinf(psi4_energies))


def calc_all_states(folder, skip_errs=False, memory=32, num_threads=8, permute_seed=None):
    """
    For all folders in the given folder, call calc_states.
    is_cleanup_run: list of bools. for every entry, iterates over all sub folders once. If True, the function will skip the folders where a calculation has already been started.
    """
    from pathlib import Path
    import random

    if permute_seed is not None:
        random.seed(permute_seed)

    log = Logger(Path(folder), print_to_screen=True)

    # if the folder itself contains the files, calculate them:
    if (Path(folder)/"positions.npy").exists():
        print(f"calculating states for {Path(folder).stem}...")
        while has_uncalculated_states(folder):
            calc_state(folder, memory=memory, num_threads=num_threads)
        return

    # iterate two times over the folders, since there might have been some errors
    
    for i in range(2):

        pdb_folders = [f for f in Path(folder).iterdir() if f.is_dir()]
        random.shuffle(pdb_folders)
        log(f"calculating states for {len(pdb_folders)} folders in a first iteration.")

        # iterate over the pdb_folders until there are no uncalculated states left:

        for i, pdb_folder in enumerate(pdb_folders):
            log("")
            log(f"calculating states for {i}, {Path(pdb_folder).stem}...")
            while has_uncalculated_states(pdb_folder):
                try:
                    calc_state(pdb_folder, memory=memory, num_threads=num_threads)
                except KeyboardInterrupt:
                    raise
                except Exception as e:
                    if not skip_errs:
                        raise
                    log(f"failed to calculate states for {i} ({Path(folder).stem}): {type(e)}\n: {e}")
                    break


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculates states for a given folder.')
    parser.add_argument('folder', type=str, help='The folder containing the PDB files.')
    parser.add_argument('--skip_errs', '-s', action='store_true', help='Skip errors.', default=False)
    parser.add_argument('--permute_seed', '-p', type=int, help='The seed to use for shuffling the folders.', default=None)
    parser.add_argument('--memory', '-m', type=int, help='The amount of memory to use.', default=32)
    parser.add_argument('--num_threads', '-t', type=int, help='The number of threads to use.', default=4)
    args = parser.parse_args()
    calc_all_states(folder=args.folder, skip_errs=args.skip_errs, memory=args.memory, num_threads=args.num_threads, permute_seed=args.permute_seed)