from ase import Atoms
from ase.calculators.psi4 import Psi4
import numpy as np

import os
from pathlib import Path


from pathlib import Path

from utils import Logger
###################
import sys
import os

###################

# define AlreadyBusy error type:
class AlreadyBusy(Exception):
    pass


def calc_states(pdb_folder, n_states=None, memory=32, num_threads=4, skip_if_busy=False):
    """
    Calculates the energies and forces for all states defined by positions.npy, atomic_numbers.npy and charge.npy in the given folder.
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
    
    if skip_if_busy:
        if (pdb_folder/Path("psi4_energies.npy")).exists() and (pdb_folder/Path("psi4_forces.npy")).exists():
            log(f"skipping {pdb_folder.stem}, already (partially) calculated")
            raise AlreadyBusy(f"skipping {pdb_folder.stem}, already (partially) calculated")

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


  
    # Calculate energies and forces using Psi4
    psi4_energies = []
    psi4_forces = []

    # load if present:
    if (pdb_folder/Path("psi4_energies.npy")).exists() and (pdb_folder/Path("psi4_forces.npy")).exists():
        psi4_energies = [e for e in np.load(str(pdb_folder/Path("psi4_energies.npy")))]
        psi4_forces = [f for f in np.load(str(pdb_folder/Path("psi4_forces.npy")))]
        

    from time import time
    start = time()

    missing_indices = range(len(psi4_energies), len(positions))

    # restroct the amount of calculations:
    if n_states is not None:
        if len(missing_indices) > n_states:
            missing_indices = missing_indices[:n_states]

    log(f"calculating {len(missing_indices)} states using the config\n")
    log(f"\tMETHOD: {METHOD}")
    log(f"\tBASIS: {BASIS}")
    log(f"\tMEMORY: {MEMORY}")
    log(f"\tNUM_THREADS: {NUM_THREADS}")
    log(f"\ttotal_charge: {total_charge}")
    log(f"\tmultiplicity: {multiplicity}")


    for num_calculated, i in enumerate(missing_indices):
        
        msg = f"calculating state number {i}/{len(positions)-1}... Progress: {num_calculated}/{len(missing_indices)-1}, time elapsed: {round((time() - start)/60., 2)} min"
        if num_calculated > 0:
            msg += f", avg time per state: {round((time() - start)/(num_calculated) / 60.,2)} min"
        log(msg)

        # Read the configuration
        atoms = Atoms(numbers=atomic_numbers, positions=positions[i])

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


        psi4_energies.append(energy)
        psi4_forces.append(forces)

        # save the energies and forces in every step:
        np.save(str(pdb_folder/Path("psi4_energies.npy")), np.array(psi4_energies))
        np.save(str(pdb_folder/Path("psi4_forces.npy")), np.array(psi4_forces))




def calc_all_states(folder, n_states=None, skip_errs=False, memory=32, num_threads=8, permute_seed=None, is_cleanup_run:list=[False, True]):
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
        calc_states(folder, n_states=n_states, memory=memory, num_threads=num_threads, skip_if_busy=False)
        return

    # iterate once with skip_if_busy=False and once with skip_if_busy=True (to finished the calculations that were potentially interrupted)
    for skip_if_busy in is_cleanup_run:

        skip_if_busy = not bool(skip_if_busy)

        # also recalculate the folders in case something was added

        pdb_folders = [f for f in Path(folder).iterdir() if f.is_dir()]
        random.shuffle(pdb_folders)
        log(f"calculating states for {len(pdb_folders)} folders in a first iteration.")

        for i, pdb_folder in enumerate(pdb_folders):
            log("")
            log(f"calculating states for {i}, {Path(pdb_folder).stem}...")
            try:
                calc_states(pdb_folder, n_states=n_states, memory=memory, num_threads=num_threads, skip_if_busy=skip_if_busy)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                if type(e) == AlreadyBusy:
                    continue
                if not skip_errs:
                    raise
                log(f"failed to calculate states for {i} ({Path(folder).stem}): {type(e)}\n: {e}")

        log(f"finished the iteration with skip_if_busy={skip_if_busy}.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculates states for a given folder.')
    parser.add_argument('folder', type=str, help='The folder containing the PDB files.')
    parser.add_argument('--n_states', '-n', type=int, help='The number of states to calculate.', default=None)
    parser.add_argument('--skip_errs', '-s', action='store_true', help='Skip errors.', default=False)
    parser.add_argument('--permute_seed', '-p', type=int, help='The seed to use for shuffling the folders.', default=None)
    parser.add_argument('--memory', '-m', type=int, help='The amount of memory to use.', default=32)
    parser.add_argument('--num_threads', '-t', type=int, help='The number of threads to use.', default=4)
    parser.add_argument('--cleanup_runs', '-c', type=bool, nargs='+', help='If True, the function will skip the folders where a calculation has already been started.', default=[False, True, True])
    args = parser.parse_args()
    calc_all_states(folder=args.folder, n_states=args.n_states, skip_errs=args.skip_errs, memory=args.memory, num_threads=args.num_threads, permute_seed=args.permute_seed, is_cleanup_run=args.cleanup_runs)