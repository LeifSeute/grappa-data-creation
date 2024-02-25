#%%
from pathlib import Path
import shutil
from typing import Callable

FILENAMES = ['positions.npy', 'atomic_numbers.npy', 'charge.npy', 'smiles.npy', 'mapped_smiles.npy', 'pep.pdb', 'openmm_forces.npy', 'openmm_energies.npy']

#%%

def add_to_ds(source_folder:Path, target_folder:Path, skip_existing:bool=True, molname_rule:Callable=lambda name: name, filenames:list=FILENAMES):
    """
    Adds all pdb files from source_folder to target_folder. If skip_existing is True, it will skip files that are already in the target_folder.
    """
    for subdir in source_folder.iterdir():
        if not subdir.is_dir():
            continue
        
        new_molname = molname_rule(subdir.stem)
        new_subdir = target_folder/new_molname
        new_subdir.mkdir(exist_ok=True, parents=True)

        for filename in filenames:
            if not (subdir/filename).exists():
                # raise FileNotFoundError(f"File {filename} does not exist in {subdir}.")
                print(f"File {filename} does not exist in {subdir}.")
                continue
            if skip_existing and (new_subdir/filename).exists():
                continue
            shutil.copy(subdir/filename, new_subdir/filename)

def check_dataset_full(target_folder, filenames=FILENAMES):
    """
    Checks whether all files are present fro every molecule in the dataset.
    """
    mols_with_missing_files = 0
    for subdir in target_folder.iterdir():
        if not subdir.is_dir():
            continue
        for filename in filenames:
            if not (subdir/filename).exists():
                print(f"File {filename} does not exist in {subdir}.")
                mols_with_missing_files += 1
                break
    print(f"\n{mols_with_missing_files} molecules have missing files.\n")

#%%
import argparse
parser = argparse.ArgumentParser(description='Merge uncapped datasets.')
parser.add_argument('--full', action='store_true', help='Merge the full dataset.')
parser.add_argument('--prepare', action='store_true', help='Pre merge, where we merge the dipeptides and the mono peptides separately.')
parser.add_argument('--full_from_prepared', action='store_true', help='Merge the full dataset.')

args = parser.parse_args()

assert not sum([args.full, args.prepare, args.full_from_prepared]) > 1, "Only one of --full, --prepare, --full_from_prepared can be True."

data_folder = Path(__file__).parent/'data'

ace_rule = lambda name: 'B'+name
nme_rule = lambda name: name+'Z'

if args.prepare:

    target_folder1 = data_folder/'uncapped_monopeptides'
    target_folder2 = data_folder/'uncapped_dipeptides'

    source_dirs1 = ['singly_capped_ace', 'singly_capped_nme']
    source_dirs2 = ['dipep_ace_singly_capped', 'dipep_nme_singly_capped']

    rules = [ace_rule, nme_rule]

    for source_folder, rule in zip([data_folder/source_dir for source_dir in source_dirs1], rules):
        add_to_ds(source_folder, target_folder1, molname_rule=rule)

    for source_folder, rule in zip([data_folder/source_dir for source_dir in source_dirs2], rules):
        add_to_ds(source_folder, target_folder2, molname_rule=rule)

    check_dataset_full(target_folder1)
    check_dataset_full(target_folder2)

elif args.full_from_prepared:
    # now merge the two prepared datasets, using FILENAMES + ['psi4_energies.npy', 'psi4_forces.npy']

    target_folder = data_folder/'uncapped_dataset'
    target_folder.mkdir(exist_ok=True, parents=True)

    source_folders = [data_folder/'uncapped_monopeptides', data_folder/'uncapped_dipeptides']

    for source_folder in source_folders:
        add_to_ds(source_folder, target_folder, skip_existing=True, filenames=FILENAMES+['psi4_energies.npy', 'psi4_forces.npy'])

    check_dataset_full(target_folder, filenames=FILENAMES+['psi4_energies.npy', 'psi4_forces.npy'])

elif args.full:
    # merge the full dataset from the source folders, using FILENAMES
    target_folder = data_folder/'uncapped_dataset'
    target_folder.mkdir(exist_ok=True, parents=True)

    source_folders = [data_folder/'singly_capped_ace', data_folder/'singly_capped_nme', data_folder/'dipep_ace_singly_capped', data_folder/'dipep_nme_singly_capped']

    rules = [ace_rule, nme_rule, ace_rule, nme_rule]

    for source_folder, rule in zip(source_folders, rules):
        add_to_ds(source_folder, target_folder, skip_existing=True, filenames=FILENAMES, molname_rule=rule)

    check_dataset_full(target_folder)
