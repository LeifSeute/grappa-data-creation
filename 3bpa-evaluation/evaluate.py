#%%
from argparse import ArgumentParser
import numpy as np
from pathlib import Path
from grappa.utils.openff_utils import get_openmm_system
from grappa.utils.openmm_utils import get_energies
from grappa import OpenmmGrappa
from grappa.data.mol_data import MolData
from parse_ds import parse_data, to_arrays, get_dihedrals
from grappa.data.molecule import Molecule
from ase import Atoms
from grappa.utils.graph_utils import get_isomorphic_permutation
from openff.toolkit.topology import Molecule as OpenFFMolecule
import json
import pandas as pd
from collections import defaultdict
from scipy.stats import bootstrap


# Setup argument parser
parser = ArgumentParser()
parser.add_argument('--force-field', '-ff', type=str, default='gaff-2.11',
                    help='Force field to evaluate (.offxml/.xml) or base force field for grappa')
parser.add_argument('--grappa-ff', '-g', type=str, default=None,
                    help='Grappa force field tag or checkpoint to overwrite the base force field topology')
parser.add_argument('--output_file', '-d', type=str, default=None,
                    help='Output file to store calculated energies and forces')
parser.add_argument('--max_confs', '-n', type=int, default=None,
                    help='Number of conformers to generate')
args = parser.parse_args()


DATASETS = [
    'test_dih_beta180.xyz',
    'test_dih_beta150.xyz',
    'test_dih_beta120.xyz',
    'test_300K.xyz',
    'test_600K.xyz',
    'test_1200K.xyz',
    'train_300K.xyz'
]

# Function to load Grappa force field if provided
def load_forcefield(ff_tag):
    if ff_tag is None or ff_tag == 'None':
        return None
    if ff_tag.endswith('.ckpt'):
        return OpenmmGrappa.from_ckpt(ff_tag)
    else:
        return OpenmmGrappa.from_tag(ff_tag)


def get_system(mapped_smiles, atomic_numbers, positions, ff):

    # Generate ASE atoms objects
    atoms = Atoms(numbers=atomic_numbers, positions=positions[0])
    atomic_numbers_openff = [openff_mol.atoms[i].atomic_number for i in range(openff_mol.n_atoms)]
    bonds_openff = [(bond.atom1_index, bond.atom2_index) for bond in openff_mol.to_topology().bonds]
    atoms_openff = Atoms(numbers=atomic_numbers_openff)

    # Convert to Grappa Molecule and get permutation
    mol_ase = Molecule.from_ase(atoms)
    mol_openff = Molecule.from_ase(atoms_openff, bonds=bonds_openff)
    g_ase = mol_ase.to_dgl()
    g_openff = mol_openff.to_dgl()
    permutation = get_isomorphic_permutation(g_openff, g_ase)

    # Get the base force field system
    system, top, _ = get_openmm_system(
        mapped_smiles=mapped_smiles,
        openff_forcefield=ff
    )

    if grappa_ff:
        system = grappa_ff.parametrize_system(system, top)

    return system, permutation



# Function to calculate energies and forces for a given dataset
def calculate_energies_forces(dataset, system, permutation, max_confs=None):
    """
    Returns 
    """
    file_path = f'dataset_3BPA/{dataset}'
    data = parse_data(file_path, num_confs=max_confs)
    atomic_numbers, positions, forces, energies = to_arrays(data)

    positions_perm = positions[:, permutation]
    forces_perm = forces[:, permutation]

    # Calculate energies and forces
    energies_calculated, forces_calculated = get_energies(
        system,
        positions_perm
    )

    return energies, forces_perm, energies_calculated, forces_calculated


# Prepare output filename
if args.output_file:
    output_filename = args.output_file
else:
    base_ff = args.force_field.split('/')[-1].replace('.offxml', '').replace('.xml', '')
    if args.grappa_ff:
        grappa_tag = '_'.join(args.grappa_ff.split('/'))
        grappa_tag = grappa_tag.replace('.ckpt', '')
        output_filename = f"{grappa_tag}"
    else:
        output_filename = f"{base_ff}"

output_dir = Path('eval_data/'+output_filename)
output_dir.mkdir(exist_ok=True, parents=True)

print(f"Saving results to {str(output_dir)}")

grappa_ff = load_forcefield(args.grappa_ff)

# Load the reference molecule from SMILES
openff_mol = OpenFFMolecule.from_smiles('C1=CC=C(C=C1)COC2=C(N=CC=C2)N')
mapped_smiles = openff_mol.to_smiles(mapped=True)


ds_0 = DATASETS[0]
file_path = f'dataset_3BPA/{ds_0}'
data = parse_data(file_path, num_confs=args.max_confs)
atomic_numbers, positions, forces, energies = to_arrays(data)

system, permutation = get_system(mapped_smiles=mapped_smiles, atomic_numbers=atomic_numbers, positions=positions, ff=args.force_field)

# load results.json:
resultspath = Path('eval_data/results.json')
if resultspath.exists():
    with open(resultspath, 'r') as f:
        results = json.load(f)
else:
    results = defaultdict(dict)

# Loop over datasets and compute energies/forces
for dataset in DATASETS:
    # Get dihedral angles
    file_path = f'dataset_3BPA/{dataset}'
    data = parse_data(file_path, num_confs=args.max_confs)
    alpha, beta, gamma = get_dihedrals(data)

    # Calculate energies and forces
    energies, forces, energies_calculated, forces_calculated = calculate_energies_forces(dataset, system, permutation, max_confs=args.max_confs)

    output_file = output_dir / f"{dataset.replace('.xyz', '')}.npz"

    # Save energies, forces, and dihedral angles
    np.savez(output_file,
             alpha=alpha,
             beta=beta,
             gamma=gamma,
             energies=energies_calculated,
             forces=forces_calculated)

    # print(f"Results for {dataset} saved to {str(output_file)}")

    # Save QM energies and forces separately
    qm_output_dir = Path('eval_data/qm')
    qm_output_dir.mkdir(exist_ok=True, parents=True)
    qm_output_filename = f"{dataset.replace('.xyz', '')}"
    qm_output_file = qm_output_dir / f"{qm_output_filename}.npz"
    np.savez(qm_output_file,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            energies=energies,
            forces=forces)
    
    # calculate centered energy rmse:
    centered_energies = energies - np.mean(energies)
    centered_energies_calculated = energies_calculated - np.mean(energies_calculated)
    
    num_confs = forces.shape[0]

    # get mean and std by bootstrapping with scipy:
    def rmse(data):
        return np.sqrt(np.mean(data**2))
    
    def bootstrap(data, func, n=1000):
        samples = []
        for i in range(n):
            idx = np.random.randint(0, len(data), len(data))
            samples.append(func(data[idx]))
        return np.mean(samples), np.std(samples)

    e_diffs = centered_energies - centered_energies_calculated
    f_diffs = forces - forces_calculated

    e_rmse, e_rmse_err = bootstrap(e_diffs, rmse)
    f_crmse, f_crmse_err = bootstrap(f_diffs, rmse)

    d = {
        output_filename: {
            'energy_rmse': e_rmse,
            'force_crmse': f_crmse,
            'energy_rmse_err': e_rmse_err,
            'force_crmse_err': f_crmse_err,
            'num_confs': num_confs
        }
    }

    results[dataset].update(d)


# Save results to json
with open(resultspath, 'w') as f:
    json.dump(results, f)

# print a energy_accuracy.txt file with pandas:
# dictionary {ds:ff:energy_rmse}:
energy_rmse = {}
for ds, ff_data in results.items():
    energy_rmse[ds] = {'num_confs': ff_data[list(ff_data.keys())[0]]['num_confs']}
    for ff, data in ff_data.items():
        energy_rmse[ds][ff] = data['energy_rmse']

force_crmse = {}
for ds, ff_data in results.items():
    force_crmse[ds] = {}
    for ff, data in ff_data.items():
        force_crmse[ds][ff] = data['force_crmse']

df = pd.DataFrame(energy_rmse)

# get a nice table:
df = df.T
df = df.round(2)

s = df.to_string()
with open('eval_data/energy_accuracy.txt', 'w') as f:
    f.write(s)
print(f"Energy RMSE results:")
print(s)

df = pd.DataFrame(force_crmse)

df = df.T
df = df.round(2)

s = df.to_string()
with open('eval_data/force_accuracy.txt', 'w') as f:
    f.write(s)

print(f"Energy RMSE results saved to eval_data/energy_accuracy.txt")