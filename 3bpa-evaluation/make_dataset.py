'''
Transform the 3bpa dataset into grappa format for fine-tuning on it.
'''

from parse_ds import parse_data, to_arrays, get_dihedrals
from grappa.utils.data_utils import get_data_path

ds = 'dataset_3BPA'
name1 = f'train_300K.xyz'
name2 = 'test_dih_beta150.xyz'

def safe_dataset(name):

    file_path = f'{ds}/{name}'

    DSPATH = get_data_path()/'datasets'

    data = parse_data(file_path)
    atomic_numbers, positions, forces, energies = to_arrays(data)

    alpha, beta, gamma = get_dihedrals(data)

    import matplotlib.pyplot as plt

    plt.hist(energies, bins=50)

    from openff.toolkit.topology import Molecule

    smiles = 'C1=CC=C(C=C1)COC2=C(N=CC=C2)N'
    openff_mol = Molecule.from_smiles(smiles)
    mapped_smiles = openff_mol.to_smiles(mapped=True)


    from ase import Atoms
    import numpy as np

    atomic_numbers_openff = [openff_mol.atoms[i].atomic_number for i in range(openff_mol.n_atoms)]
    bonds_openff = openff_mol.to_topology().bonds
    bonds_openff = [(bond.atom1_index, bond.atom2_index) for bond in bonds_openff]

    atoms_openff = Atoms(numbers=atomic_numbers_openff)


    atoms = Atoms(numbers=atomic_numbers, positions=positions[0])


    from grappa.data.molecule import Molecule

    mol_ase = Molecule.from_ase(atoms)
    mol_openff = Molecule.from_ase(atoms_openff, bonds=bonds_openff)

    g_ase = mol_ase.to_dgl()
    g_openff = mol_openff.to_dgl()

    from grappa.utils.graph_utils import get_isomorphic_permutation

    permutation = get_isomorphic_permutation(g_openff, g_ase)

    atomic_numbers_openff = atomic_numbers[permutation]
    xyz_openff = positions[:, permutation]
    forces_openff = forces[:, permutation]


    from grappa.utils.openff_utils import get_openmm_system
    from grappa.data.mol_data import MolData
    from pathlib import Path
    from grappa import OpenmmGrappa

    import pickle
    from openmm import XmlSerializer

    if not Path('gaff-system.xml').exists():
        gaff_system, top, _ = get_openmm_system(mapped_smiles=mapped_smiles, openff_forcefield='gaff-2.11')
    else:
        with open('gaff-top.pkl', 'rb') as f:
            top = pickle.load(f)

        with open('gaff-system.xml', 'r') as f:
            gaff_system = XmlSerializer.deserialize(f.read())


    # save system and top:

    with open('gaff-system.xml', 'w') as f:
        f.write(XmlSerializer.serialize(gaff_system))

    with open('gaff-top.pkl', 'wb') as f:
        pickle.dump(top, f)

    # gaff_system, top, _ = get_openmm_system(mapped_smiles=mapped_smiles, openff_forcefield='openff_unconstrained-1.2.0.offxml')




    max_energy = 80.
    energies -= energies.min()
    allconfs = np.where(energies < max_energy)[0]

    # shuffle the confs:
    # np.random.seed(42)
    np.random.seed(1)
    np.random.shuffle(allconfs)

    train, val = np.split(allconfs, [int(.9*len(allconfs))])

    if 'train' in name:

        for confs, splitname in zip([train, val], ['train', 'val']):

            energies_ = energies[confs]
            xyz_openff_ = xyz_openff[confs]
            forces_openff_ = forces_openff[confs]


            moldata_gaff = MolData.from_openmm_system(gaff_system, top, xyz=xyz_openff_, gradient=-forces_openff_, energy=energies_, mol_id=smiles, ff_name='gaff-2.11')

            p = Path(DSPATH) / f'3bpa-{splitname}'

            moldata_gaff.save(p/f'3bpa-{name.strip(".xyz")}.npz')

    else:

        NUM_CONFS = 10

        confs = allconfs[:NUM_CONFS]

        energies_ = energies[confs]
        xyz_openff_ = xyz_openff[confs]
        forces_openff_ = forces_openff[confs]

        print('alpha: [', ', '.join([f'{a:.1f}' for a in alpha[confs]]), ']')
        
        print('gamma: [', ', '.join([f'{g:.1f}' for g in gamma[confs]]), ']')



        moldata_gaff = MolData.from_openmm_system(gaff_system, top, xyz=xyz_openff_, gradient=-forces_openff_, energy=energies_, mol_id=smiles, ff_name='gaff-2.11')

        p = Path(DSPATH) / f'3bpa-scan-{NUM_CONFS}'

        moldata_gaff.save(p/f'{name.strip(".xyz")}.npz') 

if __name__ == '__main__':
    safe_dataset(name1)
    safe_dataset(name2)