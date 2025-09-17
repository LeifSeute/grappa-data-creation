#%%
from pathlib import Path
from grappa.utils.data_utils import get_moldata_path
from grappa.data import MolData
from tqdm.auto import tqdm
import numpy as np
import dgl
import torch
import os

# --- a) Construct graphs from smiles strings ---
from grappa.utils.openff_utils import get_openff_molecule
from grappa.utils.graph_utils import as_nx
from grappa.data import Molecule


N_MAX = None


pdbfiles = list(Path('spice-pdbs').glob('*.pdb'))
spicepaths = list(Path('spice-npz').glob('*.npz'))

if N_MAX is not None:
    spicepaths = spicepaths[:N_MAX]
    pdbfiles = pdbfiles[:N_MAX]

smiles_strs = []
energies = []
gradients = []
xyzs = []
atomic_numbers = []
smiles_graphs = []
sequences = []

for p in tqdm(list(spicepaths), desc="Loading spice .npz and constructing graphs"):
    md = np.load(p, allow_pickle=True)
    smiles = md['smiles'].item()
    smiles_strs.append(smiles.decode())
    energies.append(md['energy_qm_total'])
    gradients.append(md['gradient_qm_total']) # shape num_confs, num_atoms, 3
    xyzs.append(md['xyz']) # shape num_confs, num_atoms, 3
    sequence = p.stem
    sequences.append(sequence)
    atomic_numbers.append(md['atomic_numbers']) # shape num_atoms

    # mol = MolData.load(p)
    # g = mol.to_dgl()
    # smiles_graphs.append(g) # only keep atom nodes

    # Construct openff molecule and DGL graph
    offmol = get_openff_molecule(smiles)
    mol = Molecule.from_openff_molecule(offmol, partial_charges=np.zeros(offmol.n_atoms))
    smiles_graphs.append(mol.to_dgl()) # only keep atom nodes

sequences[:10]
#%%

# --- b) Construct graphs from pdbs ---
from grappa.utils.openff_utils import mol_from_pdb
from grappa.data import Molecule

pdb_graphs = []
pdb_names = []
pdb_strings = []
for pdbfile in tqdm(pdbfiles, desc="Loading pdbs and constructing graphs"):
    with open(pdbfile) as f:
        pdbstr = f.read()
    offmol = mol_from_pdb(pdbstr)
    mol = Molecule.from_openff_molecule(offmol, partial_charges=np.zeros(offmol.n_atoms))
    pdb_graphs.append(mol.to_dgl())
    pdb_names.append(pdbfile.name)
    pdb_strings.append(pdbstr)

#%%

# --- c) Find isomorphisms and store/apply atom permutation ---
from grappa.utils.graph_utils import get_isomorphic_permutation, get_isomorphisms

isomorphisms = get_isomorphisms(pdb_graphs, smiles_graphs)
# %%
print(f'found {len(isomorphisms)} isomorphisms out of {len(pdb_graphs)} pdb graphs and {len(smiles_graphs)} smiles graphs')

# %%

spice_with_pdb = Path('spice-with-pdb')
spice_with_pdb.mkdir(exist_ok=True)

from rdkit import Chem

def _permute_smiles(orig_smiles, perm, atomic_numbers_permuted):
    """
    Given an original mapped SMILES, a permutation, and the permuted atomic numbers,
    return the mapped SMILES for the permuted atom order and assert atomic numbers match.
    """
    offmol = get_openff_molecule(orig_smiles)
    rdmol = offmol.to_rdkit()
    # Assign atom map numbers in original order (1-based)
    for k, a in enumerate(rdmol.GetAtoms(), start=1):
        a.SetAtomMapNum(k)
    assert len(perm) == rdmol.GetNumAtoms(), (len(perm), rdmol.GetNumAtoms())
    # Apply permutation
    rdmol_perm = Chem.RenumberAtoms(rdmol, [int(x) for x in perm])
    # Reassign map numbers to reflect new order
    for k, a in enumerate(rdmol_perm.GetAtoms(), start=1):
        a.SetAtomMapNum(k)
    mapped_smiles_permuted = Chem.MolToSmiles(
        rdmol_perm,
        canonical=False,
        allHsExplicit=True,
    )
    # Test: create openff molecule from mapped_smiles_permuted and check atomic numbers
    offmol_test = get_openff_molecule(mapped_smiles_permuted)
    atomic_numbers_test = np.array([a.atomic_number for a in offmol_test.atoms])
    assert np.all(atomic_numbers_test == atomic_numbers_permuted), (atomic_numbers_test, atomic_numbers_permuted)
    return mapped_smiles_permuted

#%%

# PERMUTE EACH MOLECULE AND SAVE:

for (i_pdb, i_smiles) in isomorphisms:
    pdb_graph = pdb_graphs[i_pdb]
    smiles_graph = smiles_graphs[i_smiles]
    perm = get_isomorphic_permutation(pdb_graph, smiles_graph)
    gradients_permuted = gradients[i_smiles][:, perm, :]
    xyzs_permuted = xyzs[i_smiles][:, perm, :]
    atomic_numbers_permuted = atomic_numbers[i_smiles][perm]

    # Permute SMILES and test
    mapped_smiles_permuted = _permute_smiles(smiles_strs[i_smiles], perm, atomic_numbers_permuted)

    # save to npz
    data = {
        "smiles_original_wrong_permutation": smiles_strs[i_smiles],
        "smiles": mapped_smiles_permuted,
        "sequence": sequences[i_smiles],
        "pdb": pdb_strings[i_pdb],
        "energy_qm_total": energies[i_smiles],
        "gradient_qm_total": gradients_permuted,
        "xyz": xyzs_permuted,
        "atomic_numbers": atomic_numbers_permuted,
    }

    outpath = spice_with_pdb / f"{sequences[i_smiles]}.npz"
    np.savez(outpath, **data)

# %%

# TEST:

# quick check of the permute smiles function:
# --- make mapped SMILES for the permuted atom order ---
from rdkit import Chem

(i_pdb, i_smiles) = list(isomorphisms)[0]

perm = get_isomorphic_permutation(pdb_graphs[i_pdb], smiles_graphs[i_smiles])

atomic_numbers_permuted = atomic_numbers[i_smiles][perm]

orig_smiles = smiles_strs[i_smiles]

mapped_smiles_permuted = _permute_smiles(orig_smiles, perm, atomic_numbers_permuted)

# print(f"Original SMILES:        {orig_smiles}")
# print(f"Mapped SMILES permuted: {mapped_smiles_permuted}")
# %%
# simple test: create openff molecule from mapped_smiles_permuted and check if the atomic numbers match the ones obtained by permuting as above:
offmol_old = get_openff_molecule(orig_smiles)
offmol_test = get_openff_molecule(mapped_smiles_permuted)
atomic_numbers_test = np.array([a.atomic_number for a in offmol_test.atoms])
assert np.all(atomic_numbers_test == atomic_numbers_permuted), (atomic_numbers_test, atomic_numbers_permuted)
# %%
