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
    smiles_strs.append(smiles)
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
# now iterate over isomorphisms, calculate permutation that needs to be applied to gradients/xyz/atomic_numbers

spice_with_pdb = Path('spice-with-pdb')

spice_with_pdb.mkdir(exist_ok=True)

for (i_pdb, i_smiles) in isomorphisms:
    pdb_graph = pdb_graphs[i_pdb]
    smiles_graph = smiles_graphs[i_smiles]
    perm = get_isomorphic_permutation(pdb_graph, smiles_graph)
    gradients_permuted = gradients[i_smiles][:, perm, :]
    xyzs_permuted = xyzs[i_smiles][:, perm, :]
    atomic_numbers_permuted = atomic_numbers[i_smiles][perm]

    # save to npz
    data = {
        "smiles_wrong_permutation": smiles_strs[i_smiles],
        "sequence": sequences[i_smiles],
        "pdb": pdb_strings[i_pdb],
        "energy_qm_total": energies[i_smiles],
        "gradient_qm_total": gradients_permuted,
        "xyz": xyzs_permuted,
        "atomic_numbers": atomic_numbers_permuted,
    }

    # save to npz:
    outpath = spice_with_pdb / f"{sequences[i_smiles]}.npz"
    np.savez(outpath, **data)

# %%
