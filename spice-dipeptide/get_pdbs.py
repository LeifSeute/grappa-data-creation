#%%
from grappa.data import MolData
from grappa.utils.data_utils import get_moldata_path 

spicepaths = list(get_moldata_path('spice-dipeptide-amber99').glob('*.npz'))
len(spicepaths)
# %%
from tqdm.auto import tqdm
import numpy as np
from pathlib import Path

for p in tqdm(list(spicepaths)):
    md = np.load(p, allow_pickle=True)
    pdb = md['pdb'].item()
    # save as p.name.pdb:
    with open(Path('spice-pdbs')/(p.stem+'.pdb'), 'w') as f:
        f.write(pdb)
# %%
