#%%
import h5py

spicepath = "SPICE-2.0.1.hdf5"

# %%


import os
import numpy as np

dsname = "SPICE Dipeptides Single Points Dataset v1.3"
output_dir = "spice-npz"
os.makedirs(output_dir, exist_ok=True)

from tqdm.auto import tqdm
import random

with h5py.File(spicepath, "r") as f:
	group_names = list(f.keys())
	# permute:
	random.shuffle(group_names)
	for group_name in tqdm(group_names, desc="Processing SPICE groups"):
		group = f[group_name]
		subset = group.get("subset")
		subset_name = subset[()][0].decode()
		if subset_name != dsname:
			continue
		# Extract arrays as in the example
		smiles = group["smiles"][()]
		if isinstance(smiles, (bytes, np.bytes_)):
			smiles = smiles.decode()
		atomic_numbers = group["atomic_numbers"][()]
		conformations = group["conformations"][()]
		formation_energy = group["formation_energy"][()]
		dft_total_energy = group["dft_total_energy"][()]
		dft_total_gradient = group["dft_total_gradient"][()]
		# Save as npz
		outpath = os.path.join(output_dir, f"{group_name}.npz")
		np.savez(
			outpath,
			smiles=smiles,
			atomic_numbers=atomic_numbers,
			xyz=conformations,
			formation_energy=formation_energy,
			energy_qm_total=dft_total_energy,
			gradient_qm_total=dft_total_gradient,
		)
		print(f"Saved {outpath}")