#%%
from pathlib import Path
import numpy as np
from grappa.data import MolData
from grappa.utils import openmm_utils, openff_utils
from tqdm import tqdm

def to_grappa_format(path, forcefield, forcefield_type='openmm', charge_model='amber99', target_dir=None, crmse_limit=15):
    """
    Converts the psi4_energies.npy and psi4_forces.npy files in the given folder to a format that can be read by grappa. 
    """
    assert path.is_dir(), "The given path is not a directory."

    energy = np.load(path/"psi4_energies.npy")
    gradient = -np.load(path/"psi4_forces.npy")

    xyz = np.load(path/"positions.npy")
    atomic_numbers = np.load(path/"atomic_numbers.npy") # not needed

    pdbstring = open(path/"pep.pdb").read()

    sequence = path.stem

    ###### some special cases:
    if 'uncapped' in str(target_dir):
        # usually capping is implied, for uncapped, it is not...
        # we use _ for signalling a missing cap:
        if sequence.endswith("Z"): # capped
            sequence = sequence[:-1]
        else: # uncapped
            sequence = sequence+"_"
        if sequence.startswith("B"): # capped
            sequence = "_"+sequence[1:]
        else: # uncapped
            sequence = "_"+sequence
    ######

    smiles = str(np.load(path/"smiles.npy"))
    assert smiles is not None and smiles != "None", f"smiles is None for {path.stem}"

    mapped_smiles = str(np.load(path/"mapped_smiles.npy"))

    valid_idxs = np.isfinite(energy)
    valid_idxs = np.where(valid_idxs)[0]
    energy = energy[valid_idxs]
    gradient = gradient[valid_idxs]
    xyz = xyz[valid_idxs]
    mol_id = smiles

    ref_forcefield = forcefield

    if forcefield_type == 'openmm':
        # get topology:
        topology = openmm_utils.topology_from_pdb(pdbstring)
        ff = openmm_utils.get_openmm_forcefield(ref_forcefield)
        system = ff.createSystem(topology)

    elif forcefield_type == 'openff' or forcefield_type == 'openmmforcefields': 
        system, topology, _ = openff_utils.get_openmm_system(mapped_smiles=mapped_smiles, openff_forcefield=ref_forcefield)
    else:
        raise ValueError(f"forcefield_type must be either openmm, openff or openmmforcefields but is {forcefield_type}")

    # create moldata object from the system (calculate the parameters, nonbonded forces and create reference energies and gradients from that)
    moldata = MolData.from_openmm_system(openmm_system=system, openmm_topology=topology, xyz=xyz, gradient=gradient, energy=energy, mol_id=mol_id, pdb=pdbstring, sequence=sequence, allow_nan_params=True, charge_model=charge_model, ff_name=ref_forcefield, smiles=smiles, mapped_smiles=mapped_smiles)

    openmm_gradients = moldata.ff_gradient[ref_forcefield]["total"]

    # calculate the crmse:
    crmse = np.sqrt(np.mean((openmm_gradients-gradient)**2))

    if crmse > crmse_limit:
        print(f"Warning: crmse between {forcefield} and QM is {round(crmse, 3)} for {path.stem}, std of gradients is {round(np.std(gradient), 1)}")

    moldata.save(target_dir/(path.stem+'.npz'))

    return


def convert_dataset(path, forcefield, forcefield_type='openmm', charge_model='classical', target_dir=None, crmse_limit=15, skip_errs=False):
    """
    Converts the psi4_energies.npy and psi4_forces.npy files in the given folder to a format that can be read by grappa. 
    """
    assert path.is_dir(), "The given path is not a directory."

    if target_dir is None:
        target_dir = Path(path).parent/("grappa_"+Path(path).name+f"_{forcefield}")
    else:
        target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nConverting dataset {path}\nto grappa format {target_dir}\n")

    paths = list(p for p in path.iterdir() if p.is_dir())
    for i, subdir in enumerate(tqdm(paths, desc="Converting")):
        try:
            if subdir.is_dir():
                to_grappa_format(subdir, forcefield, forcefield_type, charge_model, target_dir, crmse_limit=crmse_limit)
        except Exception as e:
            print()
            if skip_errs:
                print(f"Error in {subdir}: {e}")
                continue
            else:
                raise e

    print()

    return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source_path", type=str, help="Path to the folder with subfolders that contain psi4_energies.npy and psi4_forces.npy files."
    )
    parser.add_argument(
        "--target_dir", type=str, help="Path to the target folder in which the dataset is stored as collection of npz files.", default=None
    )
    parser.add_argument(
        "--forcefield","-ff", type=str, default="amber99sbildn", help="Forcefield to use for the conversion."
    )
    parser.add_argument(
        "--forcefield_type", type=str, default="openmm", help="Type of forcefield to use for the conversion. Available: openmm, openff, openmmforcefields."
    )
    parser.add_argument(
        "--charge_model", "-cm", type=str, default="amber99", help="Charge model of the underlying forcefield. Available: amber99, charmm, am1BCC"
    )
    parser.add_argument(
        "--crmse_limit", "-crmse", type=float, default=15, help="Print warnings if the crmse between the forcefield and QM is higher than this limit. Unit: kcal/mol/Angstrom."
    )
    parser.add_argument(
        "--skip_errs", action="store_true", help="Skip errors during conversion."
    )
    args = parser.parse_args()

    source_path = Path(args.source_path)
    target_dir = Path(args.target_dir) if args.target_dir is not None else None
    forcefield = args.forcefield
    forcefield_type = args.forcefield_type
    charge_model = args.charge_model

    convert_dataset(source_path, forcefield, forcefield_type, charge_model, target_dir, crmse_limit=args.crmse_limit, skip_errs=args.skip_errs)
