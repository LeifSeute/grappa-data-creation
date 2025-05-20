#%%
import numpy as np
import re
import json

def parse_data(file_path, num_confs=None):
    data_blocks = []
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Split content into blocks based on occurrences of 'Lattice'
    blocks = re.split(r"(?=Lattice)", content)

    num_confs_ = 0

    for block in blocks:
        if not block.strip():
            continue
        if not 'Lattice' in block:
            continue
        # print(block)
        # Extract lattice parameters
        if num_confs is not None and num_confs_ >= num_confs:
            break
        num_confs_ += 1

        lattice = re.search(r'Lattice="([^"]+)"', block).group(1)
        lattice = list(map(float, lattice.split()))
        
        # Extract energy
        energy = float(re.search(r'energy=(-?\d+\.\d+)', block).group(1))
        
        # Extract dihedrals as JSON
        dihedrals = re.search(r'dihedrals="_JSON ([^"]+)"', block)
        if dihedrals:
            dihedrals = json.loads(dihedrals.group(1))
        
        # Extract periodic boundary conditions (pbc)
        pbc = re.search(r'pbc="([^"]+)"', block).group(1)
        pbc = pbc.split()
        
        # Extract atomic data (species, positions, forces)
        atoms = []
        atom_lines = block.strip().split('\n')[1:]
        for line in atom_lines:
            parts = line.split()
            if len(parts) == 7:
                species = parts[0]
                pos = list(map(float, parts[1:4]))
                forces = list(map(float, parts[4:7]))
                atoms.append({
                    'species': species,
                    'position': pos,
                    'forces': forces
                })
        
        # Store parsed block
        data_blocks.append({
            'lattice': lattice,
            'energy': energy,
            'dihedrals': dihedrals,
            'pbc': pbc,
            'atoms': atoms
        })
    
    return data_blocks

def get_dihedrals(data):
    '''
    Returns dihedrals as alpha, beta, gamma arrays.
    '''
    dihedrals = []
    for block in data:
        dihedrals.append(block['dihedrals'])
    return (np.array(dihedrals)[:,i] for i in range(3))

def to_arrays(data):
    '''
    Returns atomic numbers, positions, forces, and energies as numpy arrays with units Angstrom and kcal/mol.
    '''
    assert len(data) > 0, 'Data is empty'

    ELEMENTS = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9}
    
    atomic_numbers = [ELEMENTS[atom['species']] for atom in data[0]['atoms']]
    
    positions = []
    forces = []
    energies = []

    for block in data:
        assert [ELEMENTS[atom['species']] for atom in block['atoms']] == atomic_numbers, 'Atomic numbers do not match'
        energies.append(block['energy'])
        positions.append([atom['position'] for atom in block['atoms']])
        forces.append([atom['forces'] for atom in block['atoms']])


    # center energies:
    energies -= np.mean(energies)

    # convert from eV to kcal/mol:
    energies = np.array(energies) * 23.0605

    forces = np.array(forces) * 23.0605

    return np.array(atomic_numbers), np.array(positions), np.array(forces), np.array(energies)