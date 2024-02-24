#! /bin/bash

# generates pdbs and smiles strings, then samples states

set -e

source ~/.bashrc
conda activate pepgen

python generate_pdbs.py --folder data/uncapped_test_nme -s G A F --nme_cap
python generate_pdbs.py --folder data/uncapped_test_ace -s G --ace_cap

echo Done. Now generating states via MD...

conda activate grappa

python generate_states.py data/uncapped_test_nme -n 10 -t 300
python generate_states.py data/uncapped_test_ace -n 10 -t 300

echo Done. Now generating smiles strings...

conda activate espaloma

python smiles_string.py data/uncapped_test_nme
python smiles_string.py data/uncapped_test_ace