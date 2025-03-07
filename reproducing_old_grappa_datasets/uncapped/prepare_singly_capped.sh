#! /bin/bash

# generates pdbs and smiles strings, then samples states

set -e

# load the sequences
############################################
# Read sequences into a bash array
mapfile -t sequences < sequences/AAs.txt

# Convert the array into a space-separated string of sequences
sequences_str="${sequences[*]}"
############################################

THIS_DIR=/hits/basement/mbm/seutelf/grappa-data-creation/uncapped

source ~/.bashrc
conda activate pepgen

python generate_pdbs.py --folder $THIS_DIR/data/singly_capped_nme -s $sequences_str --nme_cap
python generate_pdbs.py --folder $THIS_DIR/data/singly_capped_ace -s $sequences_str --ace_cap




echo Done. Now generating smiles strings...

conda activate openff

python smiles_string.py $THIS_DIR/data/singly_capped_nme
python smiles_string.py $THIS_DIR/data/singly_capped_ace



echo Done. Now generating states via MD...

conda activate grappa_openmm

python generate_states.py $THIS_DIR/data/singly_capped_nme -n 50 -t 300
python generate_states.py $THIS_DIR/data/singly_capped_ace -n 50 -t 300