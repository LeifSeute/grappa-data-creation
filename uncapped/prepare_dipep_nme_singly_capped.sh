#! /bin/bash

# generates pdbs and smiles strings, then samples states

set -e

# load the sequences
############################################
# Read sequences into a bash array
# ACE IS LEFT, NME IS RIGHT, SO WE NEED TO FIX THE LEFT SIDE (this is the file with 'right' in the name)
# (we want every residue to have a missing ace cap a fixed number of times)
mapfile -t sequences < sequences/singly_capped_dipeptide_right.txt

# Convert the array into a space-separated string of sequences
sequences_str="${sequences[*]}"
############################################

THIS_DIR=/hits/basement/mbm/seutelf/grappa-data-creation/uncapped

source ~/.bashrc
conda activate pepgen

python generate_pdbs.py --folder $THIS_DIR/data/dipep_nme_singly_capped -s $sequences_str --nme_cap


echo Done. Now generating smiles strings...

conda activate espaloma

python smiles_string.py $THIS_DIR/data/dipep_nme_singly_capped



echo Done. Now generating states via MD...

conda activate grappa

python generate_states.py $THIS_DIR/data/dipep_nme_singly_capped -n 50 -t 300