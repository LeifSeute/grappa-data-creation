#! /bin/bash

# generates pdbs and smiles strings, then samples states

set -e

# load the sequences
############################################
# Read sequences into a bash array
# ACE IS LEFT, NME IS RIGHT, SO WE NEED TO FIX THE RIGHT SIDE (this is the file with 'left' in the name)
# (we want every residue to have a missing nme cap a fixed number of times)
mapfile -t sequences < sequences/singly_capped_dipeptide_left.txt

# Convert the array into a space-separated string of sequences
sequences_str="${sequences[*]}"
############################################

THIS_DIR=/hits/basement/mbm/seutelf/grappa-data-creation/uncapped

source ~/.bashrc
conda activate pepgen

python ../generate_pdbs.py --folder $THIS_DIR/data/dipep_ace_singly_capped -s $sequences_str --ace_cap


echo Done. Now generating smiles strings...

conda activate openff

python ../smiles_string.py $THIS_DIR/data/dipep_ace_singly_capped



echo Done. Now generating states via MD...

conda activate grappa_openmm

python ../generate_states.py $THIS_DIR/data/dipep_ace_singly_capped -n 50 -t 300