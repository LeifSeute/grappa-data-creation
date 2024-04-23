#! /bin/bash

# generates pdbs and smiles strings, then samples states

set -e

THIS_DIR=/hits/basement/mbm/seutelf/grappa-data-creation/uncapped

DATADIR=$THIS_DIR/data/dipeptides_300K

# load the sequences
############################################
# Read sequences into a bash array
# ACE IS LEFT, NME IS RIGHT, SO WE NEED TO FIX THE LEFT SIDE (this is the file with 'right' in the name)
# (we want every residue to have a missing ace cap a fixed number of times)
mapfile -t sequences < $THIS_DIR/sequences/dipeptides.txt

# Convert the array into a space-separated string of sequences
sequences_str="${sequences[*]}"
############################################

source ~/.bashrc
conda activate pepgen # requires pepgen, numpy

python ../generate_pdbs.py --folder $DATADIR -s $sequences_str --nme_cap --ace_cap

echo Done. Now generating smiles strings...


conda activate espaloma # requires openff, numpy

python ../smiles_string.py $DATADIR


echo Done. Now generating states via MD...

conda activate grappa # requires grappa and openmm

# python ../generate_states.py $DATADIR -n 50 -t 300
python ../generate_states.py $DATADIR -n 5 -t 300 -p