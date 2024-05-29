#! /bin/bash

IDX=${1:-0}

echo "Running for index $IDX"

# generates pdbs and smiles strings, then samples states
THIS_DIR=$(dirname "$(realpath "$0")")
DATADIR=$THIS_DIR/data/${IDX}_dipeptides_500K

# load the sequences
############################################
# Read sequences into a bash array
mapfile -t sequences < $THIS_DIR/sequences/dipeptides_${IDX}.txt

# Convert the array into a space-separated string of sequences
sequences_str="${sequences[*]}"

echo "Generating data for sequences: $sequences_str"
############################################

source ~/.bashrc


conda activate pepgen # requires pepgen, numpy

python ../generate_pdbs.py --folder $DATADIR -s $sequences_str --nme_cap --ace_cap

echo Done. Now generating smiles strings...

conda activate openff # requires openff, numpy

python ../smiles_string.py $DATADIR


echo Done. Now generating states via MD...

conda activate grappa_openmm # requires grappa and openmm

python ../generate_states.py $DATADIR -n 30 -t 500 --t_max 1000
