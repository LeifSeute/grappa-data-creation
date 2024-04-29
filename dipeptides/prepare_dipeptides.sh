#! /bin/bash

IDX=${1:-0}

echo "Running for index $IDX"

# generates pdbs and smiles strings, then samples states

THIS_DIR=/hits/fast/mbm/seutelf/software/grappa_data_creation/dipeptides
DATADIR=$THIS_DIR/data/${IDX}_dipeptides_300K

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


conda activate espaloma # requires openff, numpy

python ../smiles_string.py $DATADIR


echo Done. Now generating states via MD...

conda activate grappa # requires grappa and openmm

# python ../generate_states.py $DATADIR -n 5 -t 300 --t_max 1000 -p # test run with plotting
python ../generate_states.py $DATADIR -n 50 -t 300 --t_max 1000
