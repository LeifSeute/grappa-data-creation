#!/bin/bash

# call with `sbatch dft_agent_basement_large.sh`

#SBATCH -N 1                    
#SBATCH -p haswell.p
#SBATCH --job-name=psi4      # Job name
#SBATCH --mem=58G                       # Total memory for all tasks
#SBATCH --mincpus=16
#SBATCH -t 24:00:00
#SBATCH --output=./logs/job.o%j
#SBATCH --error=./logs/job.e%j

THIS_DIR=/hits/basement/mbm/seutelf/grappa-data-creation/dipeptides_1000K
DS=data/dipeptides_1000K

MEM=64                                   # Memory per python script
CORES=16                                 # Cores per python script
NUM_AGENTS=1                             # Number of parallel agents per node

pids=() # Array to hold process IDs

source /hits/basement/mbm/seutelf/.bashrc_user

conda activate psi4

cd $THIS_DIR

for i in $(seq 1 $NUM_AGENTS); do
    nohup python ../single_points.py "$THIS_DIR/$DS" -m $MEM -t $CORES -s &
    pids+=($!) # Capture the process ID of the background job
done

# Wait for all background processes to complete
for pid in "${pids[@]}"; do
    wait $pid
done
