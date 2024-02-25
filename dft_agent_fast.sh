#!/bin/bash


#SBATCH -N 1                    
#SBATCH -p cascade.p
#SBATCH --job-name=psi4      # Job name
#SBATCH --mem=90G                       # Total memory for all tasks
#SBATCH --mincpus=38
#SBATCH -t 24:00:00
#SBATCH --output=logfiles/job.o%j
#SBATCH --error=logfiles/job.e%j

THIS_DIR=/hits/fast/mbm/seutelf/software/grappa_data_creation
DS=data/uncapped_dataset_fast

MEM=32                                   # Memory per python script
CORES=12                                 # Cores per python script
NUM_AGENTS=3                             # Number of parallel agents

pids=() # Array to hold process IDs

source /hits/fast/mbm/seutelf/.bashrc_user

conda activate psi4

cd $THIS_DIR

for i in $(seq 1 $NUM_AGENTS); do
    nohup python single_points.py "$THIS_DIR/$DS" -m $MEM -t $CORES -s &
    pids+=($!) # Capture the process ID of the background job
done

# Wait for all background processes to complete
for pid in "${pids[@]}"; do
    wait $pid
done
