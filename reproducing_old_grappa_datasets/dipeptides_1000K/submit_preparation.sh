I_MIN=0
I_MAX=9

mkdir -p logs

# submit n jobs:
for i in $(seq $I_MIN $I_MAX); do
    echo "Submitting job $i"
    sbatch gpu_job.sh bash prepare_dipeptides.sh $i
done