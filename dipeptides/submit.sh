SUFFIX=${1:-""}

mkdir -p logfiles

N=10
for i in {1..10}
do
	sbatch dft_agent$SUFFIX.sh
done
