SUFFIX=${1:-""}
N=10

mkdir -p logfiles

for i in {1..10}
do
	sbatch dft_agent$SUFFIX.sh
done
