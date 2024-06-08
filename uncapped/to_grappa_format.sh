DS=data/uncapped_300K
THISDIR=$(dirname "$(realpath "$0")")
python $THISDIR/../to_grappa_format.py $THISDIR/$DS -ff amber99sbildn -cm amber99
#python $THISDIR/../to_grappa_format.py $THISDIR/$DS -ff charmm36 -cm charmm
python $THISDIR/../to_grappa_format.py $THISDIR/$DS -ff openff_unconstrained-1.2.0.offxml -cm am1BCC --forcefield_type openff
