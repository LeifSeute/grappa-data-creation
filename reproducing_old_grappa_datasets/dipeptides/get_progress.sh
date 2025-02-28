THIS_DIR=/hits/basement/mbm/seutelf/grappa-data-creation/dipeptides
DS=data/dipeptides_300K

pushd $THIS_DIR
python ../get_progress.py "$THIS_DIR/$DS"
popd