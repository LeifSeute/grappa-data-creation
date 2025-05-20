The 3bpa dataset is publicly available at the SI section at https://pubs.acs.org/doi/10.1021/acs.jctc.1c00647

In order to run the scripts in this directory, the zip file needs to be downloaded and, moved to this directory and unzipped into 'dataset_3BPA'.

Then, the gaff force field and grappa-1.4.0 can be evaluated by running
```bash
bash eval.sh
```
in a conda environment with openmm, openff.toolkit, and grappa installed.

For fine-tuning grappa on the 3bpa dataset, generate the two datasets used for training by running
```bash
python make_dataset.py
```