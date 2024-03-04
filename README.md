## Grappa Data Creation

Scripts used to create the peptide datasets 'uncapped_amber99sbildn', 'tripeptides_amber99sbildn', 'hyp-dop_amber99sbildn' and 'dipeptide_rad' from [grappa-1.1](https://github.com/hits-mbm-dev/grappa/releases/tag/v.1.1.0).

For all of these, the workflow was:
- Generate PDB files from sequences using the [pepgen](https://github.com/hits-mbm-dev/pepgen) package (generate_pdbs.py)
- Sample states from MD at 300K using a traditional forcefield (amber99sbildn for peptides, an early version of grappa for radicals). To reach a larger conformational space, we simulate at 500K and reduce the temperature back to 300K before sampling a state (generate_states.py)
- Calculate single point energy and gradient with DFT (single_points.py)
- Save energy, gradient, sequence, pdb-string and atomic numbers as npz file that can be read by grappas data.MolData class (to_grappa_format.py)
