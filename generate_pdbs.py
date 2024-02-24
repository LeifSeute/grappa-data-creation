import argparse
from pepgen.pepgen import generate_peptide
from pathlib import Path
import os

if __name__ == "__main__":
    # Initialize argparse
    parser = argparse.ArgumentParser(description='Pass single letter sequences to create peptide pdb files.')
    parser.add_argument('--folder', type=str, help='Output folder for the data.', required=True)
    parser.add_argument('--sequence', '-s', type=str, nargs='+', help='Single letter sequence.', default=[])
    parser.add_argument('--nme_cap', action='store_true', help='Add an NME cap to the peptide. (Just adds a Z to the right of the sequence)', default=False)
    parser.add_argument('--ace_cap', action='store_true', help='Add an ACE cap to the peptide. (Just adds a B to the left of the sequence)', default=False)
    
    args = parser.parse_args()

    print(f"Generating dataset with {args.sequence}.")

    for s in args.sequence:

        print(s)

        os.makedirs(args.folder, exist_ok=True)
        generate_peptide(code=s, dir=Path(args.folder)/s, e=False, silent=True, overwrite=True, openmm_standard=False, t=False, nme_cap=args.nme_cap, ace_cap=args.ace_cap)