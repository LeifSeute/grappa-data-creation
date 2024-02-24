import argparse
import random
from pathlib import Path

# Natural amino acids:
AALETTERS = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
SEED = 42
random.seed(SEED)

this_file_path = Path(__file__).parent.absolute()

def write_all_AAs():
    with open(this_file_path/"AAs.txt", "w") as file:
        for aa in AALETTERS:
            file.write(aa + "\n")
    print(f"Written {len(AALETTERS)} amino acids to {this_file_path/'AAs.txt'}.")

def sample_dipeptides(partner_position='left', n_partners=2):
    filename = f"singly_capped_dipeptide_{partner_position}.txt"
    with open(this_file_path/filename, "w") as file:
        for aa in AALETTERS:
            partners = random.sample(AALETTERS, n_partners)
            for partner in partners:
                if partner_position == 'left':
                    sequence = partner + aa
                else:
                    sequence = aa + partner
                file.write(sequence + "\n")

    print(f"Written {len(AALETTERS)} amino acids to {this_file_path/filename}.")

if __name__ == "__main__":
    write_all_AAs()
    sample_dipeptides('left', 2)
    sample_dipeptides('right', 2)
