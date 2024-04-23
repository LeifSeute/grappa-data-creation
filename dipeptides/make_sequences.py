import argparse
import random
from pathlib import Path

# Natural amino acids:
AALETTERS = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y', '1', '2', '3', '4', '5', '6'] # protonations states: 1: HID, 2: HIP, 3: GLH, 4: ASH, 5: LYN 6: CYX

SEED = 42
random.seed(SEED)

this_file_path = Path(__file__).parent.absolute()

def write_all_AAs():
    with open(this_file_path/"sequences"/"AAs.txt", "w") as file:
        for aa in AALETTERS:
            file.write(aa + "\n")
    print(f"Written {len(AALETTERS)} amino acids to {this_file_path/'AAs.txt'}.")

# AB != BA because of the caps!
def sample_dipeptides():
    filename = f"dipeptides.txt"
    with open(this_file_path/"sequences"/filename, "w") as file:
        for aa in AALETTERS:
            for aa2 in AALETTERS:
                file.write(aa + aa2 + "\n")

    print(f"Written {len(AALETTERS)} amino acids to {this_file_path/filename}.")

if __name__ == "__main__":
    (this_file_path/"sequences").mkdir(parents=True, exist_ok=True)
    write_all_AAs()
    sample_dipeptides()
