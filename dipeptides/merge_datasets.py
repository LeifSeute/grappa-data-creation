'''
Create a merged dataset from the ten subdatasets.
'''

from pathlib import Path
import shutil
import os

if __name__ == '__main__':
    # Define the paths to the ten subdatasets:
    this_dir = Path(__file__).parent
    new_dir = this_dir / 'data' / 'dipeptides_300K'

    new_dir.mkdir(exist_ok=False)

    N_SUBDATASETS = 10

    for i in range(N_SUBDATASETS):
        # move all children from data/i_dipeptides_300K to data/dipeptides_300K:
        subdataset_dir = this_dir / f'data/{i}_dipeptides_300K'
        # remove the log.txt file:
        log_file = subdataset_dir / 'log.txt'
        if log_file.exists():
            os.remove(log_file)
        for child in subdataset_dir.iterdir():
            shutil.move(str(child), str(new_dir))

        # remove the empty directory:
        os.rmdir(subdataset_dir)
