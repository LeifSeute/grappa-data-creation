import argparse
import shutil
from pathlib import Path
import random

def split_folder(source_folder: Path, split_ratio: float = 0.3):
    """
    Randomly splits the directories within source_folder into two subfolders,
    with a specified ratio going to  _fast folder, and the rest to a _basement folder.
    
    Args:
    - source_folder: Path to the source folder.
    - split_ratio: Fraction of directories to allocate to the _fast subfolder.
    """

    # Ensure source_folder is a Path object
    source_folder = Path(source_folder)
    
    # Define target subfolders
    fast_folder = source_folder.parent / (source_folder.name + "_fast")
    basement_folder = source_folder.parent / (source_folder.name + "_basement")

    # Create target subfolders if they do not exist
    fast_folder.mkdir(exist_ok=True, parents=True)
    basement_folder.mkdir(exist_ok=True, parents=True)

    # List all directories in the source_folder
    dirs = [d for d in source_folder.iterdir() if d.is_dir()]

    # Determine split index
    fast_count = int(len(dirs) * split_ratio)
    random.shuffle(dirs)  # Shuffle the directories to randomly select them
    
    # Split directories
    fast_dirs = dirs[:fast_count]
    basement_dirs = dirs[fast_count:]

    # Copy directories to their new locations
    for d in fast_dirs:
        shutil.copytree(d, fast_folder / d.name)
    for d in basement_dirs:
        shutil.copytree(d, basement_folder / d.name)
    
    print(f"Split {len(dirs)} directories into {len(fast_dirs)} in '_fast' and {len(basement_dirs)} in '_basement'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Randomly splits a folder into '_basement' and '_fast' subfolders.")
    parser.add_argument("folder", type=str, help="Path to the source folder.")
    parser.add_argument("--split", type=float, default=0.3, help="Split ratio for the '_fast' folder. Default is 0.3.")
    
    args = parser.parse_args()
    
    split_folder(args.folder, args.split)
