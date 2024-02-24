from pathlib import Path
import sys

# accept parent directory as command line argument
if len(sys.argv) > 1:
    parent_dir = Path(sys.argv[1])
else:
    raise ValueError("Please provide the parent directory as command line argument.")

# Find all subdirectories in the parent directory
subdirs = [item for item in parent_dir.iterdir() if item.is_dir()]

# Extract suffixes (integers after "_")
suffixes = set()
for subdir in subdirs:
    subdir_name = subdir.name
    try:
        suffix = int(subdir_name)
        suffixes.add(suffix)
    except:
        continue
    

# Find the smallest non-occurring positive integer
for i in range(1, len(suffixes) + 2):
    if i not in suffixes:
        print(i)
        break