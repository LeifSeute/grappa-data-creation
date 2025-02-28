#!/bin/bash

MAX_IDX=9

# rename the target directory to data/dipeptides_500K
mv data/0_dipeptides_500K data/dipeptides_500K

set -e

# Loop through all directories from 1 to MAX_IDX
for i in $(seq 1 $MAX_IDX); do
  # Check if the directory exists
  if [ -d "data/${i}_dipeptides_500K" ]; then
    # Move all subdirectories to data/0_dipeptides_500K/
    mv data/${i}_dipeptides_500K/* data/dipeptides_500K/
  fi
done

# Find and delete all log.txt files in the target directory
find data/dipeptides_500K/ -type f -name "log.txt" -delete
