#!/bin/bash

MAX_IDX=9

NAME=dipeptides_1000K

# rename the target directory to data/${NAME}
mv data/0_${NAME} data/${NAME}

set -e

# Loop through all directories from 1 to MAX_IDX
for i in $(seq 1 $MAX_IDX); do
  # Check if the directory exists
  if [ -d "data/${i}_${NAME}" ]; then
    # Move all subdirectories to data/0_${NAME}/
    mv data/${i}_${NAME}/* data/${NAME}/
  fi
done

# Find and delete all log.txt files in the target directory
find data/${NAME}/ -type f -name "log.txt" -delete
