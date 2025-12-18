#!/bin/bash

# Create an empty variable to store sorted BAM files for featureCounts
sorted_files=""

# Loop through all .bam files in the directory
for file in *.bam
do
  # Extract the base name of the file (e.g., "1" from "1.bam")
  base=$(basename "$file" .bam)

  # Sort the BAM file and output to a new sorted BAM file
  samtools sort "$file" -o "${base}.sorted.bam"

  # Check if the sorting was successful before deleting the original file
  if [ $? -eq 0 ]; then
    # Add the sorted BAM file to the list
    sorted_files="${sorted_files} ${base}.sorted.bam"
    
    # Remove the original unsorted BAM file
    rm "$file"
  else
    echo "Sorting failed for $file. Original file not deleted."
  fi
done
