#!/bin/bash

# Create md5sum files 
md5sum fastqc_results/raw/*fq.gz > fastqc_results/raw/md5sum_copied.txt

# Print the md5sum that are identical 
comm -12 <(awk '{print $1}' md5sum_copied.txt | sort) <(awk '{print $1}' info/md5sum_check.txt | sort) 

# Count the number of lines (i.e. files) that are identical : this should match with the number of fastq you have in your repository
comm -12 <(awk '{print $1}' md5sum_copied.txt | sort) <(awk '{print $1}' info/md5sum_check.txt | sort) | wc -l 

# The arg -23 print the non matching md5sum 
# Non-matching md5sum are sent in the loop and "grepped" in the file to print the problematic ones
# Creates error_md5sum.out a file where all the non-matching mdf5sum and the associated files are written
# NB : grep strands for global regular expression print
comm -23 <(awk '{print $1}' md5sum_copied | sort) \
         <(awk '{print $1}' info/md5sum_check.txt | sort) \
| while read md5; do
    grep "$md5" md5sum_copied >> error_md5sum.out
  done

