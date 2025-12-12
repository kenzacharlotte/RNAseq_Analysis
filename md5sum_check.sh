#!/bin/bash

# 1. Creation log file
echo "Creating md5sum md5sum_copied" >> logs/md5sum.out

# 2. Creation of md5sum of the copied fatsq
md5sum fastqc_results/raw/*fq.gz > fastqc_results/raw/md5sum/md5sum_copied.txt
echo "md5sum_copied successfully created" >> logs/md5sum.out
# 3. Counting identic files
echo "Number of similar files:" >> logs/md5sum.out
comm -12 <(awk '{print $1}' fastqc_results/raw/md5sum/md5sum_copied.txt | sort) \
         <(awk '{print $1}' fastqc_results/raw/md5sum/md5sum_check.txt | sort) \
| wc -l >> logs/md5sum.out

# 4. Listing of the non-similar files
echo "If error, non similar files:" >> md5sum.out
comm -23 <(awk '{print $1}' fastqc_results/raw/md5sum/md5sum_copied.txt | sort) \
         <(awk '{print $1}' fastqc_results/raw/md5sum/md5sum_check.txt | sort) \
| while read md5; do
    grep "$md5" fastqc_results/raw/md5sum/md5sum_copied.txt >> logs/md5sum.out
done
