#!/bin/bash

# 1. Log de création
echo "Creating md5sum md5sum_copied" >> md5sum.out

# 2. Créer le fichier md5sum
md5sum fastqc_results/raw/*fq.gz > fastqc_results/raw/md5sum_copied.txt
echo "md5sum_copied successfully created" >> md5sum.out

# 3. Compter les fichiers identiques
echo "Number of similar files:" >> md5sum.out
comm -12 <(awk '{print $1}' fastqc_results/raw/md5sum_copied.txt | sort) \
         <(awk '{print $1}' fastqc_results/raw/md5sum_check.txt | sort) \
| wc -l >> md5sum.out

# 4. Lister les fichiers non similaires
echo "If error, non similar files:" >> md5sum.out
comm -23 <(awk '{print $1}' fastqc_results/raw/md5sum_copied.txt | sort) \
         <(awk '{print $1}' fastqc_results/raw/md5sum_check.txt | sort) \
| while read md5; do
    grep "$md5" fastqc_results/raw/md5sum_copied.txt >> md5sum.out
done
