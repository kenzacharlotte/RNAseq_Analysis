# RNAseq_Analysis
Basic rRNAseq workflow 

A - PREPROCESSING
#1 - Receive the fastq from a sequencing commpagny
Create a directory containing all the fastq and compare md5sum with bash code

#2 - Quality check 
Controle the quality sequencing with fastqc

#3 - Trimming (optionnal)

B - PROCESSING
#1 - Alignment with STAR

#2 - Counting features with HTSeqcount

C - STATISTICAL ANALYSIS
