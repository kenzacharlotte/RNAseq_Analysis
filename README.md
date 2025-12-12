# RNAseq_Analysis
*Basic RNAseq analysis workflow*


## A - PREPROCESSING 
##### *Organize your files*

```py
create_repository.sh
```

project_name/  
📁 raw_data/             – fichiers bruts (.fastq.gz)  
📁 trimmed_data/         – fichiers après trimming  
📁 fastqc_results/       – résultats FastQC  
   📁 raw/               – FastQC sur fichiers bruts  
       📁 md5sum/        – md5sum files check and copied  
   📁 trimmed/           – FastQC sur fichiers trimmed  
📁 multiqc_results/      – résultats MultiQC  
📁 logs/                 – fichiers log  
📁 scripts/              – scripts bash / Python  

```py
md5sum_check.sh
```

###1 - Receive the fastq from a sequencing commpagny
Create a directory containing all the fastq and compare md5sum with bash code

###2 - Quality check #################################
Controle the quality sequencing with fastqc

###3 - Trimming (optionnal)

B - PROCESSING #################################
###1 - Alignment with STAR

###2 - Counting features with HTSeqcount

C - STATISTICAL ANALYSIS #################################
