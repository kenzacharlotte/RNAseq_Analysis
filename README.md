# RNAseq_Analysis
*Basic RNAseq analysis workflow*

## A - PREPROCESSING 
##### *Structure*

```py
create_repository.sh
```
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
```
##### *Md5sum*
```py
md5sum_check.sh
```
- creates md5sum of the copied fastq
- check if the md5sum are matching in `md5sum.out`
- writes in `md5sum.out` the files that are not matching

##### *Quality check*
```py
qualitycheck_fastqc.sh
```
- run fastqc for all raw data (*fq.gz)
- creates html files 

###3 - Trimming (optionnal)

B - PROCESSING #################################
###1 - Alignment with STAR

###2 - Counting features with HTSeqcount

C - STATISTICAL ANALYSIS #################################
