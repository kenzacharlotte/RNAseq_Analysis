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
- creates html report for each fastq with
For more informations about fastqc check https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

##### *Trimming (optionnal)*

## B - PROCESSING 
##### *Alignment - STAR*
1 - Build STAR
- gtf file needed
```py
build_STAR_index.sh
```
2 - Proceed to the alignment
```py
STAR_alignment.sh
```
##### *Counting features - HTSeqcount*

## C - Statistical analysis 
