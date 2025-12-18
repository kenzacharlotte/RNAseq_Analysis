# RNAseq_Analysis
*Basic RNAseq analysis workflow*

## A - PREPROCESSING 
##### *Structure*

```py
create_repository.sh
```
```
project_name/
📁 1.raw_data/             – raw data (.fq.gz)
   📁 md5sum/        – md5sum files check and copied
   📁 trimmed_data/         – fastq after trimming -optionnal (.fq.gz)
📁 2.fastqc_results/       – fastQC results (html)
   📁 trimmed/           – FastQC sur fichiers trimmed  
📁 3.star/
   📁 bam/               – fichiers log                 
📁 4.featurecounts/        – featurecounts results (txt)  
📁 logs/                 – fichiers log
📁 scripts/              – scripts bash 
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
##### *Counting features - featurecounts*
1 - Sort bam

```py
sort_bam.sh
```
```py
featurecounts.sh
```
## C - R analysis with Deseq2
- Exploring Dataset
- Create Deseq object
- Visualization
- DEgs identification
- Pathways identification
