# RNAseq_Analysis
Basic RNAseq workflow with the following organisation

project_name/
├── raw_data/             # tes fichiers bruts (.fastq.gz)
├── trimmed_data/         # fichiers après trimming (Trimmomatic ou autre)
├── fastqc_results/       # résultats FastQC bruts
│   ├── raw/              # FastQC sur les fichiers bruts
│   └── trimmed/          # FastQC sur les fichiers trimmed
├── multiqc_results/      # résultats MultiQC (rapports combinés)
├── logs/                 # fichiers log des différentes étapes
└── scripts/              # tes scripts bash / Python

A - PREPROCESSING #################################
###1 - Receive the fastq from a sequencing commpagny
Create a directory containing all the fastq and compare md5sum with bash code

###2 - Quality check #################################
Controle the quality sequencing with fastqc

###3 - Trimming (optionnal)

B - PROCESSING #################################
###1 - Alignment with STAR

###2 - Counting features with HTSeqcount

C - STATISTICAL ANALYSIS #################################
