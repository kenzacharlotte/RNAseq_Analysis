#!/bin/bash

STAR --runThreadN 18 \
     --runMode genomeGenerate \
     --genomeDir star_ref_mm10_Refseq \
     --genomeFastaFiles mm10.fa \
     --sjdbGTFfile mm10.ncbiRefSeq.gtf \
     --sjdbOverhang 100
