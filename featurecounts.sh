#!/bin/bash

featureCounts -T 17 -p -a //media/user/LINUX/linux/kenza_work/STAR/star_ref/mm10.ncbiRefSeq.gtf -o featurecounts.txt *.sorted.bam

