#!/bin/bash
echo "Processing fastqc..." >> logs/fastqc_results.out

fastqc $(find raw_data/ -name "*.fq.gz") -o fastqc_results/raw

echo "Successfully fasqchecked :)" >> logs/fastqc_results.out
