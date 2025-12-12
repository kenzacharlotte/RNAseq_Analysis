#!/bin/sh

#!/bin/bash

mkdir -p ~/fastqc_results
fastqc $(find . -name "*.fq.gz") -o ~/fastqc_results
