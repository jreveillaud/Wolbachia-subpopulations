#!/bin/sh
source /usr/local/genome/Anaconda2-5.1.0/etc/profile.d/conda.sh

# activate anvio environment
conda activate anvio-7

# run anvio metagenomics workflow command 
anvi-run-workflow -w metagenomics -c megahit_ovaries.json --additional-params --cluster 'qsub -V -cwd -N metagenomics_ovaries -q highmem.q -o {log} -e {log} -pe thread {threads}' --jobs 100 --latency-wait 100

# deactivate anvio environment
conda deactivate