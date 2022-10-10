#!/bin/bash
source /usr/local/genome/Anaconda3/etc/profile.d/conda.sh

# Shell to use for job execution
#$ -S /bin/bash

# Job name
#$ -N snakemake-metapangenomics

# Nom de la queue
#$ -q long.q

# Export all environment variables
#$ -V

# Run command from current repertory
#$ -cwd

# CPUs number
#$ -pe thread 2

# Hostname
#-l hostname=n3

conda activate anvio-7.1
anvi-run-workflow -w metagenomics -c config-references-mode-full.json --additional-params --cluster 'qsub -V -cwd -N {rule} -q long.q -o {log} -e {log} -pe thread {threads} -S /bin/bash' --jobscript jobscript.sh --rerun-incomplete --jobs 100 --latency-wait 100 --unlock
anvi-run-workflow -w metagenomics -c config-references-mode-full.json --additional-params --cluster 'qsub -V -cwd -N {rule} -q long.q -o {log} -e {log} -pe thread {threads} -S /bin/bash' --jobscript jobscript.sh --rerun-incomplete --jobs 100 --latency-wait 100
conda deactivate
