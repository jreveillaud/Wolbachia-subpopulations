#!/bin/sh
source /usr/local/genome/Anaconda3/etc/profile.d/conda.sh

# Shell to use for job execution
#$ -S /bin/bash

# Job name
#$ -N meren

# Nom de la queue
#$ -q long.q

# Export all environment variables
#$ -V

# Run command from current repertory
#$ -cwd

# CPUs number
#$ -pe thread 12


conda activate /home/hschrieke/work/anvio-7.1

for i in M11 O03 O07 O11 O12
do

anvi-gen-contigs-database -f Wolbachia-contigs-${i}.fa \
                          -o Wolbachia-${i}.db
                          
anvi-run-ncbi-cogs -c Wolbachia-${i}.db --num-threads 12 --cog-data-dir "/home/hschrieke/work/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/COG"

anvi-run-hmms -c Wolbachia-${i}.db

anvi-run-scg-taxonomy -c Wolbachia-${i}.db --num-threads 12 --scgs-taxonomy-data-dir "/home/hschrieke/work/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES/"


bowtie2-build Wolbachia-contigs-${i}.fa Wolbachia-${i}

bowtie2 -x Wolbachia-${i} \
        -1 Culex_${i}-QUALITY_PASSED_R1.fastq.gz \
        -2 Culex_${i}-QUALITY_PASSED_R2.fastq.gz \
        -S Wolbachia-${i}.sam \
        -p 12
        
samtools view -F 4 \
              -bS Wolbachia-${i}.sam \
              -o Wolbachia-${i}-RAW.bam


samtools sort Wolbachia-${i}-RAW.bam -o Wolbachia-${i}.bam
samtools index Wolbachia-${i}.bam

  
anvi-profile-blitz Wolbachia-${i}.bam \
                   -c Wolbachia-${i}.db \
                   -o OUTPUT-${i}.txt

anvi-profile-blitz Wolbachia-${i}.bam \
                   -c Wolbachia-${i}.db \
		               --gene-mode \
                   -o OUTPUT-genes-${i}.txt

done
                   
conda deactivate