#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1

# Create a specific folder for pangenomic
mkdir ../output/pangenomics && cd ../output/pangenomics


## Create a genome storage with internal and external genomes

### All metagenomes + references (including wAlb)
pwd

anvi-gen-genomes-storage -i ../../metadata/all-internal-metagenomes.txt \
                         -e ../../metadata/all-reference-genomes.txt \
                         -o Wolbachia-all-PAN-GENOMES.db


# Prepare genome name files

## All metagenomes and reference genomes
for s in M11 O11 O03 O07 O12 wPipPEL wPipJHB wPipMOL
do
echo "$s" >> all-genomes.txt
done

for s in M11 O11 wPipPEL
do
echo "$s" >> genomes-from-11.txt
done

# Compute pangenome for metapan and all metagenomes
anvi-pan-genome -g Wolbachia-all-PAN-GENOMES.db \
                -o Wolbachia-all-PAN \
                --use-ncbi-blast \
                --mcl-inflation 10 \
                -T 6 \
                --genome-names all-genomes.txt \
                -n Wolbachia

anvi-pan-genome -g Wolbachia-all-PAN-GENOMES.db \
                -o Wolbachia-O11-M11 \
                --use-ncbi-blast \
                --mcl-inflation 10 \
                -T 6 \
                --genome-names genomes-from-11.txt \
                -n Wolbachia

anvi-script-add-default-collection -p Wolbachia-all-PAN/Wolbachia-PAN.db
anvi-script-add-default-collection -p Wolbachia-O11-M11/Wolbachia-PAN.db

anvi-summarize -p Wolbachia-all-PAN/Wolbachia-PAN.db \
               -g Wolbachia-all-PAN-GENOMES.db \
               -C DEFAULT \
               -o Wolbachia-all-PAN-SUMMARY

anvi-summarize -p Wolbachia-O11-M11/Wolbachia-PAN.db \
               -g Wolbachia-all-PAN-GENOMES.db \
               -C DEFAULT \
               -o Wolbachia-O11-M11-PAN-SUMMARY

gzip -d Wolbachia-all-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt.gz
gzip -d Wolbachia-O11-M11-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt.gz

# Deactivate anvio environment
conda deactivate 