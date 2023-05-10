#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1


## PREPARE WOLBACHIA REFERENCE GENOMES

# create a specific folder for pangenomic
mkdir ../../output/Metapan ../../output/Metapan/PAN && cd ../../output/Metapan/PAN
mkdir 1_REFERENCE_GENOMES && cd 1_REFERENCE_GENOMES


# Download Wolbachia reference genomes

## wPipPel
wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_010981.1&rettype=fasta" > wPip-PEL.fasta

## wPipMol
curl "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/723/225/GCA_000723225.2_Wolbachia_endosymbiont_wPip_Mol_of_Culex_molestus/GCA_000723225.2_Wolbachia_endosymbiont_wPip_Mol_of_Culex_molestus_genomic.fna.gz" | gunzip > wPip-MOL.fasta

## wPipJHB
curl "https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/AB/ZA/ABZA01/ABZA01.1.fsa_nt.gz" | gunzip > wPip-JHB.fasta


# Prepare Wolbachia reference genomes for anvi'o

for ref in wPip-PEL wPip-MOL wPip-JHB
do

  # split name of reference
  name=`echo $ref.fasta | awk 'BEGIN{FS=".f"}{print $1}'`

  # simplify deflines in reference FASTA
  anvi-script-reformat-fasta $ref.fasta \
                           --simplify-names \
                           -o $name.fa

  # generate a contigs database
  anvi-gen-contigs-database -f $name.fa \
                          -o $name.db

  # run hmms, kegg, cogs and taxonomy
  anvi-run-hmms -c $name.db
  anvi-run-kegg-kofams -c $name.db --kegg-data-dir "/Users/hschrieke/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/KEGG" -T 12
  anvi-run-ncbi-cogs -c $name.db --cog-data-dir "/Users/hschrieke/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/COG" -T 12
  anvi-run-scg-taxonomy -c $name.db --scgs-taxonomy-data-dir "/Users/hschrieke/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/SCG_TAXONOMY/GTDB" -T 12
done



## PANGENOMES

# Create folder for Wolbachia pangenome output

cd ../
mkdir 2_WOLBACHIA_PAN && cd 2_WOLBACHIA_PAN

## Create a genome storage with internal and external genomes

### All metagenomes + references (including wAlb)
anvi-gen-genomes-storage -i ../../../../metadata/metapan/all-internal.txt \
                         -e ../../../../metadata/metapan/all-reference-genomes.txt \
                         -o Wolbachia-all-PAN-GENOMES.db


# Prepare genome name files

## Metapan
for s in M11 O11 wPipPel
do
echo "$s" >> genomes-list-metapan.txt
done

## All metagenomes
for s in M11 O11 O03 O07 O12
do
echo "$s" >> genomes-list-all-metagenomes.txt
done

# Each metagenomes + 3 selected Wolbachia reference genomes
for s in wPipPel wPipJHB wPipMol
do
echo "$s" >> genomes-list-references.txt
done

for s in M11 O11 O03 O07 O12
do
cp genomes-list-references.txt genomes-list-$s-references.txt
echo "$s" >> genomes-list-$s-references.txt
done

rm genomes-list-references.txt


# Compute pangenome for metapan and all metagenomes
for i in metapan all-metagenomes
do
anvi-pan-genome -g Wolbachia-all-PAN-GENOMES.db \
                -o Wolbachia-${i}-PAN \
                --use-ncbi-blast \
                --mcl-inflation 10 \
                -T 6 \
                --genome-names genomes-list-${i}.txt \
                -n Wolbachia

anvi-script-add-default-collection -p Wolbachia-${i}-PAN/Wolbachia-PAN.db

anvi-summarize -p Wolbachia-${i}-PAN/Wolbachia-PAN.db \
               -g Wolbachia-all-PAN-GENOMES.db \
               -C DEFAULT \
               -o Wolbachia-${i}-PAN-SUMMARY

gzip -d Wolbachia-${i}-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt.gz
done


### Individual metagenomes + references (without wAlb)
for i in O11 M11 O03 O07 O12
do
anvi-pan-genome -g Wolbachia-all-PAN-GENOMES.db \
                -o Wolbachia-${i}-ref-PAN \
                --use-ncbi-blast \
                --mcl-inflation 10 \
                -T 6 \
                --genome-names genomes-list-${i}-references.txt \
                -n Wolbachia

anvi-script-add-default-collection -p Wolbachia-${i}-ref-PAN/Wolbachia-PAN.db

anvi-summarize -p Wolbachia-${i}-ref-PAN/Wolbachia-PAN.db \
               -g Wolbachia-all-PAN-GENOMES.db \
               -C DEFAULT \
               -o Wolbachia-${i}-ref-PAN-SUMMARY

gzip -d Wolbachia-${i}-ref-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt.gz

done


# Deactivate anvio environment
conda deactivate 