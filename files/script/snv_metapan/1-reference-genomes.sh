#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1

# create a specific folder for pangenomic
mkdir ../output ../output/reference_genomes && cd ../output/reference_genomes


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

# Deactivate anvio environment
conda deactivate 