#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1


# Annotate all contigs db with good reference databases
for s in ../../output/03_CONTIGS_references_mode/*-contigs.db
do
anvi-run-hmms -c $s -T 12 --just-do-it
anvi-run-ncbi-cogs -c $s -T 12 --cog-data-dir /Users/hschrieke/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/COG
anvi-run-scg-taxonomy -c $s -T 12 --scgs-taxonomy-data-dir /Users/hschrieke/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/SCG_TAXONOMY/GTDB
anvi-run-kegg-kofams -c $s -T 12 --kegg-data-dir /Users/hschrieke/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/KEGG
done

#Export function tables
for s in ../../output/03_CONTIGS_references_mode/*-contigs.db
do
anvi-export-functions -c $s -o ../../output/03_CONTIGS_references_mode/$s-functions-table.txt
done

# # Migrate profile dbs for anvio-7.1
# for s in M11 O03 O07 O11 O12
# do
# anvi-migrate ../06_MERGED_references_mode/$s/PROFILE.db --migrate-dbs-safely
# done

# Add collection and bin name to profile dbs
for s in M11 O03 O07 O11 O12
do
anvi-script-add-default-collection -p ../../output/06_MERGED_references_mode/$s/PROFILE.db \
-c ../../output/03_CONTIGS_references_mode/$s-contigs.db \
-C Wolbachia \
-b $s
done

# Summarize each sample
for s in M11 O03 O07 O11 O12
do
anvi-summarize -c ../../output/03_CONTIGS_references_mode/$s-contigs.db \
-p ../../output/06_MERGED_references_mode/$s/PROFILE.db \
-C Wolbachia \
--init-gene-coverages \
-o ../../output/03_CONTIGS_references_mode/$s-SUMMARY
done



# PREPARE WOLBACHIA REFERENCE GENOMES

# create a specific folder for pangenomic
mkdir ../../output/Metapan ../../output/Metapan/PAN && cd ../../output/Metapan/PAN
mkdir 1_REFERENCE_GENOMES && cd 1_REFERENCE_GENOMES

# Prepare Wolbachia reference genomes
cd ../../../../data/reference-genomes

for ref in wPip.fasta ABZA01.1.fsa_nt CTEH01.1.fsa_nt wAlbB.fna
do
# split name of reference
name=`echo $ref | awk 'BEGIN{FS=".f"}{print $1}'`

# simplify deflines in reference FASTA
anvi-script-reformat-fasta $ref \
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

# Rename db files
mv wPip.db wPip-PEL.db
mv ABZA01.1.db wPip-JHB.db
mv CTEH01.1.db wPip-MOL.db




# PANGENOME

# Create folder for Wolbachia pangenome output
# mkdir ../../output/Metapan ../../output/Metapan/PAN && cd ../../output/Metapan/PAN
# mkdir 1_REFERENCE_GENOMES && cd 1_REFERENCE_GENOMES
# 
# cd ../../../../data/reference-genomes

cd ../../output/Metapan/PAN
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


# Compute pangenom for metapan and all metagenomes
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