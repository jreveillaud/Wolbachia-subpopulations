#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

############## wPipPel

# Cid genes
# cd /Volumes/LaCie/Hiseq/metapangenomics/output
cd ../../output
mkdir genes genes/cid genes/mlst
cd genes/cid

## Download

### Bonneau (cidA and cidB)
for i in {963..996}
  do
  wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MF444${i}&rettype=fasta" >> cidA_cidB_Bonneau.fa
  done

### Beckmann 
sed '/^[[:space:]]*$/d' cidA_cidB_Bonneau.fa > cidA_cidB_Bonneau.fasta

#### cidA
wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KF114896.1&rettype=fasta" > cidA_cidB_Beckmann.fa

#### cidB
wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KX077602.1&rettype=fasta" >> cidA_cidB_Beckmann.fa
sed '/^[[:space:]]*$/d' cidA_cidB_Beckmann.fa >> cidA_cidB_Beckmann.fasta


# Export the genes call from wpipPel
cd ..
conda activate anvio-7.1
for i in cid mlst
do
anvi-get-sequences-for-gene-calls -c ../../data/reference-genomes/wPip-PEL.db \
                                      -o ${i}/Wpip-PEL_genes.fa
done
conda deactivate


# Make blast database from genes of interest
cd cid
conda activate blast-2.12   
for i in Bonneau Beckmann
do
makeblastdb -in cidA_cidB_${i}.fasta -dbtype nucl

# blast serach (note the extended outfmt that is necessary to compute
# the percent query alignment later)
blastn -query Wpip-PEL_genes.fa \
       -db cidA_cidB_${i}.fasta \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -out cidA_cidB_${i}_in_Wpip_hits.txt


# add percent alignment of the query as a new column 
awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' cidA_cidB_${i}_in_Wpip_hits.txt > temp.txt
mv temp.txt cidA_cidB_${i}_in_Wpip_hits.txt

# generate an empty output file with the correct header:
echo -e "gene_callers_id\tcidA-cidB_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > cidA_cidB_${i}_in_Wpip_best_hits.txt


# sort the BLAST output to get the best hit for each gene call
# NOTE: the following command works on linux systems but fails on
# Mac systems
sort -k1,1 -k12,12gr -k11,11g -k3,3gr cidA_cidB_${i}_in_Wpip_hits.txt |\
    sort -u -k1,1 --merge >> cidA_cidB_${i}_in_Wpip_best_hits.txt

done


# MSLT + WSP 
cd ../mlst

cp ../../../metadata/genes/*.fas . 

cat *.fas > mlst_wsp.fasta

conda activate blast-2.12
for i in mlst_wsp
do
makeblastdb -in ${i}.fasta -dbtype nucl

# blast serach (note the extended outfmt that is necessary to compute
# the percent query alignment later)
blastn -query Wpip-PEL_genes.fa \
       -db ${i}.fasta \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -out ${i}_in_Wpip_hits.txt


# add percent alignment of the query as a new column
awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' ${i}_in_Wpip_hits.txt > temp.txt
mv temp.txt ${i}_in_Wpip_hits.txt

# generate an empty output file with the correct header:
echo -e "gene_callers_id\tmlst-wsp_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > ${i}_in_Wpip_best_hits.txt


# sort the BLAST output to get the best hit for each gene call
# NOTE: the following command works on linux systems but fails on
# Mac systems
sort -k1,1 -k12,12gr -k11,11g -k3,3gr ${i}_in_Wpip_hits.txt |\
    sort -u -k1,1 --merge >> ${i}_in_Wpip_best_hits.txt

done

conda deactivate
    
   
   
   
   
############## For each sample 

# Export the genes call from wpipPel
cd ../../genes



for s in M11 O11 O03 O07 O12
do

  conda activate anvio-7.1
  for i in cid mlst
  do
  anvi-get-sequences-for-gene-calls -c ../03_CONTIGS_references_mode/${s}-contigs.db \
                                  -G ${s} \
                                  -o ${i}/Wolbachia_${s}_genes.fa
  done
  conda deactivate
done


cd cid
conda activate blast-2.12
for s in M11 O11 O03 O07 O12
do 

  for i in Bonneau Beckmann
  do
  
  # blast search (note the extended outfmt that is necessary to compute the percent query alignment later)
  blastn -query Wolbachia_${s}_genes.fa \
       -db cidA_cidB_${i}.fasta \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -out cidA_cidB_${i}_in_${s}_hits.txt
    
  # add percent alignment of the query as a new column 
  awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' cidA_cidB_${i}_in_${s}_hits.txt > temp.txt
  mv temp.txt cidA_cidB_${i}_in_${s}_hits.txt

  # generate an empty output file with the correct header:
  echo -e "gene_callers_id\tcidA-cidB_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > cidA_cidB_${i}_in_${s}_best_hits.txt


  # sort the BLAST output to get the best hit for each gene call
  # NOTE: the following command works on linux systems but fails on
  # Mac systems
  sort -k1,1 -k12,12gr -k11,11g -k3,3gr cidA_cidB_${i}_in_${s}_hits.txt |\
    sort -u -k1,1 --merge >> cidA_cidB_${i}_in_${s}_best_hits.txt
  done

done


cd ../mlst
for s in M11 O11 O03 O07 O12
do
  #conda activate blast-2.12
  for i in mlst_wsp
  do

  # blast serach (note the extended outfmt that is necessary to compute
  # the percent query alignment later)
  blastn -query Wolbachia_${s}_genes.fa \
       -db ${i}.fasta \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -out ${i}_in_${s}_hits.txt


  # add percent alignment of the query as a new column
  awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' ${i}_in_${s}_hits.txt > temp.txt
  mv temp.txt ${i}_in_${s}_hits.txt

  # generate an empty output file with the correct header:
  echo -e "gene_callers_id\tmlst-wsp_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > ${i}_in_${s}_best_hits.txt


  # sort the BLAST output to get the best hit for each gene call
  # NOTE: the following command works on linux systems but fails on
  # Mac systems
  sort -k1,1 -k12,12gr -k11,11g -k3,3gr ${i}_in_${s}_hits.txt |\
    sort -u -k1,1 --merge >> ${i}_in_${s}_best_hits.txt

  done
done

conda deactivate