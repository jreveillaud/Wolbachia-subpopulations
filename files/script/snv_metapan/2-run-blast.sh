#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Create output directory
cd ../../output && mkdir genes genes/cid genes/mlst


# Download cid genes from Bonneau et al. 2018 (10.1038/s41467-017-02749-w) and Beckmann et al. 2013 (à vérifier)

## From Bonneau et al. 2018: https://www.nature.com/articles/s41467-017-02749-w

cd genes/cid

for i in {963..996}
do
  wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MF444${i}&rettype=fasta" >> cidA_cidB_Bonneau.fa
done


## From Beckmann et al .2013: https://pubmed.ncbi.nlm.nih.gov/23856508/

sed '/^[[:space:]]*$/d' cidA_cidB_Bonneau.fa > cidA_cidB_Bonneau.fasta
wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KF114896.1&rettype=fasta" > cidA_cidB_Beckmann.fa # cidA
wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KX077602.1&rettype=fasta" >> cidA_cidB_Beckmann.fa # cidB
sed '/^[[:space:]]*$/d' cidA_cidB_Beckmann.fa >> cidA_cidB_Beckmann.fasta


# Download MLST + wsp genes from https://pubmlst.org

cd ../mlst

for i in gatB coxA hcpA ftsZ fbpA wsp
do
curl "https://pubmlst.org/bigsdb?db=pubmlst_wolbachia_seqdef&page=downloadAlleles&locus=${i}" > ${i}.fasta
done 

cat *.fasta > mlst_wsp.fasta



# Get best hits of cid / MLST genes from the Wolbachia wPipPel reference genome

## Export gene calls from wPipPel
cd .. && conda activate anvio-7.1

for i in cid mlst
do
  anvi-get-sequences-for-gene-calls -c ../Metapan/PAN/1_REFERENCE_GENOMES/wPip-PEL.db \
                                      -o ${i}/Wpip-PEL_genes.fa
done

conda deactivate


## Make blast database from genes of interest

### cid genes

cd cid && conda activate blast-2.12   

for i in Bonneau Beckmann
do
  makeblastdb -in cidA_cidB_${i}.fasta -dbtype nucl

  # blast search (note the extended outfmt that is necessary to compute the percent query alignment later)
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
  # NOTE: the following command works on linux systems but fails on Mac systems
  sort -k1,1 -k12,12gr -k11,11g -k3,3gr cidA_cidB_${i}_in_Wpip_hits.txt |\
    sort -u -k1,1 --merge >> cidA_cidB_${i}_in_Wpip_best_hits.txt

done


### MLST + WSP genes

cd ../mlst

conda activate blast-2.12

for i in mlst_wsp
do
  makeblastdb -in ${i}.fasta -dbtype nucl

  blastn -query Wpip-PEL_genes.fa \
         -db ${i}.fasta \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
         -out ${i}_in_Wpip_hits.txt

  awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' ${i}_in_Wpip_hits.txt > temp.txt
  mv temp.txt ${i}_in_Wpip_hits.txt
  
  echo -e "gene_callers_id\tmlst-wsp_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > ${i}_in_Wpip_best_hits.txt

  sort -k1,1 -k12,12gr -k11,11g -k3,3gr ${i}_in_Wpip_hits.txt |\
    sort -u -k1,1 --merge >> ${i}_in_Wpip_best_hits.txt

done

conda deactivate
    



# Get best hits of cid / MLST genes from each Wolbachia MAG

## Export the genes call from MAGs

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


## Make blast database from genes of interest

### cid genes

cd cid && conda activate blast-2.12

for s in M11 O11 O03 O07 O12
do 

  for i in Bonneau Beckmann
  do
  
  blastn -query Wolbachia_${s}_genes.fa \
       -db cidA_cidB_${i}.fasta \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -out cidA_cidB_${i}_in_${s}_hits.txt
    
  awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' cidA_cidB_${i}_in_${s}_hits.txt > temp.txt
  mv temp.txt cidA_cidB_${i}_in_${s}_hits.txt

  echo -e "gene_callers_id\tcidA-cidB_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > cidA_cidB_${i}_in_${s}_best_hits.txt

  sort -k1,1 -k12,12gr -k11,11g -k3,3gr cidA_cidB_${i}_in_${s}_hits.txt |\
    sort -u -k1,1 --merge >> cidA_cidB_${i}_in_${s}_best_hits.txt
  done

done


### MLST + WSP genes

cd ../mlst
for s in M11 O11 O03 O07 O12
do
  
  for i in mlst_wsp
  do

  blastn -query Wolbachia_${s}_genes.fa \
       -db ${i}.fasta \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -out ${i}_in_${s}_hits.txt

  awk '{print $0 "\t" (($8 - $7) * 100.0) / $13}' ${i}_in_${s}_hits.txt > temp.txt
  mv temp.txt ${i}_in_${s}_hits.txt

  echo -e "gene_callers_id\tmlst-wsp_assignment\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpct_alignment" \
    > ${i}_in_${s}_best_hits.txt

  sort -k1,1 -k12,12gr -k11,11g -k3,3gr ${i}_in_${s}_hits.txt |\
    sort -u -k1,1 --merge >> ${i}_in_${s}_best_hits.txt

  done

done

# Deactivate conda environment
conda deactivate