#!/bin/sh

## 1- Filter out reads with MAPQ < 20 from bam files

# activate samtools environment
conda activate samtools
# samtools version 1.19.1

for m1 in M11 O03 O07 O11 O12; do
  echo "$m1"
  
  for m2 in M11 O03 O07 O11 O12; do
  
    # filter low quality mapping with the latest version of samtools 
    samtools view -bSq 20 ../output/04_MAPPING_references_mode/${m1}/${m2}.bam > 
        ../output/04_MAPPING_references_mode/${m1}/${m2}_filter20.bam

    # index the bam files        
    samtools index ../output/04_MAPPING_references_mode/${m1}/${m2}_filter20.bam
  
  done

done

conda deactivate

## 2- Create filtered profiles and merge

# activate anvi'o environment
conda activate anvio-7.1 

# profile 
for m1 in M11 O03 O07 O11 O12; do
  echo "$m1"
  
  for m2 in M11 O03 O07 O11 O12; do
  
    anvi-profile -i ../output/04_MAPPING_references_mode/${m1}/${m2}_filter20.bam  \
              -c ../output/03_CONTIGS_references_mode/${m1}-contigs.db \
              -o ../output/05_ANVIO_PROFILE_references_mode/${m1}/${m2}_filter20 \
              --sample-name ${m2}_filter20 --profile-SCVs \
              --overwrite-output-destinations \
              >> ../output/00_LOGS_references_mode/anvi-profile_${m1}_${m2}_filter20.log 2>&1
              
  done 

done

# merge
for m in M11 O03 O07 O11 O12; do
  echo "$m"
  
  anvi-merge -c ../output/03_CONTIGS_references_mode/${m}-contigs.db \
          ../output/05_ANVIO_PROFILE_references_mode/${m}/*_filter20/PROFILE.db \
          -o ../output/06_MERGED_references_mode/${m}_filtered \
          -S ${m}_filter20 \
          --overwrite-output-destinations \
          >> ../output/00_LOGS_references_mode/anvi-merge_${m}_filter20.log 2>&1
          
  anvi-script-add-default-collection -p ../output/06_MERGED_references_mode/${m}_filtered/PROFILE.db \
        -c ../output/03_CONTIGS_references_mode/${m}-contigs.db \
        -C Wolbachia -b ${m}

done 

conda deactivate


