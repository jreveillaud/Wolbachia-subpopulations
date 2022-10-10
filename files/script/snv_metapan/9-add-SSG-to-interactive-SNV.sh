#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

cd /Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-MAGs

conda activate anvio-7.1

for i in O03 O07 O11 O12 M11
do
anvi-import-misc-data variability-${i}-layers-SSG-wt-unique.txt \
                      -p ${i}_wolbachia_filtered/profile.db \
                      --target-data-table items \
                      --just-do-it

anvi-import-misc-data variability-${i}-layers.txt \
                      -p ${i}_wolbachia_filtered/profile.db \
                      --target-data-table items \
                      --just-do-it
                  
done
