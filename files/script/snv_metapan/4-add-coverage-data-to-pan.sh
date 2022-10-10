#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1

# Import additional data with coverage informations in pangenomic databases
for i in metapan
do
anvi-import-misc-data ../../output/Metapan/METAPAN/1_COVERAGE_TABLE/coverages-for-gene-clusters.txt \
                      -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-${i}-PAN/Wolbachia-PAN.db \
                      --target-data-table items \
                      --just-do-it
                      
anvi-import-state -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-${i}-PAN/Wolbachia-PAN.db \
                  -s ../../metadata/state-metapan2.json \
                  -n draft2

done

# Deactivate anvi'o environment
conda deactivate

## display metapangenome
#anvi-import-state -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-PAN/Wolbachia-PAN.db -s ../../metadata/state-metapan.json -n metapan

# anvi-display-pan -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-PAN/Wolbachia-PAN.db \
#                  -g ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-PAN-GENOMES.db \
#                  --server-only -P 8080
# #                  
# anvi-display-pan -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-meta-wpip-PAN/Wolbachia-PAN.db \
#                  -g ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-meta-wpip-PAN-GENOMES.db \
#                  --server-only -P 8080
                 
## export state
# anvi-export-state -s state \
#                   -p profile-db  \
#                   -o path/to/output

# for i in O11 M11 O03 O07 O12
# do
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft3-${i}.json \
#                   -s draft3
#                   
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft4-${i}.json \
#                   -s draft4
#                   
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft4-gene-${i}.json \
#                   -s draft4_gene
#                   
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft4-snv-shared-${i}.json \
#                   -s draft4_snv_shared
#                   
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft4-gene-snv-shared-${i}.json \
#                   -s draft4_gene_snv_shared
#                   
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft5-gene-${i}.json \
#                   -s draft5_gene
#                   
# anvi-export-state -p /Volumes/Elements/Hiseq/metapangenomics/output-25-09-22/SNV/snv-tables/intra-MAGs/${i}_wolbachia_filtered/profile.db \
#                   -o /Volumes/Elements/Hiseq/metapangenomics/metadata/state-files2/state-draft5-gene-snv-shared-${i}.json \
#                   -s draft5_gene_snv_shared
# done