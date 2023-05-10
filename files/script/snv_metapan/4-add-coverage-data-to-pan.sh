#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1

# Import additional data with coverage informations in pangenomics databases
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


# Commands to display metapangenome

# anvi-import-state -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-PAN/Wolbachia-PAN.db -s ../../metadata/state-metapan.json -n metapan

# anvi-display-pan -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-PAN/Wolbachia-PAN.db \
#                  -g ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-PAN-GENOMES.db \
#                  --server-only -P 8080

# anvi-display-pan -p ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-meta-wpip-PAN/Wolbachia-PAN.db \
#                  -g ../../output/Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-meta-wpip-PAN-GENOMES.db \
#                  --server-only -P 8080