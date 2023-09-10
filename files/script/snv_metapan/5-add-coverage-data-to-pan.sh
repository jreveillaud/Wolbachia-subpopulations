#!/bin/bash
source /Users/hschrieke/opt/miniconda3/etc/profile.d/conda.sh

# Activate anvi'o environment
conda activate anvio-7.1

# Import additional data with coverage informations in pangenomics databases

  anvi-import-misc-data ../output/metapangenomics/coverages-for-gene-clusters.txt \
                        -p ../output/pangenomics/Wolbachia-all-PAN/Wolbachia-PAN.db \
                        --target-data-table items \
                        --just-do-it

  anvi-import-misc-data ../output/metapangenomics/coverages-for-gene-clusters-O11-M11.txt \
                        -p ../output/pangenomics/Wolbachia-O11-M11/Wolbachia-PAN.db \
                        --target-data-table items \
                        --just-do-it

# Deactivate anvi'o environment
conda deactivate


# Commands to display metapangenome

# anvi-display-pan -p Wolbachia-all-PAN/Wolbachia-PAN.db \
#                   -g Wolbachia-all-PAN-GENOMES.db \
#                   --server-only -P 8080
# 
# anvi-display-pan -p Wolbachia-O11-M11/Wolbachia-PAN.db \
#                   -g Wolbachia-all-PAN-GENOMES.db \
#                   --server-only -P 8080