6-generate-variability-tables
================
Hans Schrieke, Blandine Trouche and Julie Reveillaud
2024-04-25

-   <a href="#generate-raw-variability-tables"
    id="toc-generate-raw-variability-tables">Generate raw variability
    tables</a>
-   <a href="#add-scg-info-to-the-variability-tables"
    id="toc-add-scg-info-to-the-variability-tables">Add SCG info to the
    variability tables</a>
    -   <a href="#snv-tables" id="toc-snv-tables">SNV tables</a>
    -   <a href="#scv-tables" id="toc-scv-tables">SCV tables</a>
    -   <a href="#saav-tables" id="toc-saav-tables">SAAV tables</a>

## Generate raw variability tables

``` bash
# Create output directories
mkdir -p 07_RAW_VAR_TABLES
mkdir -p 08_FILTERED_VAR_TABLES

# Generate SNV tables 
conda activate anvio-7.1

for s in O03 O07 O11 M11 O12
do

    ## SNVs
    anvi-gen-variability-profile -p 06_MERGED_references_mode/${s}_filtered/PROFILE.db  \
                                 -c 03_CONTIGS_references_mode/${s}-contigs.db \
                                 -C Wolbachia \
                                 -b ${s} \
                                       --include-split-names \
                                 --include-contig-names \
                                 --compute-gene-coverage-stats \
                                 -o 07_RAW_VAR_TABLES/SNV_${s}.txt

    ## SAAVs
    anvi-gen-variability-profile -p 06_MERGED_references_mode/${s}_filtered/PROFILE.db \
                                 -c 03_CONTIGS_references_mode/${s}-contigs.db \
                                 -C Wolbachia \
                                 -b ${s} \
                                 --engine AA \
                                 --include-split-names \
                                 -o 07_RAW_VAR_TABLES/SAAV_${s}.txt
      
    ## SCVs  
    anvi-gen-variability-profile -p 06_MERGED_references_mode/${s}_filtered/PROFILE.db \
                                 -c 03_CONTIGS_references_mode/${s}-contigs.db \
                                 -C Wolbachia \
                                 -b ${s} \
                                 --engine CDN \
                                 --include-split-names \
                                 --include-site-pnps \
                                 -o 07_RAW_VAR_TABLES/SCV_${s}.txt

done

conda deactivate
```

## Add SCG info to the variability tables

### SNV tables

``` r
require(tidyverse)
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
samples <- c("O03", "O07", "O11", "O12", "M11")
#setwd("/Users/btrouche/Dropbox/RosaLind/SNV_Meren")

gene_clusters <- read.table("pangenomics/Wolbachia-all-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt", 
                            header=TRUE, fill=TRUE, sep="\t", quote = "")
mlst <- read.table("pangenomics/coverages-for-gene-clusters.txt", 
                   header=TRUE, fill=TRUE, sep="\t", quote = "")
mlst <- mlst %>% select(c(1:9))

for(i in samples){
  snv <- read.table(paste0("07_RAW_VAR_TABLES/SNV_", i,".txt"), header=TRUE)
  
  snv_SCG <- snv
  snv_SCG$entry_id <- paste0("p_", sprintf("%08d", as.numeric(snv_SCG$unique_pos)))
  
  gene_clusters_sample <- gene_clusters[gene_clusters$genome_name==i,]
  colnames(gene_clusters_sample)[5] <- "corresponding_gene_call"
  snv_SCG <- snv_SCG %>% merge(gene_clusters_sample, by="corresponding_gene_call", all.x = TRUE)
  snv_SCG$corresponding_gene_call <- snv_SCG$corresponding_gene_call %>% as.character()
  snv_SCG <- snv_SCG %>% select(c(entry_id, everything()))
  
  # add MLST + wsp + cid info
  snv_SCG <- snv_SCG %>% merge(mlst, by="gene_cluster_id", all.x = TRUE)
  
  # save tables
   write.table(snv_SCG, paste0("07_RAW_VAR_TABLES/SNV_", i, "_SCG.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
}
```

### SCV tables

``` r
require(tidyverse)
samples <- c("O03", "O07", "O11", "O12", "M11")

# add SCG info
for(i in samples){
  scv <- read.table(paste0("07_RAW_VAR_TABLES/SCV_", i,".txt"),
                    header=TRUE, fill=TRUE, sep="\t", quote = "")
  
  # prepare layers to add in profile.db of interactive
  scv_layers <- scv
  scv_layers$entry_id <- paste0("p_", sprintf("%08d", as.numeric(scv_layers$entry_id)))
  
  # add SCG info
  gene_clusters_sample <- gene_clusters[gene_clusters$genome_name==i,]
  colnames(gene_clusters_sample)[5] <- "corresponding_gene_call"
  scv_layers <- scv_layers %>% merge(gene_clusters_sample, by="corresponding_gene_call", all.x = TRUE)
  scv_layers$corresponding_gene_call <- scv_layers$corresponding_gene_call %>% as.character()
  scv_layers <- scv_layers %>% select(c(entry_id, everything()))
  
  # add MLST + wsp + cid info
  scv_layers <- scv_layers %>% merge(mlst, by="gene_cluster_id", all.x = TRUE)
  
  # save tables
  write.table(scv_layers, paste0("07_RAW_VAR_TABLES/SCV_", i, "_SCG.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
}
```

### SAAV tables

``` r
require(tidyverse)
samples <- c("O03", "O07", "O11", "O12", "M11")

# inter-sample level
for(i in samples){
  saav <- read.table(paste0("07_RAW_VAR_TABLES/SAAV_", i,".txt"),
                     header=TRUE, fill=TRUE, sep="\t", quote = "")
  
  # prepare layers to add in profile.db of interactive
  saav_layers <- saav
  saav_layers$entry_id <- paste0("p_", sprintf("%08d", as.numeric(saav_layers$entry_id)))
  
  # add SCG info
  gene_clusters_sample <- gene_clusters[gene_clusters$genome_name==i,]
  colnames(gene_clusters_sample)[5] <- "corresponding_gene_call"
  saav_layers <- saav_layers %>% merge(gene_clusters_sample, by="corresponding_gene_call", all.x = TRUE)
  saav_layers$corresponding_gene_call <- saav_layers$corresponding_gene_call %>% as.character()
  saav_layers <- saav_layers %>% select(c(entry_id, everything()))
  
  # add MLST + wsp + cid info
  saav_layers <- saav_layers %>% merge(mlst, by="gene_cluster_id", all.x = TRUE)
  
  # save tables
  write.table(saav_layers, paste0("07_RAW_VAR_TABLES/SAAV_", i, "_SCG.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
}
```
