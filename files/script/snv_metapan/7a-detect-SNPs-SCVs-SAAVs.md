7a-detect-SNPs-SCVs-SAAVs
================
Hans Schrieke, Blandine Trouche and Julie Reveillaud
2024-04-25

-   <a href="#snp-detection" id="toc-snp-detection">SNP detection</a>
-   <a href="#filtering-scvs-and-saavs-corresponding-to-the-snps"
    id="toc-filtering-scvs-and-saavs-corresponding-to-the-snps">Filtering
    SCVs and SAAVs corresponding to the SNPs</a>

## SNP detection

``` bash
conda activate anvio-7.1

mkdir -p 10_SUMMARY

# summarizing the mapping information for all Wolbachia MAGs to extract detection (breadth of coverage)
for i in M11 O03 O07 O11 O12
do

anvi-summarize -p 06_MERGED_references_mode/${i}_filtered/PROFILE.db \
                -c 03_CONTIGS_references_mode/${i}-contigs.db \
                --init-gene-coverages \
                -C Wolbachia \
                -o 10_SUMMARY/${i}_filtered 
done

conda deactivate

# create the folder for the next step
mkdir -p 08_FILTERED_VAR_TABLES/SNP
```

``` r
require(tidyverse)
```

    ## Loading required package: tidyverse

    ## Warning: package 'ggplot2' was built under R version 4.2.3

    ## Warning: package 'dplyr' was built under R version 4.2.3

    ## Warning: package 'stringr' was built under R version 4.2.3

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
require(openxlsx)
```

    ## Loading required package: openxlsx

``` r
samples <- c("O11", "M11", "O03", "O07", "O12")

for (s in samples) { 
  
  # Load the SNV table (raw with SCG indications)
  snv <- read.table(paste0("07_RAW_VAR_TABLES/SNV_", s, "_SCG.txt"), 
                    header=TRUE, fill=TRUE, sep="\t", quote = "")
  
  # Filter down to what we are interested in: no coverage outliers, in wSCGs, with departure from reference >= 0.98
  snp_filtered <- snv[snv$cov_outlier_in_contig ==0 & snv$cov_outlier_in_split == 0, ]
  snp_filtered <- snp_filtered[snp_filtered$SCG == 1 & !is.na(snp_filtered$SCG), ]
  snp_filtered <- snp_filtered[snp_filtered$departure_from_reference >= 0.98, ]
  print(dim(snp_filtered))
  print(min(snp_filtered$departure_from_reference))
  
  # Get gene number + metagenome
  gene_df <- data.frame(snp_filtered$corresponding_gene_call, snp_filtered$sample_id)
  colnames(gene_df) <- c("gene_callers_id", "sample_id")
  
  # Load detection 
  gene_detection <- read.table(paste0("10_SUMMARY/", s, "_filtered/bin_by_bin/", s, "/", s, "-gene_detection.txt"), 
                               header=TRUE, fill=TRUE, sep="\t")
  gene_detection_long <- pivot_longer(gene_detection, cols = !"gene_callers_id", 
                                      names_to = "sample_id", values_to = "detection")
  
  # Check detection by merging the tables
  gene_df <-merge(gene_df, gene_detection_long, 
                  by=c("gene_callers_id", "sample_id"), all.x = TRUE)
  print(min(gene_df$detection))
  # minimum detection is 1, meaning that it all considered genes are very well covered by mapping
  
  # Write final SNP table
  write.table(snp_filtered, paste0("08_FILTERED_VAR_TABLES/SNP/SNP_", s, ".txt"), 
              sep="\t", quote=FALSE, row.names=FALSE) 
}
```

    ## [1] 37 65
    ## [1] 0.9875
    ## [1] 1
    ## [1] 37 65
    ## [1] 0.9875
    ## [1] 1
    ## [1] 51 65
    ## [1] 0.9805195
    ## [1] 1
    ## [1] 54 65
    ## [1] 0.9805195
    ## [1] 1
    ## [1] 45 65
    ## [1] 0.9896907
    ## [1] 1

## Filtering SCVs and SAAVs corresponding to the SNPs

``` r
# Load packages
require(tidyverse)
require(openxlsx)

# Select only SCVs generated from filtered SNPs 
sample <- c("M11", "O11", "O03", "O07", "O12")

for(i in sample){
  
  # read SNP table
  df <- read.table(paste0("08_FILTERED_VAR_TABLES/SNP/SNP_", i, ".txt"), 
                   header=TRUE, fill=TRUE, sep="\t", quote = "") 
  
  # read raw SCV table
  df2 <- read.table(paste0("07_RAW_VAR_TABLES/SCV_", i, "_SCG.txt"), 
                    header=TRUE, fill=TRUE, sep="\t", quote = "") 
  
  # create empty dataframe that will contain the results
  df3 <- data.frame()
  
  for(k in 1:nrow(df)){
    
    # select SCV with same sample id, gene call, codon number and codon order as in the SNP table
    # this takes care of the wSCG filter, as well as coverage outlier filter, but not departure
    dff <- df2[df2$sample_id == df$sample_id[k] &
                 df2$corresponding_gene_call == df$corresponding_gene_call[k] &
                 df2$codon_order_in_gene == df$codon_order_in_gene[k] &
                 df2$codon_number == df$codon_number[k], ] 
    
    # bind it to main SCV table generated from SNP
    df3 <- rbind(df3, dff)
  }
  # there could technically be less SCVs than SNPs, but not more
  # One more check can be done manually by looking at the SNP and the codon transition in the SCV

  # filter on departure from reference
  df3 <- df3[df3$departure_from_reference >= 0.98, ]
  
  # write table
  write.table(df3, paste0("08_FILTERED_VAR_TABLES/SCV/SCV_corresponding_to_SNP_", i, ".txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  # assign tables for each sample
  assign(paste0("snp_", i), df)
  assign(paste0("scv_raw_", i), df2)
  assign(paste0("scv_", i), df3)
}


# Select SAAVs generated from SCVs generated from filtered SNPs
for(i in sample){
  df <- eval(parse(text = paste0("snp_", i)))
  df3 <- read.table(paste0("07_RAW_VAR_TABLES/SAAV_", i, "_SCG.txt"),
                    header=TRUE, fill=TRUE, sep="\t", quote = "") 
  
  df4 <- data.frame()

  for(k in 1:nrow(df)){
    
    # select SAAV with same sample id, gene call, codon number and codon order as in the SNP table
    dff <- df3[df3$sample_id == df$sample_id[k] &
                 df3$corresponding_gene_call == df$corresponding_gene_call[k] &
                 df3$codon_order_in_gene == df$codon_order_in_gene[k] &
                 df3$codon_number == df$codon_number[k], ] 
    
    # bind it to main SCV table generated from SNP
    df4 <- rbind(df4, dff)
  }
  
  # filter on departure from reference
  df4 <- df4[df4$departure_from_reference >= 0.98, ]
  
  write.table(df4, paste0("08_FILTERED_VAR_TABLES/SAAV/SAAV_corresponding_to_SNP_", i, ".txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  assign(paste0("saav_raw_", i), df3)
  assign(paste0("saav_", i), df4)
}
```
