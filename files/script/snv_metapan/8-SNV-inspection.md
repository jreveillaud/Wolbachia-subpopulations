8-SNV-inspection
================
Hans Schrieke, Blandine Trouche and Julie Reveillaud
2024-04-25

-   <a href="#figure-s8" id="toc-figure-s8">Figure S8</a>
-   <a href="#figures-s9-14-a-to-c" id="toc-figures-s9-14-a-to-c">Figures
    S9-14 (A to C)</a>
-   <a href="#figures-s9-14-d-coverage-correlation-plots"
    id="toc-figures-s9-14-d-coverage-correlation-plots">Figures S9-14 (D);
    coverage correlation plots</a>
    -   <a href="#mag-m11-gene-301" id="toc-mag-m11-gene-301">MAG M11: gene
        301</a>
    -   <a href="#2b--mag-o11-gene-755" id="toc-2b--mag-o11-gene-755">2b- MAG
        O11: gene 755</a>
    -   <a href="#2c--mag-o03-gene-589" id="toc-2c--mag-o03-gene-589">2c- MAG
        O03: gene 589</a>
    -   <a href="#2d--mag-o07-gene-230" id="toc-2d--mag-o07-gene-230">2d- MAG
        O07: gene 230</a>
    -   <a href="#2e--mag-o12-gene-423-gene-640"
        id="toc-2e--mag-o12-gene-423-gene-640">2e- MAG O12: gene 423, gene
        640</a>

``` r
library(tidyverse)
library(data.table)
library(ggplot2)
```

## Figure S8

``` r
samples_loop <- c("O11", "M11", "O03", "O07", "O12")

df_all_MAGs <- data.frame()

# filter table for each MAG and bind them in an overall table
for (s in samples_loop) { 
  # load the full SNV table with wSCG info
  snv_MAG <- read.table(paste0("07_RAW_VAR_TABLES/SNV_", s, "_SCG.txt"), 
                        header=TRUE, fill=TRUE, sep="\t", quote = "")

  # only intra comparison here
  snv_MAG <- snv_MAG[snv_MAG$sample_id == paste0(s, "_filter20"), ]
  
  ### alluvial plot
  df <- data.frame("id" = snv_MAG$unique_pos_identifier, 
                   "SCG" = snv_MAG$SCG, "outlier_split" = snv_MAG$cov_outlier_in_split,
                   "outlier_contig" = snv_MAG$cov_outlier_in_contig, 
                   "departure_from_consensus" = snv_MAG$departure_from_consensus, 
                   "entropy" = snv_MAG$entropy)
  
  df[is.na(df$SCG), "SCG"] <- "Out of gene"
  df$SCG <- gsub(1, "wSCG", df$SCG)
  df$SCG <- gsub(0, "not wSCG", df$SCG)
  
  df[df$outlier_split == 1 | df$outlier_contig == 1, "outlier"] <- "outlier"
  df[df$outlier_split != 1 & df$outlier_contig != 1, "outlier"] <- "not outlier"
  
  df[df$departure_from_consensus >= 0.2 & df$entropy >= 0.2, "departure_filter"] <- "over 0.2"
  df[is.na(df$departure_filter), "departure_filter"] <- "under 0.2"
  
  df_short <- df[, c("SCG", "outlier", "departure_filter")]
  df_short <- df_short %>% group_by(SCG, outlier, departure_filter) %>% 
    summarise(total_count=n(),.groups = 'drop') %>%
    as.data.frame()
  
  df_short[df_short$SCG == "wSCG" & df_short$outlier == "not outlier" & 
             df_short$departure_filter == "over 0.2", "Status"] <- "Retained"
  df_short[df_short$SCG == "not wSCG", "Status"] <- "Discarded (not in wSCG)"
  df_short[df_short$SCG == "Out of gene", "Status"] <- "Not considered (out of gene)"
  df_short[df_short$SCG == "wSCG" & (df_short$outlier == "outlier" | 
                                       df_short$departure_filter == "under 0.2"), "Status"] <- "Filtered out (coverage outliers or below departure threshold)"
  df_short$MAG <- s
  
  # full dataframe
  df_all_MAGs <- rbind(df_all_MAGs, df_short)
}

# define explicit categories for plotting clearly
df_all_MAGs2 <- df_all_MAGs
df_all_MAGs2[df_all_MAGs2$Status == "Discarded (not in wSCG)", "SNV_type"] <- "SNVs in multi-copy genes (discarded)"
df_all_MAGs2[df_all_MAGs2$Status == "Not considered (out of gene)", "SNV_type"] <- "SNVs outside of genes (not considered)"

df_all_MAGs2[df_all_MAGs2$SCG == "wSCG" & df_all_MAGs2$outlier == "outlier", "SNV_type"] <- 
  "Coverage outlier SNVs in wSCGs (filtered out)"
df_all_MAGs2[df_all_MAGs2$SCG == "wSCG" & df_all_MAGs2$outlier != "outlier" & 
               df_all_MAGs2$departure_filter == "under 0.2", "SNV_type"] <- 
  "SNVs in wSCGs below 0.2 departure from consensus and entropy thresholds (filtered out)"
df_all_MAGs2[df_all_MAGs2$SCG == "wSCG" & df_all_MAGs2$outlier != "outlier" & 
               df_all_MAGs2$departure_filter == "over 0.2", "SNV_type"] <- 
  "SNVs in wSCGs passing filters (retained)"

# summarize to obtain frequencines for the barplot  
df_all_MAGs2 <- df_all_MAGs2 %>%
  group_by(SCG, SNV_type, MAG) %>%
  summarize(Freq = sum(total_count))
```

    ## `summarise()` has grouped output by 'SCG', 'SNV_type'. You can override using
    ## the `.groups` argument.

``` r
# define factor order for clear plotting
df_all_MAGs2$SNV_type <- factor(df_all_MAGs2$SNV_type, 
                         levels = c("SNVs in wSCGs passing filters (retained)",
                                    "SNVs in wSCGs below 0.2 departure from consensus and entropy thresholds (filtered out)",
                                    "Coverage outlier SNVs in wSCGs (filtered out)",
                                    "SNVs outside of genes (not considered)",
                                    "SNVs in multi-copy genes (discarded)")) 

## stacked version of the plot faceted by MAG
q <- ggplot(df_all_MAGs2, aes(x = SCG, y = Freq, fill = SNV_type)) + 
  geom_col(position = "stack") +
  facet_wrap(~MAG, strip.position = "bottom", ncol = 5) +
  geom_text(aes(label=Freq)) +
  scale_fill_manual(name = "SNV type",
                    values = c("SNVs outside of genes (not considered)" = "gray88",
                               "SNVs in wSCGs passing filters (retained)" = "brown3",
                               "SNVs in multi-copy genes (discarded)" = "aliceblue",
                               "SNVs in wSCGs below 0.2 departure from consensus and entropy thresholds (filtered out)" = "cadetblue2",
                               "Coverage outlier SNVs in wSCGs (filtered out)" = "gold")) +
  ylab("Number of variable positions") + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 12)) 
q
```

![](8-SNV-inspection_files/figure-gfm/all%20MAGs-1.png)<!-- -->

``` r
# save as svg for touch-ups in affinity designer
svg("11_snv_figures/Fig1_barplot_stacked_by_MAG_nolegend.svg", width=8, height=6)
q + theme(legend.position = "none")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_snv_figures/Fig1_barplot_stacked_by_MAG.svg", width=12, height=6)
q 
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Figures S9-14 (A to C)

``` zsh
# Get split coverage and SNVs from selected split of Wolbachia O11 MAG
# activate conda environment
conda activate anvio-7.1

# move to work directory
cd /Users/btrouche/Dropbox/RosaLind/SNV_Meren

# create specific folder for outputs
mkdir -p 11_snv_figures

# Split from O11 
# get coverage from splits of interest of O11 for SNV visualization
anvi-get-split-coverages -p 06_MERGED_references_mode/O11_filtered/PROFILE.db \
                         -c 03_CONTIGS_references_mode/O11-contigs.db \
                         --split-name Culex_O11_000000098617_split_00001 \
                         -o 11_snv_figures/Culex_O11_000000098617_split_00001_coverage.txt

anvi-get-split-coverages -p 06_MERGED_references_mode/O11_filtered/PROFILE.db \
                         -c 03_CONTIGS_references_mode/O11-contigs.db \
                         --gene-caller-id 755 \
                         -o 11_snv_figures/Culex_O11_000000098617_split_00001_gene_coverage.txt


# generate SNV table from splits of interests of O11
# this has to be in a file
echo Culex_O11_000000098617_split_00001 > 11_snv_figures/split_of_interest.txt

anvi-gen-variability-profile -p 06_MERGED_references_mode/O11_filtered/PROFILE.db \
                             -c 03_CONTIGS_references_mode/O11-contigs.db \
                             --splits-of-interest 11_snv_figures/split_of_interest.txt \
                             --include-split-names \
                             --include-contig-names \
                             -o 11_snv_figures/Culex_O11_000000098617_split_00001_SNVs.txt

anvi-gen-variability-profile -p 06_MERGED_references_mode/O11_filtered/PROFILE.db \
                             -c 03_CONTIGS_references_mode/O11-contigs.db \
                             --gene-caller-ids 755 \
                             --include-split-names \
                             --include-contig-names \
                             -o 11_snv_figures/Culex_O11_000000098617_split_00001_gene_SNVs.txt
```

``` r
# Filter SNVs
# packages
require(tidyverse)

# import SNV table
splits <- c("Culex_O11_000000098617_split_00001")
  
for(i in splits){
  snv <- read.table(paste0("11_snv_figures/", i, "_SNVs.txt"), header=TRUE)
  cov <- read.table(paste0("11_snv_figures/", i, "_coverage.txt"), header=TRUE)

  # filter the SNV table by removing the outliers
  snv_wt_outlier <- snv[snv$cov_outlier_in_contig==0 & snv$cov_outlier_in_split==0, ]

  # filter the SNV table by entropy and departure from consensus
  snv_entropy <- snv_wt_outlier[snv_wt_outlier$entropy>=0.2,]
  snv_dep <- snv_entropy[snv_entropy$departure_from_consensus>=0.2,]

  # export snv tables
  write.table(snv_wt_outlier, paste0("11_snv_figures/", i, "_SNVs_wt_outliers.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_dep, paste0("11_snv_figures/", i, "_SNVs_filtered.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(cov, paste0("11_snv_figures/", i, "_coverage.txt"), sep="\t", quote=FALSE, row.names=FALSE)
}


splits <- c("Culex_O11_000000098617_split_00001")

for(i in splits){
  snv <- read.table(paste0("11_snv_figures/", i, "_gene_SNVs.txt"), header=TRUE)
  cov <- read.table(paste0("11_snv_figures/", i, "_gene_coverage.txt"), header=TRUE)
  
  cov$gene_caller_id <- i
  colnames(cov)[3] <- "split_name"

  # filter the SNV table by removing the outliers
  snv_wt_outlier <- snv[snv$cov_outlier_in_contig==0 & snv$cov_outlier_in_split==0, ]

  # filter the SNV table by entropy and departure from consensus
  snv_entropy <- snv_wt_outlier[snv_wt_outlier$entropy>=0.2,]
  snv_dep <- snv_entropy[snv_entropy$departure_from_consensus>=0.2,]

  # export snv tables
  write.table(snv_wt_outlier, paste0("11_snv_figures/", i, "_gene_SNVs_wt_outliers.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_dep, paste0("11_snv_figures/", i, "_gene_SNVs_filtered.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(cov, paste0("11_snv_figures/", i, "_gene_coverage.txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
```

``` zsh
# Generate figures with filtered SNV tables
# activate conda environment
conda activate anvio-7.1

# move to work directory
cd /Users/btrouche/Dropbox/RosaLind/SNV_Meren

echo O11_filter20 > 11_snv_figures/sample_of_interest.txt
# don't forget to add header "sample_name"

# for split in Culex_O11_000000098617_split_00001_gene
for split in Culex_O11_000000098617_split_00001 Culex_O11_000000098617_split_00001_gene 
do

# generate pdf of split of interest (without SNVs)                        
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect.pdf \
                                      -m 600

# generate pdf of splits of interest with SNVs
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect_SNV.pdf \
                                      --snv-data 11_snv_figures/${split}_SNVs.txt \
                                      -s 11_snv_figures/sample_of_interest.txt \
                                      -m 600

# generate pdf of splits of interest with SNVs all samples
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect_SNV_all_samples.pdf \
                                      --snv-data 11_snv_figures/${split}_SNVs.txt \
                                      -m 600

# generate pdf of splits of interest with non-outliers SNVs
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect_SNV_wt_outliers.pdf \
                                      --snv-data 11_snv_figures/${split}_SNVs_wt_outliers.txt \
                                      -s 11_snv_figures/sample_of_interest.txt \
                                      -m 600

# generate pdf of splits of interest with non-outliers SNVs all samples
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect_SNV_wt_outliers_all_samples.pdf \
                                      --snv-data 11_snv_figures/${split}_SNVs_wt_outliers.txt \
                                      -m 600
                                      
# generate pdf of splits of interest with filtered SNVs
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect_SNV_filtered.pdf \
                                      --snv-data 11_snv_figures/${split}_SNVs_filtered.txt \
                                      -s 11_snv_figures/sample_of_interest.txt \
                                      -m 600

# generate pdf of splits of interest with filtered SNVs all samples
anvi-script-visualize-split-coverages -i 11_snv_figures/${split}_coverage.txt \
                                      -o 11_snv_figures/${split}_inspect_SNV_filtered_all_samples.pdf \
                                      --snv-data 11_snv_figures/${split}_SNVs_filtered.txt \
                                      -m 600

done
```

## Figures S9-14 (D); coverage correlation plots

### MAG M11: gene 301

``` r
library(ggpmisc)

# M11
snv_M11 <- read.table("07_RAW_VAR_TABLES/SNV_M11.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_M11_301 <- snv_M11[snv_M11$corresponding_gene_call == 301, ]

# à faire pour chaque gène wSCG qui montre un/des SNVs
r <- round(cor(snv_M11_301$departure_from_consensus, snv_M11_301$coverage), 2)
p <- cor.test(snv_M11_301$departure_from_consensus, snv_M11_301$coverage)$p.value

ggplot(snv_M11_301, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig))) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 140, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 140, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  ggtitle("GENE 301 in MAG M11")
```

### 2b- MAG O11: gene 755

``` r
# O11
snv_O11 <- read.table("07_RAW_VAR_TABLES/SNV_O11.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O11_755 <- snv_O11[snv_O11$corresponding_gene_call == 755, ]

# à faire pour chaque gène wSCG qui montre un/des SNVs
library(ggpmisc)

r <- round(cor(snv_O11_755$departure_from_consensus, snv_O11_755$coverage), 2)
p <- cor.test(snv_O11_755$departure_from_consensus, snv_O11_755$coverage)$p.value

ggplot(snv_O11_755, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = cov_outlier_in_contig)) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 420, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 420, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  geom_abline(slope = 0, intercept = 0.2, color = "red") +
  ggtitle("GENE 755 in MAG O11")

ggplot(snv_O11_755, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = cov_outlier_in_contig)) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 420, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 420, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  ggtitle("GENE 755 in MAG O11")
```

### 2c- MAG O03: gene 589

``` r
# O03
snv_O03 <- read.table("07_RAW_VAR_TABLES/SNV_O03.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O03_589 <- snv_O03[snv_O03$corresponding_gene_call == 589, ]

# à faire pour chaque gène wSCG qui montre un/des SNVs
r <- round(cor(snv_O03_589$departure_from_consensus, snv_O03_589$coverage), 2)
p <- cor.test(snv_O03_589$departure_from_consensus, snv_O03_589$coverage)$p.value

ggplot(snv_O03_589, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig))) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 200, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 200, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  ggtitle("GENE 589 in MAG O03")
```

### 2d- MAG O07: gene 230

``` r
# O07
snv_O07 <- read.table("07_RAW_VAR_TABLES/SNV_O07.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O07_230 <- snv_O07[snv_O07$corresponding_gene_call == 230, ]

# à faire pour chaque gène wSCG qui montre un/des SNVs
r <- round(cor(snv_O07_230$departure_from_consensus, snv_O07_230$coverage), 2)
p <- cor.test(snv_O07_230$departure_from_consensus, snv_O07_230$coverage)$p.value

ggplot(snv_O07_230, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig))) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 130, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 130, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  ggtitle("GENE 230 in MAG O07")
```

### 2e- MAG O12: gene 423, gene 640

``` r
# O12
snv_O12 <- read.table("07_RAW_VAR_TABLES/SNV_O12.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O12_423 <- snv_O12[snv_O12$corresponding_gene_call == 423, ]

# à faire pour chaque gène wSCG qui montre un/des SNVs
r <- round(cor(snv_O12_423$departure_from_consensus, snv_O12_423$coverage), 2)
p <- cor.test(snv_O12_423$departure_from_consensus, snv_O12_423$coverage)$p.value

ggplot(snv_O12_423, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig))) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 130, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 130, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  ggtitle("GENE 423 in MAG O12")
```

``` r
# O12
snv_O12 <- read.table("07_RAW_VAR_TABLES/SNV_O12.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O12_640 <- snv_O12[snv_O12$corresponding_gene_call == 640, ]

# à faire pour chaque gène wSCG qui montre un/des SNVs
r <- round(cor(snv_O12_640$departure_from_consensus, snv_O12_640$coverage), 2)
p <- cor.test(snv_O12_640$departure_from_consensus, snv_O12_640$coverage)$p.value

ggplot(snv_O12_640, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig))) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 130, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 130, label=paste0("p-value = ", p), hjust = 0) +
  theme_classic() +
  ggtitle("GENE 640 in MAG O12")
```
