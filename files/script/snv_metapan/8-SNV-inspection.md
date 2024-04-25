8-SNV-inspection
================
Hans Schrieke, Blandine Trouche and Julie Reveillaud
2024-04-25

-   <a href="#figure-s8-illustration-of-snv-filtering"
    id="toc-figure-s8-illustration-of-snv-filtering">Figure S8: illustration
    of SNV filtering</a>
-   <a href="#figures-s9-14-visualization-of-snvs-in-gene-context"
    id="toc-figures-s9-14-visualization-of-snvs-in-gene-context">Figures
    S9-14: visualization of SNVs in gene context</a>
-   <a href="#figures-s9-14-coverage-correlation-plots"
    id="toc-figures-s9-14-coverage-correlation-plots">Figures S9-14:
    coverage correlation plots</a>
    -   <a href="#mag-m11-gene-301" id="toc-mag-m11-gene-301">MAG M11: gene
        301</a>
    -   <a href="#mag-o11-gene-755" id="toc-mag-o11-gene-755">MAG O11: gene
        755</a>
    -   <a href="#mag-o03-gene-589" id="toc-mag-o03-gene-589">MAG O03: gene
        589</a>
    -   <a href="#mag-o07-gene-230" id="toc-mag-o07-gene-230">MAG O07: gene
        230</a>
    -   <a href="#mag-o12-gene-423" id="toc-mag-o12-gene-423">MAG O12: gene
        423</a>
    -   <a href="#mag-o12-gene-640" id="toc-mag-o12-gene-640">MAG O12: gene
        640</a>

``` r
library(tidyverse)
library(data.table)
library(ggplot2)
```

## Figure S8: illustration of SNV filtering

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
q + theme(legend.position = "bottom", legend.direction="vertical")
```

![](8-SNV-inspection_files/figure-gfm/all%20MAGs-1.png)<!-- -->

``` r
# save as svg for touch-ups in affinity designer
svg("11_SNV_figures/Fig1_barplot_stacked_by_MAG_nolegend.svg", width=8, height=6)
q + theme(legend.position = "none")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/Fig1_barplot_stacked_by_MAG.svg", width=12, height=6)
q 
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Figures S9-14: visualization of SNVs in gene context

List of splits / genes: - MAG M11, Culex_M11_000000026618_split_00001,
gene 301  
- MAG O11, Culex_O11_000000098617_split_00001, gene 755  
- MAG O03, Culex_O03_000000128791_split_00001, gene 589  
- MAG O07, Culex_O07_000000082262_split_00001, gene 230  
- MAG O12, Culex_O12_000000005804_split_00001, gene 423  
- MAG O12, Culex_O12_000000099075_split_00001, gene 640

``` zsh
# Get split coverage and SNVs from selected split of Wolbachia O11 MAG
# activate conda environment
conda activate anvio-7.1

# move to work directory
# cd your/path/here

# create specific folder for outputs
mkdir -p 11_SNV_figures

# gene level information extraction
for split in Culex_M11_000000026618_split_00001_gene_301 \
              Culex_O11_000000098617_split_00001_gene_755 \
              Culex_O03_000000128791_split_00001_gene_589 \
              Culex_O07_000000082262_split_00001_gene_230 \
              Culex_O12_000000005804_split_00001_gene_423 \
              Culex_O12_000000099075_split_00001_gene_640
do

  MAG=$(echo ${split} | cut -d '_' -f 2)
  gene=$(echo ${split} | cut -d '_' -f 7)
  
  mkdir -p 11_SNV_figures/${MAG}
  
  # get coverage from splits of interest for SNV visualization
  anvi-get-split-coverages -p 06_MERGED_references_mode/${MAG}_filtered/PROFILE.db \
                         -c 03_CONTIGS_references_mode/${MAG}-contigs.db \
                         --gene-caller-id ${gene} \
                         -o 11_SNV_figures/${MAG}/${split}_coverage.txt
  
  # generate variability table from genes of interests 
  anvi-gen-variability-profile -p 06_MERGED_references_mode/${MAG}_filtered/PROFILE.db \
                             -c 03_CONTIGS_references_mode/${MAG}-contigs.db \
                             --gene-caller-ids ${gene} \
                             --include-split-names \
                             --include-contig-names \
                             -o 11_SNV_figures/${MAG}/${split}_SNVs.txt
  
done 
```

``` r
# Filter SNVs to produce the different levels of plot
# packages
require(tidyverse)

splits <- c("Culex_M11_000000026618_split_00001_gene_301",
              "Culex_O11_000000098617_split_00001_gene_755",
              "Culex_O03_000000128791_split_00001_gene_589",
              "Culex_O07_000000082262_split_00001_gene_230",
              "Culex_O12_000000005804_split_00001_gene_423",
              "Culex_O12_000000099075_split_00001_gene_640")

for(i in splits){
  MAG <- sapply(strsplit(i, "_"), "[", 2)
  gene <- sapply(strsplit(i, "_"), "[", 7)
  
  snv <- read.table(paste0("11_SNV_figures/", MAG, "/", i, "_SNVs.txt"), header=TRUE)
  cov <- read.table(paste0("11_SNV_figures/", MAG, "/", i, "_coverage.txt"), header=TRUE)
  
  cov$gene_caller_id <- gsub(pattern = "_gene_.*", replacement = "", i)
  colnames(cov)[3] <- "split_name"

  # filter the SNV table by removing the outliers
  snv_no_outlier <- snv[snv$cov_outlier_in_contig==0 & snv$cov_outlier_in_split==0, ]

  # filter the SNV table by entropy and departure from consensus
  snv_entropy <- snv_no_outlier[snv_no_outlier$entropy >= 0.2, ]
  snv_dep <- snv_entropy[snv_entropy$departure_from_consensus >= 0.2, ]

  # export snv tables
  write.table(snv_no_outlier, paste0("11_SNV_figures/", MAG, "/",  i, "_SNVs_no_outliers.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_dep, paste0("11_SNV_figures/", MAG, "/", i, "_SNVs_filtered.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(cov, paste0("11_SNV_figures/", MAG, "/", i, "_coverage.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
}
```

``` zsh
# Generate figures with filtered SNV tables
# activate conda environment
conda activate anvio-7.1

for split in Culex_M11_000000026618_split_00001_gene_301 \
              Culex_O11_000000098617_split_00001_gene_755 \
              Culex_O03_000000128791_split_00001_gene_589 \
              Culex_O07_000000082262_split_00001_gene_230 \
              Culex_O12_000000005804_split_00001_gene_423 \
              Culex_O12_000000099075_split_00001_gene_640
do

  MAG=$(echo ${split} | cut -d '_' -f 2)
  gene=$(echo ${split} | cut -d '_' -f 7)
  
  echo sample_name"\n"${MAG}_filter20 > 11_SNV_figures/sample_of_interest.txt

  # generate pdf of splits of interest with SNVs
  anvi-script-visualize-split-coverages -i 11_SNV_figures/${MAG}/${split}_coverage.txt \
                                      -o 11_SNV_figures/${MAG}/${split}_inspect_SNV.pdf \
                                      --snv-data 11_SNV_figures/${MAG}/${split}_SNVs.txt \
                                      -s 11_SNV_figures/sample_of_interest.txt \
                                      -m 600

  # generate pdf of splits of interest with non-outliers SNVs
  anvi-script-visualize-split-coverages -i 11_SNV_figures/${MAG}/${split}_coverage.txt \
                                      -o 11_SNV_figures/${MAG}/${split}_inspect_SNV_no_outliers.pdf \
                                      --snv-data 11_SNV_figures/${MAG}/${split}_SNVs_no_outliers.txt \
                                      -s 11_SNV_figures/sample_of_interest.txt \
                                      -m 600

  # generate pdf of splits of interest with filtered SNVs
  anvi-script-visualize-split-coverages -i 11_SNV_figures/${MAG}/${split}_coverage.txt \
                                      -o 11_SNV_figures/${MAG}/${split}_inspect_SNV_filtered.pdf \
                                      --snv-data 11_SNV_figures/${MAG}/${split}_SNVs_filtered.txt \
                                      -s 11_SNV_figures/sample_of_interest.txt \
                                      -m 600

done
```

## Figures S9-14: coverage correlation plots

### MAG M11: gene 301

``` r
library(ggpmisc)

# M11
snv_M11 <- read.table("07_RAW_VAR_TABLES/SNV_M11.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

# we want to look at a specific gene
snv_M11_301 <- snv_M11[snv_M11$corresponding_gene_call == 301, ]
# we are only interested in intra-sample comparison
snv_M11_301 <- snv_M11_301[snv_M11_301$sample_id == "M11_filter20",]

# compute correlation info
r <- round(cor(snv_M11_301$departure_from_consensus, snv_M11_301$coverage), 2)
p <- cor.test(snv_M11_301$departure_from_consensus, snv_M11_301$coverage)$p.value

plot_noinfo <- ggplot(snv_M11_301, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_full <- ggplot(snv_M11_301, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 140, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 140, label=paste0("p-value = ", p), hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot_full
```

![](8-SNV-inspection_files/figure-gfm/coverage%20correlation%20M11-1.png)<!-- -->

``` r
# save plots as svg for touch-ups in Affinity designer
svg("11_SNV_figures/M11_gene301_noinfo.svg", width=8, height=6)
plot_noinfo
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/M11_gene301.svg", width=8, height=6)
plot_full
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### MAG O11: gene 755

``` r
# O11
snv_O11 <- read.table("07_RAW_VAR_TABLES/SNV_O11.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O11_755 <- snv_O11[snv_O11$corresponding_gene_call == 755, ]
# we are only interested in intra-sample comparison
snv_O11_755 <- snv_O11_755[snv_O11_755$sample_id == "O11_filter20",]

r <- round(cor(snv_O11_755$departure_from_consensus, snv_O11_755$coverage), 2)
p <- cor.test(snv_O11_755$departure_from_consensus, snv_O11_755$coverage)$p.value

plot_noinfo <- ggplot(snv_O11_755, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_full <- ggplot(snv_O11_755, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 420, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 420, label=paste0("p-value = ", p), hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot_full
```

![](8-SNV-inspection_files/figure-gfm/coverage%20correlation%20O11-1.png)<!-- -->

``` r
# save plots as svg for touch-ups in Affinity designer
svg("11_SNV_figures/O11_gene755_noinfo.svg", width=8, height=6)
plot_noinfo
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/O11_gene755.svg", width=8, height=6)
plot_full
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### MAG O03: gene 589

``` r
# O03
snv_O03 <- read.table("07_RAW_VAR_TABLES/SNV_O03.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O03_589 <- snv_O03[snv_O03$corresponding_gene_call == 589, ]
# we are only interested in intra-sample comparison
snv_O03_589 <- snv_O03_589[snv_O03_589$sample_id == "O03_filter20",]

r <- round(cor(snv_O03_589$departure_from_consensus, snv_O03_589$coverage), 2)
p <- cor.test(snv_O03_589$departure_from_consensus, snv_O03_589$coverage)$p.value

plot_noinfo <- ggplot(snv_O03_589, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_full <- ggplot(snv_O03_589, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 200, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 200, label=paste0("p-value = ", p), hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot_full
```

![](8-SNV-inspection_files/figure-gfm/coverage%20correlation%20O03-1.png)<!-- -->

``` r
# save plots as svg for touch-ups in Affinity designer
svg("11_SNV_figures/O03_gene589_noinfo.svg", width=8, height=6)
plot_noinfo
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/O03_gene589.svg", width=8, height=6)
plot_full
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### MAG O07: gene 230

``` r
# O07
snv_O07 <- read.table("07_RAW_VAR_TABLES/SNV_O07.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O07_230 <- snv_O07[snv_O07$corresponding_gene_call == 230, ]
snv_O07_230 <- snv_O07_230[snv_O07_230$sample_id == "O07_filter20",]

r <- round(cor(snv_O07_230$departure_from_consensus, snv_O07_230$coverage), 2)
p <- cor.test(snv_O07_230$departure_from_consensus, snv_O07_230$coverage)$p.value

plot_noinfo <- ggplot(snv_O07_230, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_full <- ggplot(snv_O07_230, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 130, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 130, label=paste0("p-value = ", p), hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot_full
```

![](8-SNV-inspection_files/figure-gfm/coverage%20correlation%20O07-1.png)<!-- -->

``` r
# save plots as svg for touch-ups in Affinity designer
svg("11_SNV_figures/O07_gene230_noinfo.svg", width=8, height=6)
plot_noinfo
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/O07_gene230.svg", width=8, height=6)
plot_full
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### MAG O12: gene 423

``` r
# O12
snv_O12 <- read.table("07_RAW_VAR_TABLES/SNV_O12.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O12_423 <- snv_O12[snv_O12$corresponding_gene_call == 423, ]
snv_O12_423 <- snv_O12_423[snv_O12_423$sample_id == "O12_filter20",]

r <- round(cor(snv_O12_423$departure_from_consensus, snv_O12_423$coverage), 2)
p <- cor.test(snv_O12_423$departure_from_consensus, snv_O12_423$coverage)$p.value

plot_noinfo <- ggplot(snv_O12_423, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_full <- ggplot(snv_O12_423, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 130, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 130, label=paste0("p-value = ", p), hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot_full
```

![](8-SNV-inspection_files/figure-gfm/coverage%20correlation%20O12%20first-1.png)<!-- -->

``` r
# save plots as svg for touch-ups in Affinity designer
svg("11_SNV_figures/O12_gene423_noinfo.svg", width=8, height=6)
plot_noinfo
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/O12_gene423.svg", width=8, height=6)
plot_full
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### MAG O12: gene 640

``` r
# O12
snv_O12 <- read.table("07_RAW_VAR_TABLES/SNV_O12.txt", 
                      header=TRUE, fill=TRUE, sep="\t", quote = "")

snv_O12_640 <- snv_O12[snv_O12$corresponding_gene_call == 640, ]
snv_O12_640 <- snv_O12_640[snv_O12_640$sample_id == "O12_filter20",]

r <- round(cor(snv_O12_640$departure_from_consensus, snv_O12_640$coverage), 2)
p <- cor.test(snv_O12_640$departure_from_consensus, snv_O12_640$coverage)$p.value


plot_noinfo <- ggplot(snv_O12_640, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_full <- ggplot(snv_O12_640, aes(y=departure_from_consensus, x=coverage)) + 
  geom_smooth(method="lm", col="black") + 
  geom_point(aes(color = as.character(cov_outlier_in_contig)), size = 3) +
  labs(x = "Coverage", y = "Departure from consensus", color = "Coverage outlier") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "left", label.y = "top", size = 4) +
  annotate("text", y = 0.4, x = 130, label=paste0("correlation = ", r), hjust = 0) +
  annotate("text", y = 0.35, x = 130, label=paste0("p-value = ", p), hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot_full
```

![](8-SNV-inspection_files/figure-gfm/coverage%20correlation%20O12%20second-1.png)<!-- -->

``` r
# save plots as svg for touch-ups in Affinity designer
svg("11_SNV_figures/O12_gene640_noinfo.svg", width=8, height=6)
plot_noinfo
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("11_SNV_figures/O12_gene640.svg", width=8, height=6)
plot_full
dev.off()
```

    ## quartz_off_screen 
    ##                 2
