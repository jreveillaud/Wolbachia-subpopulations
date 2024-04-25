7b-SNP-figures
================
Hans Schrieke, Blandine Trouche and Julie Reveillaud
2024-04-25

-   <a href="#figure-1a-summary-of-snp-positions"
    id="toc-figure-1a-summary-of-snp-positions">Figure 1A: Summary of SNP
    positions</a>
-   <a href="#figure-1b-visualization-of-snps-in-gene-276-of-mag-o07"
    id="toc-figure-1b-visualization-of-snps-in-gene-276-of-mag-o07">Figure
    1B: Visualization of SNPs in gene 276 of MAG O07</a>
-   <a
    href="#figures-s5-s7-visualization-of-snps-in-gene-993-41-nd-611-of-mag-o12"
    id="toc-figures-s5-s7-visualization-of-snps-in-gene-993-41-nd-611-of-mag-o12">Figures
    S5-S7: Visualization of SNPs in gene 993, 41 nd 611 of MAG O12</a>

## Figure 1A: Summary of SNP positions

``` r
library(ggplot2)
library(patchwork)
library(forcats)

sample <- c("M11", "O11", "O03", "O07", "O12")

for(i in sample){
  
  # read SNP table
  snp <- read.table(paste0("08_FILTERED_VAR_TABLES/SNP/SNP_", i, ".txt"), 
                   header=TRUE, fill=TRUE, sep="\t", quote = "") 
  
  # extract SNP identification information
  gcs <- snp[, c("gene_cluster_id", "base_pos_in_codon", "codon_number", "sample_id")]
  gcs$MAG <- paste0("MAG ", i)

 assign(paste0("gene_clusters_", i), gcs)
}

# join tables
gcs <- rbind(gene_clusters_M11,
         gene_clusters_O03,
         gene_clusters_O07,
         gene_clusters_O11,
         gene_clusters_O12)

# create a new name sumarising gene cluser + codon number
gcs$GC_pos <- paste0(gcs$gene_cluster_id, "_", gcs$codon_number)


p <- ggplot(gcs, aes(x = factor(1), fill = sample_id)) +
  geom_bar(aes(fill = sample_id)) + coord_polar("y") +
  facet_grid(MAG~GC_pos, switch = "y", labeller = label_wrap_gen(3)) +
  scale_fill_manual(values = c("M11_filter20" = "gold", "O03_filter20" = "red", 
                               "O07_filter20" = "darkblue", "O11_filter20" = "darkgreen",
                               "O12_filter20" = "pink")) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(face = "bold", angle = 90),
        strip.text.y = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p
```

![](7b-SNP-figures_files/figure-gfm/plot%20of%20shared-1.png)<!-- -->

``` r
svg("12_SNP_figures/Fig1A_SNP_positions.svg", width=12, height=6)
p 
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Figure 1B: Visualization of SNPs in gene 276 of MAG O07

``` zsh
# Get split coverage and SNPs from selected splits of Wolbachia MAGs
# activate conda environment
conda activate anvio-7.1

# move to work directory
# cd your/path/here

# create specific folder for outputs
mkdir -p 12_SNP_figures

# gene level
for split in Culex_O07_000000124766_split_00001_gene_276
do

  MAG=$(echo ${split} | cut -d '_' -f 2)
  gene=$(echo ${split} | cut -d '_' -f 7)
  
  anvi-get-split-coverages -p 06_MERGED_references_mode/${MAG}_filtered/PROFILE.db \
                         -c 03_CONTIGS_references_mode/${MAG}-contigs.db \
                         --gene-caller-id ${gene} \
                         -o 12_SNP_figures/${split}_coverage.txt
  
  # generate variability table from splits of interests 
  # this has to be in a file
  echo ${split} > 12_SNP_figures/split_of_interest.txt
  
  anvi-gen-variability-profile -p 06_MERGED_references_mode/${MAG}_filtered/PROFILE.db \
                             -c 03_CONTIGS_references_mode/${MAG}-contigs.db \
                             --gene-caller-ids ${gene} \
                             --include-split-names \
                             --include-contig-names \
                             -o 12_SNP_figures/${split}_SNPs.txt
  
done 
```

``` zsh
# activate conda environment
conda activate anvio-7.1

for split in Culex_O07_000000124766_split_00001_gene_276
do

# generate pdf of splits of interest with variability in all samples
# Here it is possible that SNVs appear since we have not filtered the variability table 
# on coverage outliers, departure, etc... But it is interesting to consider the SNPs in context
anvi-script-visualize-split-coverages -i 12_SNP_figures/${split}_coverage.txt \
                                      -o 12_SNP_figures/${split}_inspect_SNP_all_samples.pdf \
                                      --snv-data 12_SNP_figures/${split}_SNPs.txt \
                                      -m 600

done
```

## Figures S5-S7: Visualization of SNPs in gene 993, 41 nd 611 of MAG O12

List of splits / genes:  
- MAG O12, Culex_O12_000000008388, gene 993  
- MAG O12, Culex_O12_000000068699,gene 41  
- MAG O12, Culex_O12_000000221252, gene 611

``` zsh
# Get split coverage and SNPs from selected splits of Wolbachia MAGs
# activate conda environment
conda activate anvio-7.1

# move to work directory
# cd your/path/here

# create specific folder for outputs
mkdir -p 12_SNP_figures

# gene level
for split in Culex_O12_000000008388_split_00001_gene_993 Culex_O12_000000068699_split_00002_gene_41 Culex_O12_000000221252_split_00001_gene_611
do

  MAG=$(echo ${split} | cut -d '_' -f 2)
  gene=$(echo ${split} | cut -d '_' -f 7)
  
  anvi-get-split-coverages -p 06_MERGED_references_mode/${MAG}_filtered/PROFILE.db \
                         -c 03_CONTIGS_references_mode/${MAG}-contigs.db \
                         --gene-caller-id ${gene} \
                         -o 12_SNP_figures/${split}_coverage.txt
  
  # generate variability table from splits of interests of O12
  # this has to be in a file
  echo ${split} > 12_SNP_figures/split_of_interest.txt
  
  anvi-gen-variability-profile -p 06_MERGED_references_mode/${MAG}_filtered/PROFILE.db \
                             -c 03_CONTIGS_references_mode/${MAG}-contigs.db \
                             --gene-caller-ids ${gene} \
                             --include-split-names \
                             --include-contig-names \
                             -o 12_SNP_figures/${split}_SNPs.txt
  
done 
```

``` zsh
# Generate figures with filtered SNP tables
# activate conda environment
conda activate anvio-7.1

for split in Culex_O12_000000008388_split_00001_gene_993 Culex_O12_000000068699_split_00002_gene_41 Culex_O12_000000221252_split_00001_gene_611
do

# generate pdf of splits of interest with SNPs all samples
anvi-script-visualize-split-coverages -i 12_SNP_figures/${split}_coverage.txt \
                                      -o 12_SNP_figures/${split}_inspect_SNP_all_samples.pdf \
                                      --snv-data 12_SNP_figures/${split}_SNPs.txt \
                                      -m 600
                                    
done
```
