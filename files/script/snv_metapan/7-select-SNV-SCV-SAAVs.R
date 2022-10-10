require(tidyverse)
require(openxlsx)

path <- "/Volumes/Elements/Hiseq/draft/tables/V2/V3/V4"
path <- "/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables"
setwd(path)

sample <- c("M11", "O11", "O03", "O07", "O12")
# i <- "M11"
# j <- 20
# test <- df[df$corresponding_gene_call==j,]

for(i in sample){
  df <- read.table(paste0("intra-MAGs/variability-", i, "-filtered-2.txt"), header=TRUE) # filtered SNVs (intra-MAGs)
  df2 <- read.table(paste0("SCVs/SCV_variability_", i, "_intra-MAGs_filtered.txt"), header=TRUE) # filtered SCVs (intra-MAGs)
  
  df_test <- data.frame()
  test <- data.frame()
  
  for(j in unique(df$corresponding_gene_call)){
    codon_numbers <- df[df$corresponding_gene_call==j, "codon_number"]
    codon_orders <- df[df$corresponding_gene_call==j, "codon_order_in_gene"]
    
    for(k in 1:length(codon_orders)){
      dff <- df2[df2$codon_number == codon_numbers[k] & 
                   df2$codon_order_in_gene == codon_orders[k] &
                   df2$corresponding_gene_call==j,] 
      df_test <- rbind(df_test, dff)
    }
  }
  df_test <- unique(df_test)
  write.table(df_test, paste0("SCVs/SCV_variability_", i, "_intra-MAGs_filtered_from_filtered_SNVs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  assign(paste0("snv_", i), df)
  assign(paste0("scv_raw_", i), df2)
  assign(paste0("scv_", i), df_test)
}

# i <- "M11"
# j <- 20
# test <- df2[df2$corresponding_gene_call==j,]
for(i in sample){
  # df2 <- read.table(paste0("SCVs/SCV_variability_", i, "_intra-MAGs_filtered.txt"), header=TRUE) # filtered SCVs (intra-MAGs)
  df2 <- eval(parse(text = paste0("scv_", i)))
  df3 <- read.table(paste0("SAAVs/SAAV_variability_", i, "_intra-MAGs_filtered.txt"), header=TRUE) # filtered SAAVs (intra-MAGs)
  
  # df2 <- read.xlsx("Supplementary Table 10.xlsx", sheet=i)
  # df3 <- read.xlsx("Supplementary Table BONUS 2.xlsx", sheet=i)
  
  df_test <- data.frame()
  test <- data.frame()
  
  for(j in unique(df2$corresponding_gene_call)){
    codon_numbers <- df2[df2$corresponding_gene_call==j, "codon_number"]
    codon_orders <- df2[df2$corresponding_gene_call==j, "codon_order_in_gene"]
    
    for(k in 1:length(codon_orders)){
      dff <- df3[df3$codon_number == codon_numbers[k] & 
                   df3$codon_order_in_gene == codon_orders[k] &
                   df3$corresponding_gene_call==j,] 
      df_test <- rbind(df_test, dff)
    }
  }
  write.table(df_test, paste0("SAAVs/SAAV_variability_", i, "_intra-MAGs_filtered_from_filtered_SCVs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  assign(paste0("saav_raw_", i), df3)
  assign(paste0("saav_", i), df_test)
}