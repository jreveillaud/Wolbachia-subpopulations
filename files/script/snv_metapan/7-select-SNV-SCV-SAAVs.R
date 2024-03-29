# Load packages
require(tidyverse)
require(openxlsx)

# Set path
path <- "../output/variability"
setwd(path)

# Select only SCVs generated from filtered SNVs 

sample <- c("M11", "O11", "O03", "O07", "O12")

# sample <- c("M11")
# i <- "M11"
# j <- 20
# k <- 1
# test <- df[df$corresponding_gene_call==j,]

### TO NOTE : much more SCV than SNV because the SCV table does not have a "cov_outlier" and therefore has a filter only on entropy/departure ###

for(i in sample){
  
  # read filtered SNV table 
  df <- read.table(paste0("SNV_", i, "_intra_filtered_SCG.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "") # filtered SNVs (intra-samples)
  
  # reads filtered SCV table
  df2 <- read.table(paste0("SCV_", i, "_intra_filtered_SCG.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "") # filtered SCVs (intra-samples)
  
  # create empty dataframes that will contain the scv generated by filtered snvs
  df3 <- data.frame()
  
  # For each gene caller id
  for(j in unique(df$corresponding_gene_call)){
    
    # list codon numbers and codon orders from the filtered SNVs table
    codon_numbers <- df[df$corresponding_gene_call==j, "codon_number"]
    codon_orders <- df[df$corresponding_gene_call==j, "codon_order_in_gene"]
    
    # for each codon order in gene caller id 
    for(k in 1:length(codon_orders)){
      
      # select SCV with same codon number and codon order as in the filtered SNV table
      dff <- df2[df2$codon_number == codon_numbers[k] & 
                   df2$codon_order_in_gene == codon_orders[k] &
                   df2$corresponding_gene_call==j,] 
      
      # bind it to main SCV table generated from SNV
      df3 <- rbind(df3, dff)
    }
  }
  
  # remove duplicated SCVs
  df3 <- unique(df3)
  
  # write table
  write.table(df3, paste0("SCV_", i, "_intra_filtered_from_filtered_SNVs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
  # assign tables for each sample
  assign(paste0("snv_", i), df)
  assign(paste0("scv_raw_", i), df2)
  assign(paste0("scv_", i), df3)
}


# Select SAAVs generated from SCVs generated themselves by filtered SNVs

# i <- "M11"
# j <- 20
# test <- df2[df2$corresponding_gene_call==j,]

for(i in sample){
  df2 <- eval(parse(text = paste0("scv_", i)))
  df3 <- read.table(paste0("SAAV_", i, "_intra_filtered_SCG.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "") # filtered SAAVs (intra-samples)
  
  df4 <- data.frame()
  
  for(j in unique(df2$corresponding_gene_call)){
    codon_numbers <- df2[df2$corresponding_gene_call==j, "codon_number"]
    codon_orders <- df2[df2$corresponding_gene_call==j, "codon_order_in_gene"]
    
    for(k in 1:length(codon_orders)){
      dff <- df3[df3$codon_number == codon_numbers[k] & 
                   df3$codon_order_in_gene == codon_orders[k] &
                   df3$corresponding_gene_call==j,] 
      df4 <- rbind(df4, dff)
    }
  }
  
  write.table(df4, paste0("SAAV_", i, "_intra_filtered_from_filtered_SCVs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  assign(paste0("saav_raw_", i), df3)
  assign(paste0("saav_", i), df4)
}
