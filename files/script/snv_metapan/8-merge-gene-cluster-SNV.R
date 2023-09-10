setwd("../output/")

require(tidyverse)
require(openxlsx)

'%!in%' <- function(x,y)!('%in%'(x,y))

samples_loop <- c("O11", "M11", "O03", "O07", "O12")

# s <- "O11"
# s <- "M11"

SNV_dup_SCG_all <- data.frame()

# IDENTIFY SNV SHARED IN METAGENOMES AND EGGS #

for (s in samples_loop) { 
  
  # Import SNV raw table with sample of interest as reference (to have all samples as SNVs)
  snv <- read.table(paste0("variability/SNV_", s, "_filtered_SCG.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "")
  snv_filtered <- snv
  
  # Select splits
  split_selection <- snv_filtered[, "split_name"] %>% unique() %>% as.character()
  
  # Filter and save snv with that are at the same position at least 2 times (and 3 times) in each split
  snv_dup <- data.frame() # all duplicated SNVs
  
  # Create columns to indicate if SNVs are shared at the same position and in which sample
  snv_filtered$snv_shared_in <- 0
  snv_filtered$nb_samples_where_snv_shared <- 0
  
  # For each split
  for(i in split_selection){
    
    # Select filtered SNVs occuring in the selected split
    sbs <- snv_filtered[snv_filtered$split_name==i,]
    
    # For each unique position
    for(j in unique(sbs$pos)){
      
      # Select SNVs occuring at the same position in the selected split
      sbs_2 <- sbs[sbs$pos==j,]
      
      # Select SNVs occuring at least two times at the same position
      sbs_dup <- sbs_2[sbs_2$pos %in% names(which(table(sbs_2$pos) >= 2)),]
      
      # If at least 2 SNVs are shared at the same positions, annotate them to "shared"
      if(nrow(sbs_dup)!=0){
        shared <- paste(sbs_dup[, "sample_id"] %>% as.character(), collapse=", ")
        sbs_dup[, "snv_shared_in"] <- shared
        sbs_dup$nb_samples_where_snv_shared <- nrow(sbs_dup)
      }
      
      # Store shared SNVs by adding them into empty dataframe
      snv_dup <- snv_dup %>% rbind(sbs_dup)
    }
  }
  
  # Select shared SNV from SCG
  snv_dup_SCG <- snv_dup[!is.na(snv_dup$SCG), ]
  snv_dup_SCG <- snv_dup_SCG[snv_dup_SCG$SCG==1, ]
  
  snv_dup_SCG_reduced <- snv_dup_SCG %>% select(c(corresponding_gene_call, gene_cluster_id, SCG, COG20_FUNCTION, entry_id, pos, split_name, 
                                                  sample_id, in_coding_gene_call, base_pos_in_codon, coverage,departure_from_consensus))
  
  snv_dup_SCG_reduced$genome <- s
  
  snv_dup_SCG_reduced <- snv_dup_SCG_reduced %>% select(c(genome, everything()))
  
  SNV_dup_SCG_all <- SNV_dup_SCG_all %>% rbind(snv_dup_SCG_reduced)
  
  # Select splits, gene and position where SNVs fall
  splits_SCG <- snv_dup_SCG[, "split_name"] %>% unique() %>% as.character()
  genes_SCG <- snv_dup_SCG[, "corresponding_gene_call"] %>% unique() %>% as.character()
  pos_SCG <- snv_dup_SCG[, "pos_in_contig"] %>% unique() %>% as.character()
  
  
  # EGGS 
  
  # import raw variability table including egg metagenomes
  snv_eggs <- read.table(paste0("eggs/variability/SNV_", s, "_raw_SCG_eggs.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "")
  
  # select SNV occuring only in split, gene, pos from filtered SNV table
  snv_eggs_select <- snv_eggs[snv_eggs$split_name %in% splits_SCG, ]
  snv_eggs_select <- snv_eggs_select[snv_eggs_select$corresponding_gene_call %in% genes_SCG, ]
  snv_eggs_select <- snv_eggs_select[snv_eggs_select$pos_in_contig %in% pos_SCG, ]
  
  # filter table
  snv_eggs_select_f <- snv_eggs_select[snv_eggs_select$entropy>=0.2 & snv_eggs_select$departure_from_consensus>=0.2,]
  
  # create columns indicating shared SNV from egg table
  snv_dup_eggs <- data.frame()
  split_selection <- splits_SCG
  snv_filtered <- snv_eggs_select
  
  # For each split
  for(i in split_selection){
    
    # Select filtered SNVs occuring in the selected split
    sbs <- snv_filtered[snv_filtered$split_name==i,]
    
    # For each unique position
    for(j in unique(sbs$pos)){
      
      # Select SNVs occuring at the same position in the selected split
      sbs_2 <- sbs[sbs$pos==j,]
      
      # Select SNVs occuring at least two times at the same position
      sbs_dup <- sbs_2[sbs_2$pos %in% names(which(table(sbs_2$pos) >= 2)),]
      
      # If at least 2 SNVs are shared at the same positions, annotate them to "shared"
      if(nrow(sbs_dup)!=0){
        shared <- paste(sbs_dup[, "sample_id"] %>% as.character(), collapse=", ")
        sbs_dup[, "snv_shared_in"] <- shared
        sbs_dup$nb_samples_where_snv_shared <- nrow(sbs_dup)
      }
      
      # Store shared SNVs by adding them into empty dataframe
      snv_dup_eggs <- snv_dup_eggs %>% rbind(sbs_dup)
    }
  }
  
  
  # Select unique position from the SNV tables
  snv_dup_SCG_pos <- snv_dup_SCG[!duplicated(snv_dup_SCG$pos_in_contig),]
  snv_dup_eggs_pos <- snv_dup_eggs[!duplicated(snv_dup_eggs$pos_in_contig),]

  # Select columns of interest
  
  snv_final_SCG <- snv_dup_SCG_pos %>% select(c(
    "entry_id",
    "pos",
    "pos_in_contig",
    "split_name",
    "corresponding_gene_call",
    "in_coding_gene_call",
    "COG20_FUNCTION",
    "snv_shared_in",
    "nb_samples_where_snv_shared",
    "base_pos_in_codon",
    "codon_order_in_gene",
  ))
  
  snv_final_eggs_SCG <- snv_dup_eggs_pos %>% select(c(
    "entry_id",
    "pos",
    "pos_in_contig",
    "split_name",
    "corresponding_gene_call",
    "in_coding_gene_call",
    "COG20_FUNCTION",
    "snv_shared_in",
    "nb_samples_where_snv_shared",
    "base_pos_in_codon",
    "codon_order_in_gene",
  ))
  
  write.table(snv_dup, paste0("shared_snv_stats/SNV_", s, "_dup.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_dup_SCG, paste0("shared_snv_stats/SNV_SCG_", s, "_dup.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_dup_SCG_reduced, paste0("shared_snv_stats/SNV_SCG_", s, "_dup_reduced.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_final_SCG, paste0("shared_snv_stats/SNV_SCG_", s, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_final_eggs_SCG, paste0("shared_snv_stats/SNV_SCG_EGGS_", s, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
