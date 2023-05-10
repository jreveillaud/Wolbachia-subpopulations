setwd("/Volumes/Elements/Hiseq/metapangenomics/output/")

require(tidyverse)
require(openxlsx)

'%!in%' <- function(x,y)!('%in%'(x,y))

samples_loop <- c("O11", "M11", "O03", "O07", "O12")

s <- "O11"
# s <- "M11"


for (s in samples_loop) { 
  
  
  #### 1) SNV table and duplicates ####
  
  # Import SNV raw table with sample of interest as reference (to have all samples as SNVs)
  snv <- read.table(paste0("SNV/snv-tables/", s, "-as-reference/variability-", s, "-all.txt"), header=TRUE, fill=TRUE)
  
  # Filter SNVs by removing cov outliers and keep SNVs with 20% of entropy and departure from consensus
  snv_filtered <- snv[snv$cov_outlier_in_contig==0 & snv$cov_outlier_in_split==0 & snv$entropy>=0.2 & snv$departure_from_consensus>=0.2,]
  
  # Import SNV raw table from intra-samples in order to have the correct entry_id for SNV layers (later)
  snv2 <- read.table(paste0("SNV/snv-tables/intra-samples/variability_", s, ".txt"), header=TRUE, fill=TRUE)
  snv2_filtered <- snv2[snv2$cov_outlier_in_contig==0 & snv2$cov_outlier_in_split==0 & snv2$entropy>=0.2 & snv2$departure_from_consensus>=0.2,]
  
  # Add pos, split_name and entry_id from full SNVs table to SNVs with sample of interest using as reference in order to correct entry_id if needed
  snv_filtered <- snv_filtered %>% merge(snv2_filtered %>% select(c(pos, split_name, entry_id)), by=c("pos","split_name"), all.x=TRUE, suffixes=c("","2"))
  snv_filtered <- snv_filtered %>% select(c(pos, split_name, entry_id, entry_id2, everything()))
  colnames(snv_filtered)[3:4] <- c("entry_id_ref", "entry_id")
  
  # Select splits
  split_selection <- snv_filtered[, "split_name"] %>% unique() %>% as.character()
  
  # Filter and save snv with that are at the same position at least 2 times (and 3 times) in each split
  snv_dup <- data.frame() # all duplicated SNVs
  snv_dup_3 <- data.frame() # SNVs duplicated at least 3 times
  snv_dup_4 <- data.frame() # SNVs duplicated at least 4 times
  snv_dup_5 <- data.frame() # SNVs duplicated at least 5 times (in all the metagenomes)
  
  # Create columns to indicate if SNVs are shared at the same position and in which sample
  snv_filtered$snv_shared_in <- 0
  snv_filtered$nb_samples_where_snv_shared <- 0
  
  #i <- "Culex_O11_000000010761_split_00001"
  #j <- 148
  
  #i <- "Culex_O11_000000168902_split_00001"
  #j <- 112
  
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
  
  # Select filtered SNVs and shared SNVs occuring in the sample
  snv_filtered_sample <- snv_filtered[snv_filtered$sample_id==s,]
  snv_dup_sample <- snv_dup[snv_dup$sample_id==s,]
  
  # Annotate shared SNVs in the filtered table by merging it with the shared SNVs table (intra-sample level)
  snv_to_merge <- snv_filtered_sample %>% merge(snv_dup_sample %>% select(c(entry_id_ref, snv_shared_in, nb_samples_where_snv_shared)), 
                                                by="entry_id_ref", all.x = TRUE, suffixes=c("_to_remove",""))
  snv_to_merge <- snv_to_merge %>% select(-c(snv_shared_in_to_remove, nb_samples_where_snv_shared_to_remove))
  snv_to_merge <- snv_to_merge %>% select(c(pos, split_name, entry_id_ref, entry_id, everything()))
  snv_to_merge[is.na(snv_to_merge$snv_shared_in), c("snv_shared_in", "nb_samples_where_snv_shared")] <- "SNV not shared"
  
  # Annotate shared SNVs in the filtered table by merging it with the shared SNVs table (inter-sample level)
  snv_to_merge_inter <- snv_filtered %>% merge(snv_dup %>% select(c(entry_id_ref, snv_shared_in, nb_samples_where_snv_shared)), 
                                                by="entry_id_ref", all.x = TRUE, suffixes=c("_to_remove",""))
  snv_to_merge_inter <- snv_to_merge_inter %>% select(-c(snv_shared_in_to_remove, nb_samples_where_snv_shared_to_remove))
  snv_to_merge_inter <- snv_to_merge_inter %>% select(c(entry_id_ref, unique_pos_identifier, pos, pos_in_contig, contig_name, split_name, sample_id, everything()))
  snv_to_merge_inter[is.na(snv_to_merge_inter$snv_shared_in), c("snv_shared_in", "nb_samples_where_snv_shared")] <- "SNV not shared"
  colnames(snv_to_merge_inter)[1] <- "entry_id"
  
  
  
  #### 2) Merge SNV dup to SNV layer ####
  
  # Load SNV layer table
  snv_layer <- read.table(paste0("SNV/snv-tables/intra-samples/variability-", s, "-layers.txt"), header=TRUE, fill=TRUE)
  snv_layer <- snv_layer %>% select(c(pos, split_name, entry_id, everything()))
  
  # Edit label of SNV entry_id to merge
  snv_to_merge$entry_id <- paste0("p_", sprintf("%08d",snv_to_merge$entry_id))
  
  # Merge SNV layer table and shared SNV table
  snv_to_merge <- snv_to_merge %>% merge(snv_layer %>% select(c(entry_id, Bonneau, Bonneau_pct_alignment, Beckmann, Beckmann_pct_alignment,
                                                                MLST_WSP, MLST_WSP_pct_alignment)), by="entry_id")
  
  
  
  
  #### 3) Make wSCG table with shared wSCGs ####
  
  
  ### 3.1) wSCG-metagenomes ###
  
  # Load gene clusters from pangenomics with all metagenomes
  gene_clusters <- read.table("Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-all-metagenomes-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt", header=TRUE, fill=TRUE, sep="\t", quote = "")
  gene_clusters_draft <- gene_clusters
  
  # Select gene clusters from sample of the loop
  gene_clusters_sample <- gene_clusters[gene_clusters$genome_name==s,]
  
  # Edit colname of gene calls
  colnames(gene_clusters_sample)[5] <- "corresponding_gene_call"
  
  # Create empty dataframes that will contain the wSCGs, the wSCG but not present in the sample and the multi-copy genes 
  wSCG <- data.frame()
  wSCG_but_not_in_this_sample <- data.frame()
  multi_copy_genes <- data.frame()

  # List the unique gene clusters ids
  gene_cluster_ids <- gene_clusters$gene_cluster_id %>% as.character() %>% unique()

  # For test
  # cluster <- "GC_00000001"
  # cluster <- "GC_00001054"
  # cluster <- "GC_00001020"
  
  # For each gene cluster (GC)
  for(cluster in gene_cluster_ids){
    
    # Select it
    df <- gene_clusters[gene_clusters$gene_cluster_id==cluster,]
    
    # Create a column that will indicate if the GC is a wSCG
    df$wSCG_detailed <- "wSCG"
    genome_names <- df$genome_name %>% as.character() %>% unique()
    
    # Identify and store wSCGs in the wSCG table :
    
    # If there is only one gene in the GC and the number of gene is inferior or equal to the number of samples
    if(nrow(df[df$genome_name==s,])==1 && nrow(df)<=length(genome_names)){
      
      # Create a column that will include the wSCG and the samples with wich they are shared
      df2 <- df
      df2$wSCG_shared_in_metagenomes <- ""
      
      # If there is only one gene in the GC, tag it to "Only in sample"
      if(nrow(df2)==1){
        df2$wSCG_detailed <- paste0("Only in ", s)
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_metagenomes"] <- shared
        
      # If there is two genes in the GC, tag it to "Shared by sample and 1 MAG"
      } else if(nrow(df2)==2){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 1 MAG")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_metagenomes"] <- shared
        
      # If there is three genes in the GC, tag it to "Shared by sample and 2 MAGs"
      } else if(nrow(df2)==3){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 2 MAGs")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_metagenomes"] <- shared
        
      # If there is four genes in the GC, tag it to "Shared by sample and 3 MAGs"
      } else if(nrow(df2)==4){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 3 MAGs")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_metagenomes"] <- shared
        
      # If there is five genes in the GC, tag it to "Shared by sample and 4 MAGs"
      } else if(nrow(df2)==5){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 4 MAGs")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_metagenomes"] <- shared
      }
      
      # Add all the wSCGs identified in the wSCG table
      wSCG <- wSCG %>% rbind(df2)
    }
    
    # If there is no gene in the GC from the sample, tag the GC as "wSCG but not present in sample"
    if(nrow(df[df$genome_name==s,])==0 && nrow(df)<=length(genome_names)){
      df3 <- df
      if(nrow(df)!=1){
        df3$wSCG_detailed <- paste0("wSCG but not present in ", s)
        
        shared <- paste(df3[, "genome_name"] %>% as.character(), collapse=", ")
        df3[, "wSCG_shared_in_metagenomes"] <- shared
        
        #df3$wSCG_shared_in_metagenomes <- ""
      } else {
        df3$wSCG_detailed <- ""
        df3$wSCG_shared_in_metagenomes <- ""
      }
      wSCG_but_not_in_this_sample <- wSCG_but_not_in_this_sample %>% rbind(df3)
    }
    
    # If there is at least 2 genes in the GC from one of the sample, tag the GC as "multicopy genes"
    if(nrow(df[df$genome_name=="O11",])>1 | nrow(df[df$genome_name=="M11",])>1 | nrow(df[df$genome_name=="O03",])>1 |
       nrow(df[df$genome_name=="O07",])>1 | nrow(df[df$genome_name=="O12",])>1){
      df4 <- df
      
      # If there is no gene in the GC from the sample and genes in multi-copy from another sample, tag the GC as "multicopy genes but not present in sample"
      if(nrow(df[df$genome_name==s,])==0){
        df4$wSCG_detailed <- paste0("multicopy genes but not present in ", s)
        
        shared <- paste(df4[, "genome_name"] %>% as.character(), collapse=", ")
        df4[, "wSCG_shared_in_metagenomes"] <- shared
        
        #df4$wSCG_shared_in_metagenomes <- df$genome_name
        #df4$wSCG_shared_in_metagenomes <- ""
        
      # If there is genes in the GC from the sample, tag the GC as "multicopy genes"
      } else if(nrow(df[df$genome_name==s,])!=0){
        df4$wSCG_detailed <- "multicopy genes"
        
        shared <- paste(df4[, "genome_name"] %>% as.character(), collapse=", ")
        df4[, "wSCG_shared_in_metagenomes"] <- shared
        
        #df4$wSCG_shared_in_metagenomes <- df$genome_name
        #df4$wSCG_shared_in_metagenomes <- ""
      }
      multi_copy_genes <- multi_copy_genes %>% rbind(df4)
      #multi_copy_genes <- multi_copy_genes[duplicated(multi_copy_genes$unique_id),]
      #multi_copy_genes <- multi_copy_genes %>% unique("unique_id")
    }
    
  }
  
  # List of unique wSCGs identified
  test1 <- unique(wSCG$gene_cluster_id %>% droplevels()) %>% as.character()
  
  # List of unique wSCGs identified but not present in the sample
  test2 <- unique(wSCG_but_not_in_this_sample$gene_cluster_id %>% droplevels()) %>% as.character()
  
  # List of unique multicopy genes identified
  test3 <- unique(multi_copy_genes$gene_cluster_id %>% droplevels()) %>% as.character()
  
  # Merge all the wSCG, wSCG but not in the sample and multicopy genes in a single table
  wSCG_final <- rbind(wSCG, multi_copy_genes, wSCG_but_not_in_this_sample)
  
  # Check if we have the good total of unique GC in the final table
  if(unique(wSCG_final$gene_cluster_id %>% droplevels()) %>% as.character() %>% length() == gene_cluster_ids %>% length){
    print(paste0("The total of GCs in wSCG table is OK for ", s, " !"))
  } else {
    print(paste0("We have a problem with the total of GCs in ", s))
  }
    
  
  
  ### 3.2) wSCG-references ###
  
  # Load gene clusters from pangenomics with the MAG corresponding to the sample and the 3 Wolbachia reference genomes
  gene_clusters_ref <- read.table(paste0("Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-", s, "-ref-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "")
  gene_clusters_ref_draft <- gene_clusters_ref
  
  # Select gene clusters from sample of the loop
  gene_cluster_ids <- gene_clusters_ref$gene_cluster_id %>% as.character() %>% unique()
  
  # Create empty dataframes that will contain the wSCGs, the wSCG but not present in the sample and the multi-copy genes 
  wSCG_ref <- data.frame()
  wSCG_ref_not_in_this_sample <- data.frame()
  multi_copy_genes_ref <- data.frame()
  
  # Identify and store wSCGs in the wSCG ref table :
  
  cluster <- "GC_00001151"
  # For each gene cluster (GC)
  for(cluster in gene_cluster_ids){
    
    # Select it
    df <- gene_clusters_ref[gene_clusters_ref$gene_cluster_id==cluster,]
    
    # Create a column that will indicate if the GC is a wSCG
    df$wSCG_detailed <- "wSCG"
    genome_names <- df$genome_name %>% as.character() %>% unique()
    
    # If there is only one gene in the GC and the number of gene is inferior or equal to the number of sample + reference genomes
    if(nrow(df[df$genome_name==s,])==1 && nrow(df)<=length(genome_names)){
      
      # Create a column that will include the wSCG and the samples with wich they are shared
      df2 <- df
      df2$wSCG_shared_in_references <- ""
      
      # If there is only one gene in the GC, tag it to "Only in sample"
      if(nrow(df2)==1){
        df2$wSCG_detailed <- paste0("Only in ", s)
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_references"] <- shared
        
      # If there is two genes in the GC, tag it to "Shared by sample and 1 reference genome"
      } else if(nrow(df2)==2){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 1 reference genome")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_references"] <- shared
        
      # If there is three genes in the GC, tag it to "Shared by sample and 2 reference genomes"
      } else if(nrow(df2)==3){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 2 reference genomes")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_references"] <- shared
        
      # If there is three genes in the GC, tag it to "Shared by sample and 3 reference genomes"
      } else if(nrow(df2)==4){
        df2$wSCG_detailed <- paste0("Shared by ", s, " and 3 reference genomes")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "wSCG_shared_in_references"] <- shared
      }
      
      # Add all the wSCGs identified in the wSCG table
      wSCG_ref <- wSCG_ref %>% rbind(df2)
    }
    
    # If there is no gene in the GC from the sample, tag the GC as "wSCG but not present in sample"
    if(nrow(df[df$genome_name==s,])==0 && nrow(df)<=length(genome_names)){
      df3 <- df
      if(nrow(df)!=1){
        df3$wSCG_detailed <- paste0("wSCG but not present in ", s)
        
        shared <- paste(df3[, "genome_name"] %>% as.character(), collapse=", ")
        df3[, "wSCG_shared_in_references"] <- shared
        
        #df3$wSCG_shared_in_references<- ""
      } else {
        df3$wSCG_detailed <- ""
        df3$wSCG_shared_in_references <- ""
      }
      wSCG_ref_not_in_this_sample <- wSCG_ref_not_in_this_sample %>% rbind(df3)
    }
    
    # If there is at least 2 genes in the GC from one of the sample, tag the GC as "multicopy genes"
    if(nrow(df[df$genome_name==s,])>1 | nrow(df[df$genome_name=="wPipJHB",])>1 | nrow(df[df$genome_name=="wPipMol",])>1 |
       nrow(df[df$genome_name=="wPipPel",])>1){
      df4 <- df
      
      # If there is no gene in the GC from the sample and genes in multi-copy from one reference genome, tag the GC as "multicopy genes but not present in sample"
      if(nrow(df[df$genome_name==s,])==0){
        df4$wSCG_detailed <- paste0("multicopy but not present in ", s)
        shared <- paste(df4[, "genome_name"] %>% as.character(), collapse=", ")
        df4[, "wSCG_shared_in_references"] <- shared
        
        #df4$wSCG_shared_in_references <- ""
        
      # If there is genes in the GC from the sample, tag the GC as "multicopy genes"
      } else if(nrow(df[df$genome_name==s,])!=0){
        df4$wSCG_detailed <- "multicopy genes"
        shared <- paste(df4[, "genome_name"] %>% as.character(), collapse=", ")
        df4[, "wSCG_shared_in_references"] <- shared
        #df4$wSCG_shared_in_references <- ""
      }
      multi_copy_genes_ref <- multi_copy_genes_ref %>% rbind(df4)
      #multi_copy_genes_ref <- multi_copy_genes_ref[duplicated(multi_copy_genes_ref$unique_id),]
      #multi_copy_genes_ref <- multi_copy_genes_ref %>% unique("unique_id")
    }
  }
  
  # List of unique wSCG_ref identified
  test1 <- unique(wSCG_ref$gene_cluster_id %>% droplevels()) %>% as.character()
  
  # List of unique wSCG_ref identified but not present in the sample
  test2 <- unique(wSCG_ref_not_in_this_sample$gene_cluster_id %>% droplevels()) %>% as.character()
  
  # List of unique multicopy genes identified
  test3 <- unique(multi_copy_genes_ref$gene_cluster_id %>% droplevels()) %>% as.character()
  
  # Merge all the wSCG, wSCG but not in the sample and multicopy genes in a single table
  wSCG_ref_final <- rbind(wSCG_ref, multi_copy_genes_ref, wSCG_ref_not_in_this_sample)
  
  # Check if we have the good total of unique GC in the final table
  if(unique(wSCG_ref_final$gene_cluster_id %>% droplevels()) %>% as.character() %>% length() == gene_cluster_ids %>% length){
    print(paste0("The total of GCs in wSCG-ref table is OK for ", s, " !"))
  } else {
    print(paste0("We have a problem with the total of GCs in ", s))
  }
  
  
  ### 3.3) Make proper GC tables including wSCG-MAGs and wSCG-ref infos ###
  
  # GC for MAGs
  gene_clusters <- wSCG_final[wSCG_final$genome_name==s,]
  gene_clusters <- gene_clusters %>% select(c(gene_callers_id, wSCG_detailed, c(13:27)))
  colnames(gene_clusters)[1:2] <- c("corresponding_gene_call", "wSCG")
  gene_clusters = gene_clusters[!duplicated(gene_clusters$corresponding_gene_call),]
  
  # GC for reference genomes
  gene_clusters_ref <- wSCG_ref_final[wSCG_ref_final$genome_name==s,]
  gene_clusters_ref <- gene_clusters_ref %>% select(c(gene_callers_id, wSCG_detailed, c(13:27)))
  colnames(gene_clusters_ref)[1:2] <- c("corresponding_gene_call", "wSCG")
  gene_clusters_ref = gene_clusters_ref[!duplicated(gene_clusters_ref$corresponding_gene_call),]    
  
  
  # Nice tables for draft
  gene_clusters_draft <- gene_clusters_draft %>% merge(wSCG_final %>% select(c(unique_id, wSCG_detailed, wSCG_shared_in_metagenomes)), 
                                                       by="unique_id", all.x=TRUE)
  colnames(gene_clusters_draft)[26:27] <- c(paste0("wSCG_", s), paste0("wSCG_", s, "_shared_in_metagenomes"))
  
  gene_clusters_ref_draft <- gene_clusters_ref_draft %>% merge(wSCG_ref_final %>% select(c(unique_id, wSCG_detailed, wSCG_shared_in_references)), 
                                                               by="unique_id", all.x=TRUE)
  colnames(gene_clusters_ref_draft)[26] <- paste0("wSCG_references")
  
  gene_clusters_draft2 <- gene_clusters_draft %>% select(c(gene_cluster_id, SCG, paste0("wSCG_", s), paste0("wSCG_", s, "_shared_in_metagenomes")))
  gene_clusters_draft2 <- unique(gene_clusters_draft2)
  
  gene_clusters_ref_draft2 <- gene_clusters_ref_draft %>% select(c(gene_cluster_id, SCG, wSCG_references, wSCG_shared_in_references))
  gene_clusters_ref_draft2 <- unique(gene_clusters_ref_draft2)  
  
  write.table(gene_clusters_draft2, paste0("SNV/snv-tables/intra-samples/GC-wSCG-",s, "-metagenomes-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(gene_clusters_ref_draft2, paste0("SNV/snv-tables/intra-samples/GC-wSCG-",s, "-references-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
  
  
  #### 4) Merge wSCG informations ####
  merged_gene_clusters <- gene_clusters %>% merge(gene_clusters_ref %>% select(c(corresponding_gene_call, wSCG, wSCG_shared_in_references)), 
                                                  by="corresponding_gene_call", all.x = TRUE, all.y=FALSE)
  
  merged_gene_clusters <- merged_gene_clusters %>% select(c(corresponding_gene_call, wSCG.x, wSCG.y, everything()))
  colnames(merged_gene_clusters)[2:3] <- c(paste0("wSCG_", s), "wSCG_ref")
  
  
  
  #### 5) Merge SCG to variability tables ####
  
  #### 5.1) Add functions to gene in SNV table
  
  snv_to_merge <- snv_to_merge %>% merge(gene_clusters_sample %>% select(c("corresponding_gene_call", 13:25)), by="corresponding_gene_call", all.x=TRUE)
  
  
  
  #### 5.1) Merge SCV and wSCG tables
  
  scv_filtered <- read.table(paste0("SNV/snv-tables/SCVs/SCV_variability_", s, "_intra-samples_filtered_from_filtered_SNVs.txt"), header=TRUE, fill=TRUE)
  
  test <- unique(snv_to_merge %>% select(c(corresponding_gene_call, in_coding_gene_call, in_noncoding_gene_call, Bonneau, Bonneau_pct_alignment, MLST_WSP, MLST_WSP_pct_alignment)))
  scv_filtered <- scv_filtered %>% merge(test, by="corresponding_gene_call", all.x = TRUE)
  
  scv_scg <- scv_filtered %>% merge(merged_gene_clusters %>% select(c("corresponding_gene_call",
                                                                      "wSCG_shared_in_metagenomes", 
                                                                      "wSCG_shared_in_references",
                                                                      "wSCG_ref", 
                                                                      paste0("wSCG_", s))), 
                                    by="corresponding_gene_call", all.x=TRUE)
  
  
  
  #### 5.2) Merge SAAV and wSCG tables
  
  saav_filtered <- read.table(paste0("SNV/snv-tables/SAAVs/SAAV_variability_", s, "_intra-samples_filtered_from_filtered_SCVs.txt"), header=TRUE, fill=TRUE)
  
  test <- unique(snv_to_merge %>% select(c(corresponding_gene_call, in_coding_gene_call, in_noncoding_gene_call, Bonneau, Bonneau_pct_alignment, MLST_WSP, MLST_WSP_pct_alignment)))
  saav_filtered <- saav_filtered %>% merge(test, by="corresponding_gene_call", all.x = TRUE)
  
  saav_scg <- saav_filtered %>% merge(merged_gene_clusters %>% select(c("corresponding_gene_call",
                                                                        "wSCG_shared_in_metagenomes", 
                                                                        "wSCG_shared_in_references",
                                                                        "wSCG_ref", 
                                                                        paste0("wSCG_", s))), 
                                      by="corresponding_gene_call", all.x=TRUE)
  
  
  
  #### 5.3) Merge SNV and wSCG tables ####
  
  snv_scg <- snv_to_merge %>% merge(merged_gene_clusters %>% select(c("corresponding_gene_call", 
                                                                      "wSCG_shared_in_metagenomes", 
                                                                      "wSCG_shared_in_references",
                                                                      "wSCG_ref", 
                                                                      paste0("wSCG_", s))), 
                                    by="corresponding_gene_call", all.x=TRUE)
  
  # Tag SNV that not occuring in gene as "out of gene"
  snv_scg[is.na(snv_scg[,paste0("wSCG_", s)]), paste0("wSCG_", s)] <- "out of gene"
  snv_scg$wSCG_ref[is.na(snv_scg$wSCG_ref)] <- "out of gene"
  snv_scg$wSCG_shared_in_metagenomes[is.na(snv_scg$wSCG_shared_in_metagenomes)] <- "out of gene"
  snv_scg$wSCG_shared_in_references[is.na(snv_scg$wSCG_shared_in_references)] <- "out of gene"
  
  
  snv_scg_reduced <- snv_scg %>% select(c(entry_id,
                                          pos,
                                          pos_in_contig,
                                          split_name,
                                          contig_name,
                                          in_coding_gene_call,
                                          corresponding_gene_call,
                                          paste0("wSCG_", s), 
                                          wSCG_shared_in_metagenomes,
                                          wSCG_ref,
                                          wSCG_shared_in_references,
                                          nb_samples_where_snv_shared,
                                          snv_shared_in,
                                          Bonneau,
                                          Bonneau_pct_alignment,
                                          Beckmann, 
                                          Beckmann_pct_alignment,
                                          MLST_WSP,
                                          MLST_WSP_pct_alignment,
                                          COG20_PATHWAY_ACC,
                                          COG20_PATHWAY,
                                          COG20_FUNCTION_ACC,
                                          COG20_FUNCTION,
                                          COG20_CATEGORY_ACC,
                                          COG20_CATEGORY,
                                          KOfam_ACC,
                                          KOfam,
                                          KEGG_Module_ACC,
                                          KEGG_Module,
                                          KEGG_Class_ACC,
                                          KEGG_Class
  ))
  
  snv_scg_complete <- snv_scg %>% select(c(entry_id,
                                           entry_id_ref,
                                           pos,
                                           pos_in_contig,
                                           split_name,
                                           contig_name,
                                           in_coding_gene_call,
                                           corresponding_gene_call,
                                           paste0("wSCG_", s), 
                                           wSCG_shared_in_metagenomes,
                                           wSCG_ref,
                                           wSCG_shared_in_references,
                                           nb_samples_where_snv_shared,
                                           snv_shared_in,
                                           Bonneau,
                                           Bonneau_pct_alignment,
                                           Beckmann, 
                                           Beckmann_pct_alignment,
                                           MLST_WSP,
                                           MLST_WSP_pct_alignment,
                                           COG20_PATHWAY_ACC,
                                           COG20_PATHWAY,
                                           COG20_FUNCTION_ACC,
                                           COG20_FUNCTION,
                                           COG20_CATEGORY_ACC,
                                           COG20_CATEGORY,
                                           KOfam_ACC,
                                           KOfam,
                                           KEGG_Module_ACC,
                                           KEGG_Module,
                                           KEGG_Class_ACC,
                                           KEGG_Class,
                                           everything()
  ))
  
  snv_scg_reduced$SNV_shared_in_metagenomes <- "SNV not shared"
  snv_scg_reduced[snv_scg_reduced$nb_samples_where_snv_shared!="SNV not shared", "SNV_shared_in_metagenomes"] <- "SNV shared"
  
  snv_scg_complete$SNV_shared_in_metagenomes <- "SNV not shared"
  snv_scg_complete[snv_scg_complete$nb_samples_where_snv_shared!="SNV not shared", "SNV_shared_in_metagenomes"] <- "SNV shared"
  
  
  #### 6) Create layers for SNV interactive figure  ####
  
  # Select infos to add to the SNV interactive plot
  snv_scg_layer <- snv_scg_complete %>% select(c(entry_id, paste0("wSCG_", s), wSCG_shared_in_metagenomes, 
                                                 wSCG_ref, wSCG_shared_in_references, 
                                                 nb_samples_where_snv_shared, snv_shared_in, SNV_shared_in_metagenomes,
                                                 in_coding_gene_call))
  
  # Give the number 3 as label for SNV occuring outside genes
  snv_scg_layer[snv_scg_layer[, paste0("wSCG_", s)]=="out of gene", "in_coding_gene_call"] <- 3
  snv_scg_layer[snv_scg_layer$wSCG_ref=="out of gene", "in_coding_gene_call"] <- 3
  
  # Change name of column "in_coding_gene_call" to avoid conflicts with existing layer in SNV interactive plot
  colnames(snv_scg_layer)[ncol(snv_scg_layer)] <- "genes"

  # Change label of coding gene call (1 -> A), non-coding gene call (0 -> B) and not in gene (3 -> C) to easily order this layer in the plot
  snv_scg_layer[snv_scg_layer$genes==1, "genes"] <- "A"
  snv_scg_layer[snv_scg_layer$genes==0, "genes"] <- "B"
  snv_scg_layer[snv_scg_layer$genes==3, "genes"] <- "C"


  #### 6) Make the table for the draft (with reduced and coherent columns)

  selected_cols <- c("entry_id",
                     "unique_pos_identifier",
                     "pos",
                     "pos_in_contig",
                     "split_name",
                     "sample_id",
                     "corresponding_gene_call",
                     "in_noncoding_gene_call",
                     "in_coding_gene_call",
                     "base_pos_in_codon",
                     "codon_order_in_gene",
                     "codon_number",
                     "gene_length",
                     "coverage",
                     "cov_outlier_in_split",
                     "cov_outlier_in_contig",
                     "A", "C", "G", "N", "T",
                     "reference",
                     "consensus",
                     "competing_nts",
                     "departure_from_reference",
                     "departure_from_consensus",
                     "n2n1ratio",
                     "entropy",
                     paste0("wSCG_", s),
                     "wSCG_shared_in_metagenomes",
                     "wSCG_ref",
                     "wSCG_shared_in_references",
                     "nb_samples_where_snv_shared",
                     "snv_shared_in",
                     "SNV_shared_in_metagenomes",
                     "Bonneau",
                     "Bonneau_pct_alignment",
                     "MLST_WSP",
                     "MLST_WSP_pct_alignment",
                     "COG20_PATHWAY_ACC",
                     "COG20_PATHWAY",
                     "COG20_FUNCTION_ACC",
                     "COG20_FUNCTION",
                     "COG20_CATEGORY_ACC",
                     "COG20_CATEGORY",
                     "KOfam_ACC",
                     "KOfam",
                     "KEGG_Module_ACC",
                     "KEGG_Module",
                     "KEGG_Class_ACC",
                     "KEGG_Class",
                     "aa_sequence"
  )

  snv_scg_for_draft <- snv_scg_complete %>% select(all_of(selected_cols))
  colnames(snv_scg_for_draft)[33] <- "Number_of_metagenomes_where_SNV_are_shared"
  # colnames(snv_scg_for_draft)[34] <- "SNV_shared_in_metagenomes"
  colnames(snv_scg_for_draft)[36] <- "cid_genes_from_Bonneau_et_al._2018"
  colnames(snv_scg_for_draft)[38] <- "MLST_and_wsp_genes_from_PubMLST"
  
  # Write table 
  write.table(snv_to_merge_inter, paste0("SNV/snv-tables/", s, "-as-reference/variability-", s, "-snv-shared.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_scg_for_draft, paste0("SNV/snv-tables/intra-samples/SNV-wSCG-",s, "-complete.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_scg_layer, paste0("SNV/snv-tables/intra-samples/variability-",s, "-layers-wSCG.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(scv_scg, paste0("SNV/snv-tables/intra-samples/SCV-",s , "-filtered-SCG-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(saav_scg, paste0("SNV/snv-tables/intra-samples/SAAV-",s , "-filtered-SCG-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
#   
#   
}
