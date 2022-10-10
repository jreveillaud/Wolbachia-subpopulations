setwd("/Volumes/Elements/Hiseq/metapangenomics/output/")

require(tidyverse)
require(openxlsx)


samples_loop <- c("O11", "M11", "O03", "O07", "O12")

#s <- "O11"
#s <- "M11"

for (s in samples_loop) { 
  
  #### 1) SNV table and duplicates ####
  
  # Import SNV raw table with sample of interest as reference (to have all samples as SNVs)
  snv <- read.table(paste0("SNV/snv-tables/", s, "-as-reference/variability-", s, "-all.txt"), header=TRUE, fill=TRUE)
  snv_filtered <- snv[snv$cov_outlier_in_contig==0 & snv$cov_outlier_in_split==0 & snv$entropy>=0.2 & snv$departure_from_consensus>=0.2,]
  
  # Import SNV raw table from intra-MAGs in order to have the correct entry_id for SNV layers (later)
  snv2 <- read.table(paste0("SNV/snv-tables/intra-MAGs/variability_", s, ".txt"), header=TRUE, fill=TRUE)
  snv2_filtered <- snv2[snv2$cov_outlier_in_contig==0 & snv2$cov_outlier_in_split==0 & snv2$entropy>=0.2 & snv2$departure_from_consensus>=0.2,]
  
  snv_filtered <- snv_filtered %>% merge(snv2_filtered %>% select(c(pos, split_name, entry_id)), by=c("pos","split_name"), all.x=TRUE, suffixes=c("","2"))
  snv_filtered <- snv_filtered %>% select(c(pos, split_name, entry_id, entry_id2, everything()))
  colnames(snv_filtered)[3:4] <- c("entry_id_ref", "entry_id")
  
  # Select only gene that are SSG in the snv table 
  split_selection <- snv_filtered[, "split_name"] %>% unique() %>% as.character()
  
  # Filter and save snv with that are at the same position at least 2 times (and 3 times) in each split
  snv_dup <- data.frame()
  snv_dup_3 <- data.frame()
  snv_dup_4 <- data.frame()
  snv_dup_5 <- data.frame()
  
  snv_filtered$snv_shared_in <- 0
  snv_filtered$nb_samples_where_snv_shared <- 0
  
  #i <- "Culex_O11_000000010761_split_00001"
  #j <- 148
  
  #i <- "Culex_O11_000000168902_split_00001"
  #j <- 112
  
  for(i in split_selection){
    sbs <- snv_filtered[snv_filtered$split_name==i,]
    
    for(j in unique(sbs$pos)){
      sbs_2 <- sbs[sbs$pos==j,]
      
      sbs_dup <- sbs_2[sbs_2$pos %in% names(which(table(sbs_2$pos) >= 2)),]
      
      if(nrow(sbs_dup)!=0){
        shared <- paste(sbs_dup[, "sample_id"] %>% as.character(), collapse=", ")
        sbs_dup[, "snv_shared_in"] <- shared
        sbs_dup$nb_samples_where_snv_shared <- nrow(sbs_dup)
      }
      
      snv_dup <- snv_dup %>% rbind(sbs_dup)
      
      if(nrow(sbs_dup)>=3){
        snv_dup_3 <- snv_dup_3 %>% rbind(sbs_dup)
      }
      
      if(nrow(sbs_dup)>=4){
        snv_dup_4 <- snv_dup_4 %>% rbind(sbs_dup)
      }
      
      if(nrow(sbs_dup)>=5){
        snv_dup_5 <- snv_dup_5 %>% rbind(sbs_dup)
      }
    }
  }
  
  snv_filtered_sample <- snv_filtered[snv_filtered$sample_id==s,]
  snv_dup_sample <- snv_dup[snv_dup$sample_id==s,]
  snv_to_merge <- snv_filtered_sample %>% merge(snv_dup_sample %>% select(c(entry_id_ref, snv_shared_in, nb_samples_where_snv_shared)), 
                                                by="entry_id_ref", all.x = TRUE, suffixes=c("_to_remove",""))
  snv_to_merge <- snv_to_merge %>% select(-c(snv_shared_in_to_remove, nb_samples_where_snv_shared_to_remove))
  snv_to_merge <- snv_to_merge %>% select(c(pos, split_name, entry_id_ref, entry_id, everything()))
  
  
  
  #### 2) Merge SNV dup to SNV layer ####
  
  snv_layer <- read.table(paste0("SNV/snv-tables/intra-MAGs/variability-", s, "-layers.txt"), header=TRUE, fill=TRUE)
  snv_layer <- snv_layer %>% select(c(pos, split_name, entry_id, everything()))
  
  snv_to_merge$entry_id <- paste0("p_", sprintf("%08d",snv_to_merge$entry_id))
  snv_to_merge <- snv_to_merge %>% merge(snv_layer %>% select(c(entry_id, Bonneau, Bonneau_pct_alignment, Beckmann, Beckmann_pct_alignment,
                                                                MLST_WSP, MLST_WSP_pct_alignment)), by="entry_id")
  
  
  
  #### 3) Add functions to gene in SNV table ####
  
  # Import gene cluster table (pangenomics)
  gene_clusters <- read.table("Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-all-metagenomes-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt", header=TRUE, fill=TRUE, sep="\t", quote = "")
  gene_clusters_draft <- gene_clusters
  gene_clusters_sample <- gene_clusters[gene_clusters$genome_name==s,]
  colnames(gene_clusters_sample)[5] <- "corresponding_gene_call"
  snv_to_merge <- snv_to_merge %>% merge(gene_clusters_sample %>% select(c("corresponding_gene_call", 13:25)), by="corresponding_gene_call", all.x=TRUE)
  
  
  #### 4) Make SSG table with shared SSGs ####
  
#  gene_clusters_ref <- read.table(paste0("Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-", s, "-ref-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "")
  
  # Select gene cluster with at least the sample of interest as unique gene within each gene cluster
  gene_cluster_filt <- data.frame()
  main_sample <- s
  gene_cluster_ids <- gene_clusters$gene_cluster_id %>% as.character() %>% unique()
  
  # cluster <- "GC_00000001"
  
  for(cluster in gene_cluster_ids){
    df <- gene_clusters[gene_clusters$gene_cluster_id==cluster,]
    df$SSG_detailed <- "SSG"
    genome_names <- df$genome_name %>% as.character() %>% unique()
    
    if(nrow(df[df$genome_name==main_sample,])==1 && nrow(df)<=length(genome_names)){
      df2 <- df
      df2$SSG_shared_in_metagenomes <- ""
      if(nrow(df2)==1){
        df2$SSG_detailed <- paste0("Only in ", s)
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_metagenomes"] <- shared
      } else if(nrow(df2)==2){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 1 MAG")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_metagenomes"] <- shared
      } else if(nrow(df2)==3){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 2 MAGs")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_metagenomes"] <- shared
      } else if(nrow(df2)==4){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 3 MAGs")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_metagenomes"] <- shared
      } else if(nrow(df2)==5){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 4 MAGs")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_metagenomes"] <- shared
      }
      gene_cluster_filt <- gene_cluster_filt %>% rbind(df2)
    }
  }
  
  # Select gene cluster with at least i as unique gene within each gene cluster (reference)
  
  gene_clusters_ref <- read.table(paste0("Metapan/PAN/2_WOLBACHIA_PAN/Wolbachia-", s, "-ref-PAN-SUMMARY/Wolbachia_gene_clusters_summary.txt"), header=TRUE, fill=TRUE, sep="\t", quote = "")
  gene_clusters_ref_draft <- gene_clusters_ref
  
  gene_cluster_filt_ref <- data.frame()
  main_sample <- s
  gene_cluster_ids <- gene_clusters_ref$gene_cluster_id %>% as.character() %>% unique()
  
  for(cluster in gene_cluster_ids){
    df <- gene_clusters_ref[gene_clusters_ref$gene_cluster_id==cluster,]
    df$SSG_detailed <- "SSG"
    genome_names <- df$genome_name %>% as.character() %>% unique()
    if(nrow(df[df$genome_name==main_sample,])==1 && nrow(df)<=length(genome_names)){
      df2 <- df
      df2$SSG_shared_in_references <- ""
      if(nrow(df2)==1){
        df2$SSG_detailed <- paste0("Only in ", s)
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_references"] <- shared
      } else if(nrow(df2)==2){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 1 reference genome")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_references"] <- shared
      } else if(nrow(df2)==3){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 2 reference genomes")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_references"] <- shared
      } else if(nrow(df2)==4){
        df2$SSG_detailed <- paste0("Shared by ", s, " and 3 reference genomes")
        shared <- paste(df2[, "genome_name"] %>% as.character(), collapse=", ")
        df2[, "SSG_shared_in_references"] <- shared
      }
      gene_cluster_filt_ref <- gene_cluster_filt_ref %>% rbind(df2)
    }
  }
  
  gene_clusters <- gene_cluster_filt[gene_cluster_filt$genome_name==s,]
  gene_clusters <- gene_clusters %>% select(c(gene_callers_id, SSG_detailed, c(13:27)))
  colnames(gene_clusters)[1:2] <- c("corresponding_gene_call", "SSG")
  gene_clusters = gene_clusters[!duplicated(gene_clusters$corresponding_gene_call),]
  
  
  gene_clusters_ref <- gene_cluster_filt_ref[gene_cluster_filt_ref$genome_name==s,]
  gene_clusters_ref <- gene_clusters_ref %>% select(c(gene_callers_id, SSG_detailed, c(13:27)))
  colnames(gene_clusters_ref)[1:2] <- c("corresponding_gene_call", "SSG")
  gene_clusters_ref = gene_clusters_ref[!duplicated(gene_clusters_ref$corresponding_gene_call),]
  
  
  # Make nice tables for GC / SSG (draft)
  gene_clusters_draft <- gene_clusters_draft %>% merge(gene_cluster_filt %>% select(c(unique_id, SSG_detailed, SSG_shared_in_metagenomes)), 
                                                        by="unique_id", all.x=TRUE)
  colnames(gene_clusters_draft)[26:27] <- c(paste0("SSG_", s), paste0("SSG_", s, "_shared_in_metagenomes"))
  
  gene_clusters_ref_draft <- gene_clusters_ref_draft %>% merge(gene_cluster_filt_ref %>% select(c(unique_id, SSG_detailed, SSG_shared_in_references)), 
                                                        by="unique_id", all.x=TRUE)
  colnames(gene_clusters_ref_draft)[26] <- paste0("SSG_references")
  
  gene_clusters_draft2 <- gene_clusters_draft %>% select(c(gene_cluster_id, SCG, paste0("SSG_", s), paste0("SSG_", s, "_shared_in_metagenomes")))
  gene_clusters_draft2 <- unique(gene_clusters_draft2)
  
  gene_clusters_ref_draft2 <- gene_clusters_ref_draft %>% select(c(gene_cluster_id, SCG, SSG_references, SSG_shared_in_references))
  gene_clusters_ref_draft2 <- unique(gene_clusters_ref_draft2)
  
  # Merge SSG informations
  merged_gene_clusters <- gene_clusters %>% merge(gene_clusters_ref %>% select(c(corresponding_gene_call, SSG, SSG_shared_in_references)), 
                                                  by="corresponding_gene_call", all.x = TRUE, all.y=FALSE)
  
  merged_gene_clusters <- merged_gene_clusters %>% select(c(corresponding_gene_call, SSG.x, SSG.y, everything()))
  colnames(merged_gene_clusters)[2:3] <- c(paste0("SSG_", s), "SSG_ref")
  
  
  
  #### 5.1) Merge SCV and SSG tables
  
  scv_filtered <- read.table(paste0("SNV/snv-tables/SCVs/SCV_variability_", s, "_intra-MAGs_filtered_from_filtered_SNVs.txt"), header=TRUE, fill=TRUE)
  
  test <- unique(snv_to_merge %>% select(c(corresponding_gene_call, Bonneau, Bonneau_pct_alignment, MLST_WSP, MLST_WSP_pct_alignment)))
  scv_filtered <- scv_filtered %>% merge(test, by="corresponding_gene_call", all.x = TRUE)
  
  scv_scg <- scv_filtered %>% merge(merged_gene_clusters %>% select(c("corresponding_gene_call", 
                                                                      "SSG_shared_in_metagenomes", 
                                                                      "SSG_shared_in_references",
                                                                      "SSG_ref", 
                                                                      paste0("SSG_", s))), 
                                    by="corresponding_gene_call", all.x=TRUE)
  
  
  
  #### 5.2) Merge SAAV and SSG tables
  
  saav_filtered <- read.table(paste0("SNV/snv-tables/SAAVs/SAAV_variability_", s, "_intra-MAGs_filtered_from_filtered_SCVs.txt"), header=TRUE, fill=TRUE)
  
  test <- unique(snv_to_merge %>% select(c(corresponding_gene_call, Bonneau, Bonneau_pct_alignment, MLST_WSP, MLST_WSP_pct_alignment)))
  saav_filtered <- saav_filtered %>% merge(test, by="corresponding_gene_call", all.x = TRUE)
  
  saav_scg <- saav_filtered %>% merge(merged_gene_clusters %>% select(c("corresponding_gene_call", 
                                                                      "SSG_shared_in_metagenomes", 
                                                                      "SSG_shared_in_references",
                                                                      "SSG_ref", 
                                                                      paste0("SSG_", s))), 
                                    by="corresponding_gene_call", all.x=TRUE)
  
  
  
  #### 5.3) Merge SNV and SSG tables ####
  
  snv_scg <- snv_to_merge %>% merge(merged_gene_clusters %>% select(c("corresponding_gene_call", 
                                                                      "SSG_shared_in_metagenomes", 
                                                                      "SSG_shared_in_references",
                                                                      "SSG_ref", 
                                                                      paste0("SSG_", s))), 
                                    by="corresponding_gene_call", all.x=TRUE)
  
  snv_scg[is.na(snv_scg[,paste0("SSG_", s)]), paste0("SSG_", s)] <- "No SSG"
  snv_scg$SSG_ref[is.na(snv_scg$SSG_ref)] <- "No SSG"
  snv_scg$SSG_shared_in_metagenomes[is.na(snv_scg$SSG_shared_in_metagenomes)] <- "No SSG"
  snv_scg$SSG_shared_in_references[is.na(snv_scg$SSG_shared_in_references)] <- "No SSG"
  snv_scg$nb_samples_where_snv_shared[is.na(snv_scg$nb_samples_where_snv_shared)] <- "SNV not shared"
  snv_scg$snv_shared_in[is.na(snv_scg$snv_shared_in)] <- "SNV not shared"
  
  snv_scg$SNV_shared_in_metagenomes <- snv_scg$nb_samples_where_snv_shared
  snv_scg[snv_scg$SNV_shared_in_metagenomes!="SNV not shared", "SNV_shared_in_metagenomes"] <- "SNV shared"
  
  snv_scg_reduced <- snv_scg %>% select(c(entry_id,
                                          pos,
                                          pos_in_contig,
                                          split_name,
                                          contig_name,
                                          in_coding_gene_call,
                                          corresponding_gene_call,
                                          paste0("SSG_", s), 
                                          SSG_shared_in_metagenomes,
                                          SSG_ref,
                                          SSG_shared_in_references,
                                          nb_samples_where_snv_shared,
                                          snv_shared_in,
                                          SNV_shared_in_metagenomes,
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
                                           paste0("SSG_", s), 
                                           SSG_shared_in_metagenomes,
                                           SSG_ref,
                                           SSG_shared_in_references,
                                           nb_samples_where_snv_shared,
                                           snv_shared_in,
                                           SNV_shared_in_metagenomes,
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
  
  
  ## SNV dup
  snv_scg_dup_reduced <- snv_scg_reduced[snv_scg_reduced$nb_samples_where_snv_shared!="SNV not shared",]
  snv_scg_dup_complete <- snv_scg_complete[snv_scg_complete$nb_samples_where_snv_shared!="SNV not shared",]
  
  ## dup SNV only SSG (all)
  snv_only_scg_reduced <- snv_scg_dup_reduced[snv_scg_dup_reduced[, paste0("SSG_",s)]!="No SSG",]
  snv_only_scg_complete <- snv_scg_dup_complete[snv_scg_dup_complete[, paste0("SSG_",s)]!="No SSG",]
  
  ## dup SNV only SSG (samples)
  snv_only_scg_reduced_samples <- snv_only_scg_reduced[snv_only_scg_reduced$SSG_ref!="No SSG",]
  snv_only_scg_complete_samples <- snv_only_scg_complete[snv_only_scg_complete$SSG_ref!="No SSG",]
  
  ## dup SNV no scg
  snv_no_scg_reduced <- snv_scg_dup_reduced[snv_scg_dup_reduced$SSG_ref=="No SSG",]
  snv_no_scg_reduced <- snv_scg_dup_reduced[snv_scg_dup_reduced[, paste0("SSG_",s)]=="No SSG",]
  snv_no_scg_complete <- snv_scg_dup_complete[snv_scg_dup_complete$SSG_ref=="No SSG",]
  snv_no_scg_complete <- snv_scg_dup_complete[snv_scg_dup_complete[, paste0("SSG_",s)]=="No SSG",]
  
  ## snv_scg_layer to import in anvi-interactive
  snv_scg_layer <- snv_scg_complete %>% select(c(entry_id, paste0("SSG_", s), SSG_shared_in_metagenomes, 
                                                 SSG_ref, SSG_shared_in_references, 
                                                 nb_samples_where_snv_shared, snv_shared_in, SNV_shared_in_metagenomes))
  
  ## make an snv_scg_layer table with SSG unique in Wolbachia MAG for entropy figure
  snv_scg_layer_wt_unique_SSG <- snv_scg_layer
  snv_scg_layer_wt_unique_SSG[snv_scg_layer_wt_unique_SSG$SSG_ref==paste0("Only in ",s), c(paste0("SSG_", s), "SSG_shared_in_metagenomes", 
                                                                                              "SSG_ref", "SSG_shared_in_references", 
                                                                                              "nb_samples_where_snv_shared", "snv_shared_in", "SNV_shared_in_metagenomes")] <- "No SSG"
  snv_scg_layer_wt_unique_SSG[snv_scg_layer_wt_unique_SSG$SSG_shared_in_metagenomes==s, c(paste0("SSG_", s), "SSG_shared_in_metagenomes", 
                                                                                             "SSG_ref", "SSG_shared_in_references", 
                                                                                             "nb_samples_where_snv_shared", "snv_shared_in", "SNV_shared_in_metagenomes")] <- "No SSG"
  snv_scg_layer_wt_unique_SSG[snv_scg_layer_wt_unique_SSG$nb_samples_where_snv_shared=="No SSG", c("nb_samples_where_snv_shared", "snv_shared_in", "SNV_shared_in_metagenomes")] <- "SNV not shared"
  
  ## contingency table for SNV in coding gene call and on base positions
  
  ### keep SNV in coding gene call
  contingency <- table(in_coding_gene_call=snv_scg_complete$in_coding_gene_call, base_pos_in_codon=snv_scg_complete$base_pos_in_codon) %>% as.data.frame()
  contingency <- contingency[contingency$Freq!=0,]
  
  
  contingency2 <- table(SSG_metagenome=snv_scg_complete[, paste0("SSG_",s)]) %>% as.data.frame()
  str(contingency2)
  
  contingency3 <- table(SSG_ref=snv_scg_complete$SSG_ref) %>% as.data.frame()
  
  contingency4 <- table(SSG_metagenome=snv_scg_complete[, paste0("SSG_",s)], Bonneau=snv_scg_complete$Bonneau) %>% as.data.frame()
  
  contingency5 <- table(SSG_metagenome=snv_scg_complete[, paste0("SSG_",s)], MLST_WSP=snv_scg_complete$MLST_WSP) %>% as.data.frame()
  
  contingency6 <- table(base_pos_in_codon=snv_scg_complete$base_pos_in_codon) %>% as.data.frame()
  
  contingency7 <- table(SSG_metagenome=snv_scg_complete[, paste0("SSG_",s)], base_pos_in_codon=snv_scg_complete$base_pos_in_codon) %>% as.data.frame()
  
  contingency8 <- table(SSG_ref=snv_scg_complete$SSG_ref, base_pos_in_codon=snv_scg_complete$base_pos_in_codon) %>% as.data.frame()
  
  contingency9 <- table(Bonneau=snv_scg_complete$Bonneau, base_pos_in_codon=snv_scg_complete$base_pos_in_codon) %>% as.data.frame()
  
  
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
                    paste0("SSG_", s),
                    "SSG_shared_in_metagenomes",
                    "SSG_ref",
                    "SSG_shared_in_references", 
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
  colnames(snv_scg_for_draft)[34] <- "SNV_shared_in_metagenomes"
  colnames(snv_scg_for_draft)[36] <- "cid_genes_from_Bonneau_et_al._2018"
  colnames(snv_scg_for_draft)[38] <- "MLST_and_wsp_genes_from_PubMLST"
  
  #### 7) Make nice tables to export (for Excel) ####
  
  wb2 <- createWorkbook()
  addWorksheet(wb2, "All filtered SNVs")
  writeData(wb2, "All filtered SNVs", snv_scg_complete, startRow = 1, startCol = 1)
  addWorksheet(wb2, "Only duplicated SNVs")
  writeData(wb2, "Only duplicated SNVs", snv_scg_dup_complete, startRow = 1, startCol = 1)
  addWorksheet(wb2, "Dup SNVs in SSGs")
  writeData(wb2, "Dup SNVs in SSGs", snv_only_scg_complete, startRow = 1, startCol = 1)
  addWorksheet(wb2, "Dup SNVs in sample SSGs")
  writeData(wb2, "Dup SNVs in sample SSGs", snv_only_scg_complete_samples, startRow = 1, startCol = 1)
  addWorksheet(wb2, "Dup SNVs not in SSGs")
  writeData(wb2, "Dup SNVs not in SSGs", snv_no_scg_complete, startRow = 1, startCol = 1)
  addWorksheet(wb2, "contingency SNVs in coding")
  writeData(wb2, "contingency SNVs in coding", contingency, startRow = 1, startCol = 1)
  writeData(wb2, "contingency SNVs in coding", contingency2, startRow = 1, startCol = 5)
  writeData(wb2, "contingency SNVs in coding", contingency3, startRow = 1, startCol = 10)
  writeData(wb2, "contingency SNVs in coding", contingency4, startRow = 1, startCol = 15)
  writeData(wb2, "contingency SNVs in coding", contingency5, startRow = 1, startCol = 20)
  writeData(wb2, "contingency SNVs in coding", contingency6, startRow = 1, startCol = 25)
  writeData(wb2, "contingency SNVs in coding", contingency7, startRow = 1, startCol = 30)
  writeData(wb2, "contingency SNVs in coding", contingency8, startRow = 1, startCol = 35)
  writeData(wb2, "contingency SNVs in coding", contingency9, startRow = 1, startCol = 40)
  saveWorkbook(wb2, file = paste0("SNV/snv-tables/intra-MAGs/SNV_SSG_table_", s, "-complete.xlsx"), overwrite = TRUE)
  
  
  wb3 <- createWorkbook()
  addWorksheet(wb3, "All filtered SNVs")
  writeData(wb3, "All filtered SNVs", snv_scg_reduced, startRow = 1, startCol = 1)
  addWorksheet(wb3, "Only duplicated SNVs")
  writeData(wb3, "Only duplicated SNVs", snv_scg_dup_reduced, startRow = 1, startCol = 1)
  addWorksheet(wb3, "Dup SNVs in SSGs")
  writeData(wb3, "Dup SNVs in SSGs", snv_only_scg_reduced, startRow = 1, startCol = 1)
  addWorksheet(wb3, "Dup SNVs in sample SSGs")
  writeData(wb3, "Dup SNVs in sample SSGs", snv_only_scg_reduced_samples, startRow = 1, startCol = 1)
  addWorksheet(wb3, "Dup SNVs not in SSGs")
  writeData(wb3, "Dup SNVs not in SSGs", snv_no_scg_reduced, startRow = 1, startCol = 1)
  saveWorkbook(wb3, file = paste0("SNV/snv-tables/intra-MAGs/SNV_SSG_table_", s, "-reduced.xlsx"), overwrite = TRUE)
  
  wb4 <- createWorkbook()
  addWorksheet(wb4, "Filtered SNVs, SSG, functions")
  writeData(wb4, "Filtered SNVs, SSG, functions", snv_scg_for_draft, startRow = 1, startCol = 1)
  saveWorkbook(wb4, file = paste0("SNV/snv-tables/intra-MAGs/SNV_SSG_table_", s, "-for-draft.xlsx"), overwrite = TRUE)
  
  write.table(snv_scg_layer, paste0("SNV/snv-tables/intra-MAGs/variability-",s, "-layers-SSG.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_scg_layer_wt_unique_SSG, paste0("SNV/snv-tables/intra-MAGs/variability-",s, "-layers-SSG-wt-unique.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
  snv_dup2 <- snv_dup %>% select(-c(entry_id))
  snv_dup2 <- snv_dup2 %>% select(c(entry_id_ref, unique_pos_identifier, pos, pos_in_contig, contig_name, split_name, sample_id, corresponding_gene_call, everything()))
  snv_dup2 <- snv_dup2[order(snv_dup2$entry_id_ref),]
  colnames(snv_dup2)[1] <- "entry_id"
  write.table(snv_dup2, paste0("SNV/snv-tables/intra-MAGs/SNV-",s, "-filtered-shared.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
  
  write.table(gene_clusters_draft2, paste0("SNV/snv-tables/intra-MAGs/GC-SSG-",s, "-metagenomes-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(gene_clusters_ref_draft2, paste0("SNV/snv-tables/intra-MAGs/GC-SSG-",s, "-references-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
  write.table(scv_scg, paste0("SNV/snv-tables/intra-MAGs/SCV-",s , "-filtered-SCG-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(saav_scg, paste0("SNV/snv-tables/intra-MAGs/SAAV-",s , "-filtered-SCG-table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
}
