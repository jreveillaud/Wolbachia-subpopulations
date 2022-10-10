require(tidyverse)
samples <- c("O03", "O07", "O11", "O12", "M11")

# O03, O12 -> no hit
#i <- "O11"

for(i in samples){
  
  ## import filtered snv table
  snv_dep <- read.table(paste0("/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-MAGs/variability-", i, "-filtered-2.txt"), header=TRUE)
  
  ## add cid from Bonneau et al. 2018 assignments
  cid_bonneau <- read.table(paste0("/Volumes/Elements/Hiseq/metapangenomics/output/genes/cid/cidA_cidB_Bonneau_in_", i, "_best_hits.txt"), header=TRUE)
  cid_bonneau[,"pct_alignment"] <- gsub(",", ".", cid_bonneau [,"pct_alignment"])
  cid_bonneau$pct_alignment <- cid_bonneau$pct_alignment %>% as.numeric()
  
  cid_bonneau <- cid_bonneau %>% select(c(gene_callers_id, cidA.cidB_assignment, pct_alignment))
  colnames(cid_bonneau) <- c("corresponding_gene_call", "Bonneau", "Bonneau_pct_alignment")
  
  snv_layers_full <- snv_dep %>% merge(cid_bonneau, by="corresponding_gene_call", all.x=TRUE)
  if(nrow(cid_bonneau)!=0){
    snv_layers_full$Bonneau <- droplevels(snv_layers_full$Bonneau)
  }
  
  if(i=="O11" | i=="M11"){
    levels(snv_layers_full$Bonneau) <- "cidA_II(alpha/1)"
  } 
  
  if(i=="O07"){
    levels(snv_layers_full$Bonneau) <- c("cidB_I(a/1)", "cidB_IV(b/1)")
  }
  
  
  ## add cid from Beckmann et al. 2013 assignments
  cid_beckmann <- read.table(paste0("/Volumes/Elements/Hiseq/metapangenomics/output/genes/cid/cidA_cidB_Beckmann_in_", i, "_best_hits.txt"), header=TRUE)
  cid_beckmann[,"pct_alignment"] <- gsub(",", ".", cid_beckmann [,"pct_alignment"])
  cid_beckmann$pct_alignment <- cid_beckmann$pct_alignment %>% as.numeric()
  
  cid_beckmann <- cid_beckmann %>% select(c(gene_callers_id, cidA.cidB_assignment, pct_alignment))
  colnames(cid_beckmann) <- c("corresponding_gene_call", "Beckmann", "Beckmann_pct_alignment")
  
  snv_layers_full <- snv_layers_full %>% merge(cid_beckmann, by="corresponding_gene_call", all.x=TRUE)
  
  if(nrow(cid_beckmann)!=0){
    snv_layers_full$Beckmann <- droplevels(snv_layers_full$Beckmann)
  }
  
  if(i=="O07"){
    levels(snv_layers_full$Beckmann) <- "cidB"
  }
  
  if(i=="O11" | i=="M11"){
    levels(snv_layers_full$Beckmann) <- "cidA"
  }
  
  ## add MLST and wsp from PubMLST assignments
  mlst <- read.table(paste0("/Volumes/Elements/Hiseq/metapangenomics/output/genes/mlst/mlst_wsp_in_", i, "_best_hits.txt"), header=TRUE)
  mlst[,"pct_alignment"] <- gsub(",", ".", mlst [,"pct_alignment"])
  mlst$pct_alignment <- mlst$pct_alignment %>% as.numeric()
  
  mlst <- mlst %>% select(c(gene_callers_id, mlst.wsp_assignment, pct_alignment))
  colnames(mlst) <- c("corresponding_gene_call", "MLST_WSP", "MLST_WSP_pct_alignment")
  
  snv_layers_full <- snv_layers_full %>% merge(mlst, by="corresponding_gene_call", all.x=TRUE)
  snv_layers_full$MLST_WSP <- droplevels(snv_layers_full$MLST_WSP)
  
  ## prepare layers to add in profile.db of interactive
  snv_layers <- snv_layers_full %>% select(c(entry_id, 
                                  split_name,
                                  unique_pos_identifier, 
                                  pos, 
                                  pos_in_contig,
                                  in_coding_gene_call, 
                                  coverage, 
                                  entropy, 
                                  Bonneau, 
                                  Bonneau_pct_alignment, 
                                  Beckmann, 
                                  Beckmann_pct_alignment, 
                                  MLST_WSP,
                                  MLST_WSP_pct_alignment))
  snv_layers$entry_id <- paste0("p_", sprintf("%08d", as.numeric(snv_layers$entry_id)))
  
  snv_layers_full$entry_id <- paste0("p_", sprintf("%08d", as.numeric(snv_layers_full$entry_id)))
  
  ## write snv layer tables
  write.table(snv_layers_full, paste0("/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-MAGs/variability-", i, "-layers-full.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_layers, paste0("/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-MAGs/variability-", i, "-layers.txt"), sep="\t", quote=FALSE, row.names=FALSE)
}