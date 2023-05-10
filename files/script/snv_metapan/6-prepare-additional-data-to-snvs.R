require(tidyverse)
samples <- c("O03", "O07", "O11", "O12", "M11")

# O03, O12 -> no hit
# i <- "O11"

# Make table of cid IDs from Bonneau et al. 2018: https://www.nature.com/articles/s41467-017-02749-w
cid_variants <- c("cidA_I(alpha/1)", "cidA_I(gamma/1)", "cidA_I(gamma/2)", "cidA_I(beta/2)", "cidA_II(alpha/1)", "cidA_II(alpha/2)",
                  "cidA_II(beta/2)", "cidA_III(alpha/1)", "cidA_III(beta/2)", "cidA_III(beta/1)", "cidA_III(beta/3)", "cidA_IV(alpha/1)", 
                  "cidA_IV(alpha/2)","cidA_IV(gamma/1)", "cidA_IV(gamma/2)", "cidA_IV(delta/1)", "cidA_IV(delta/2)", "cidA_IV(beta/1)", 
                  "cidA_IV(beta/2)",
                  "cidB_I(a/1)", "cidB_I(a/2)", "cidB_I(b/1)", "cidB_I(b/2)", "cidB_II(a/2)", "cidB_II(a/1)", "cidB_III(c/1)", "cidB_III(a/1)", 
                  "cidB_III(b/1)", "cidB_IV(a/1)", "cidB_IV(a/2)", "cidB_IV(b/1)", "cidB_IV(b/2)", "cidB_IV(a/3)", "cidB_IV(b/3)")
acc_numbers <- c("MF444963", "MF444964", "MF444965", "MF444966", "MF444967", "MF444968", "MF444969", "MF444970", "MF444971", "MF444972", "MF444973",
                 "MF444974", "MF444975", "MF444976", "MF444977", "MF444978", "MF444979", "MF444980", "MF444981",
                 "MF444982", "MF444983", "MF444984", "MF444985", "MF444986", "MF444987", "MF444988", "MF444989", "MF444990", "MF444991", "MF444992",
                 "MF444993", "MF444994", "MF444995", "MF444996")

cid_var <- c("cidA", "cidB")
acc_num <- c("KF114896", "KX077602")


# Make table of cid IDs from Beckmann et al. 2013: https://pubmed.ncbi.nlm.nih.gov/23856508/
cids <- data.frame(cid_variants, Bonneau=acc_numbers)
cids$Bonneau <- cids$Bonneau %>% paste0(".1")

cids_Beckmann <- data.frame(cid_var, Beckmann=acc_num)
cids_Beckmann$Beckmann <- cids_Beckmann$Beckmann %>% paste0(".1")


# Add cid and MSLT annotations to SNV tables (intra-samples level)

for(i in samples){
  
  ## import filtered snv table
  snv_dep <- read.table(paste0("/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-samples/variability-", i, "-filtered.txt"), header=TRUE)
  
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
  
  snv_layers_full <- snv_layers_full %>% merge(cids, by="Bonneau", all.x=TRUE)
  snv_layers_full <- snv_layers_full %>% select(-c(Bonneau))
  colnames(snv_layers_full)[ncol(snv_layers_full)] <- "Bonneau"
  
  snv_layers_full$Bonneau <- droplevels(snv_layers_full$Bonneau)
  levels(snv_layers_full$Bonneau)
  
  
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
  
  snv_layers_full <- snv_layers_full %>% merge(cids_Beckmann, by="Beckmann", all.x=TRUE)
  snv_layers_full <- snv_layers_full %>% select(-c(Beckmann))
  colnames(snv_layers_full)[ncol(snv_layers_full)] <- "Beckmann"
  
  snv_layers_full$Beckmann <- droplevels(snv_layers_full$Beckmann)
  levels(snv_layers_full$Beckmann)
  
  
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
  write.table(snv_layers_full, paste0("/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-samples/variability-", i, "-layers-full.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snv_layers, paste0("/Volumes/Elements/Hiseq/metapangenomics/output/SNV/snv-tables/intra-samples/variability-", i, "-layers.txt"), sep="\t", quote=FALSE, row.names=FALSE)
}