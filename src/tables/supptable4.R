library(httr)
library(jsonlite)
suppressPackageStartupMessages(library(UniProt.ws))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(biomaRt))

obtain.protein.desc <- function(protein.names) {
  # Same protein names have ; select everything before ;
  protein.split <- strsplit(protein.names, split=";")
  for (i in 1:length(protein.split)) {
    protein.split[[i]] <- protein.split[[i]][1]
  }
  protein.names <- protein.split %>% unlist
  
  # ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  # protein_table <- getBM(attributes = c("external_gene_name","uniprotswissprot", "description","hgnc_symbol","gene_biotype"),#,"entrezgene" 
  #                        filters = "external_gene_name", 
  #                        values = protein.names, 
  #                        mart = ensembl)
  # # print(protein_table)
  # # Only obtain where protein_names$uniprotswissprot is not na
  # protein_table <- protein_table[which(protein_table$uniprotswissprot!=""),]
  # print(protein_table)
  # stop()
  
  
  # Convert protein names to uniprot ID's
  id_table <- bitr(protein.names,
                   fromType="SYMBOL",
                   toType="UNIPROT",
                   OrgDb = org.Hs.eg.db)
  # Obtain first mapping
  id_table <- id_table %>%
    group_by(SYMBOL) %>%
    dplyr::slice(1) %>%
    ungroup()
  # Obtain proteins that have no match
  no.match.protein <- setdiff(protein.names, id_table$SYMBOL)
  # stop()
  up <- UniProt.ws(taxId=9606)
  result <- select(
    x = up,
    keys = id_table$UNIPROT,
    columns = c("id", "protein_name", "gene_primary", "organism_id"),
    keytype = "UniProtKB"
  )
  colnames(result) <- c("from_id", "uniprot_id", "protein_name", "protein_description", "primary_gene_id", "organism_id")
  # print()
  result <- result[!duplicated(result$from_id),]
  result <- result %>% 
    dplyr::select(all_of(c("uniprot_id", "protein_name", "protein_description")))
  result$primary_gene_id <- id_table$SYMBOL

  # Append the no matches on
  result <- rbind(result, data.frame("uniprot_id"=NA, "protein_name"=NA,
                                     "protein_description"=NA, "primary_gene_id"=no.match.protein))
  result
}
save.dir <- file.path("src", "tables", "supp_table_4")
R.utils::mkdirs(save.dir)
# Load GBMVMET
gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
gbm.v.met.protein.info <- obtain.protein.desc(rownames(gbm.v.met.processed))
# print(gbm.v.met.protein.info)
# Cbind gbm.v.met.protein.info based on matches in the primary_gene_id
match.id <- match(gbm.v.met.protein.info$primary_gene_id,
                  rownames(gbm.v.met.processed))
gene.names <- rownames(gbm.v.met.processed)
new.gbm.v.met.data <- cbind(gbm.v.met.protein.info[match.id,],
                            gbm.v.met.processed)
new.gbm.v.met.data$primary_gene_id <- gene.names

# new.gbm.v.met.data <- cbind(gbm.v.met.protein.info, gbm.v.met.processed)

write.xlsx(new.gbm.v.met.data, file.path(save.dir, "gbm_v_met_supp_4.xlsx"))


