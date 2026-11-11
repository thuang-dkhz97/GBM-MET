suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library(EnhancedVolcano))


# base_dir <- "/home/abhivij/UNSW/VafaeeLab/GBMPlasmaEV"
# setwd(base_dir)

# source("scripts/R/utils.R")
# source("scripts/R/plot_data.R")



gbm.v.met.RNA <- openxlsx::read.xlsx("data/GBM-BMET.xlsx", rowNames=T)
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)
match.i <- match(colnames(gbm.v.met.RNA), gbm.v.met.labels$Sample)
gbm.v.met.labels <- gbm.v.met.labels[match.i,]
colours <- c("#34a1e6", "#27d83a", "#b5b5b5")
# create_volcano_plot(gbm.v.met.RNA, gbm.v.met.labels, c("MET", "GBM"),
#                     "GBM compared to MET", "GBMvsMET-RNA-volcano.png",
#                     save.dir,
#                     logFC_cutoff=log2(1.2),
#                     use_p_val=T, colours=colours, k=10, RNA=T)
gbm.v.met.labels$Label <- factor(gbm.v.met.labels$Label)
dim(gbm.v.met.RNA)
dim(gbm.v.met.labels)

all(gbm.v.met.RNA == round(gbm.v.met.RNA), na.rm = TRUE)
dds <- DESeqDataSetFromMatrix(countData=gbm.v.met.RNA,
                              colData = gbm.v.met.labels,
                              design=~Label)


# umi_counts <- read.csv("Data/RNA/umi_counts.csv", row.names = 1)
# metadata <- read.csv("Data/metadata.csv") %>%
#   column_to_rownames("SUBJECT_ORIGINAL")
# umi_counts <- umi_counts[, rownames(metadata)]
# 
# all(rownames(metadata) == colnames(umi_counts))
# 
# dds <- DESeqDataSetFromMatrix(countData = umi_counts,
#                               colData = metadata,
#                               design = ~ GROUP_Q1to6)
# dds
# 
# print(dim(counts(dds)))
# 
# # keep <- rowSums(counts(dds)) >= 30
# # dds <- dds[keep,]
# 
# print(dim(counts(dds)))
# 
# dds <- DESeq(dds)
# 
# contrast_vec <- c("MET", "HC")
# comparison_num <- 2
# get_comparison_results <- function(contrast_vec, comparison_num, umi_counts, metadata){
#   res <- results(dds, contrast = c("GROUP_Q1to6", contrast_vec[1], contrast_vec[2]), alpha = 0.05)
#   res
#   print(summary(res))
#   
#   result <- data.frame(res) %>%
#     filter(!is.na(padj)) %>%
#     rownames_to_column("rna") %>%
#     select(rna, log2FoldChange, pvalue, padj) %>%
#     dplyr::rename(Molecule = rna, adjPVal = padj, pVal = pvalue, logFC = log2FoldChange)
#   
#   title <- paste(contrast_vec, collapse = " Vs ")
#   
#   file_name <- paste("de", "comp", comparison_num, 
#                      paste(contrast_vec, collapse = "_"), sep = "_")
#   file_name <- paste(file_name, "png", sep = ".")
#   
#   
#   plot_transcriptomic_data(umi_counts[(result %>% filter(adjPVal < 0.05))$Molecule, ], 
#                            metadata, contrast_vec, "Adj pval < 0.05 transcripts",
#                            file_name = "adjpval")
#   plot_transcriptomic_data(umi_counts[(result %>% filter(adjPVal < 0.05))$Molecule, ], 
#                            metadata, contrast_vec, "Adj pval < 0.05 transcripts", dim_red = "umap",
#                            file_name = "adjpval")
#   
#   plot_transcriptomic_data(umi_counts[(result %>% filter(pVal < 0.05))$Molecule, ], 
#                            metadata, contrast_vec, "pval < 0.05 transcripts",
#                            file_name = "pval")
#   plot_transcriptomic_data(umi_counts[(result %>% filter(pVal < 0.05))$Molecule, ], 
#                            metadata, contrast_vec, "pval < 0.05 transcripts", dim_red = "umap",
#                            file_name = "pval")
#   
#   
#   create_volcano_plot(result, 
#                       title = title, 
#                       file_name = paste("volcano", file_name, sep = "_"), 
#                       dir_path = "plots/de/transcriptomic/deseq2", logFC_cutoff = 1)
#   create_volcano_plot(result, 
#                       title = title, 
#                       file_name = paste("volcano_pval", file_name, sep = "_"), 
#                       dir_path = "plots/de/transcriptomic/deseq2", logFC_cutoff = 1,
#                       use_p_val = TRUE)  
#   
# }
# 
# 
# 
# get_comparison_results(c("PREOPE", "MET"), 2, umi_counts, metadata)
# get_comparison_results(c("PREOPE", "HC"), 2, umi_counts, metadata)
# get_comparison_results(c("MET", "HC"), 2, umi_counts, metadata)
# 
# 
# 
# 
# plot_transcriptomic_data <- function(data, metadata, categories, title, dim_red = "pca",
#                                      file_name = ""){
#   metadata <- metadata %>%
#     filter(GROUP_Q1to6 %in% categories)
#   print("metadata size")
#   print(dim(metadata))
#   
#   data <- data[, rownames(metadata)]
#   print(all(rownames(metadata) == colnames(data)))
#   groups <- metadata$GROUP_Q1to6
#   
#   file_name <- paste("umi_counts", file_name, paste(categories, collapse = "_"), sep = "_")
#   file_name <- paste(file_name, "png", sep = ".")
#   
#   plot_data(t(data), file_name, title, groups = groups, dim_red = dim_red,
#             colour_label = "Labels")
# }
# 
# plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "MET"), "All transcripts")
# plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "MET"), "All transcripts", dim_red = "umap")
# 
# plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "HC"), "All transcripts")
# plot_transcriptomic_data(umi_counts, metadata, c("PREOPE", "HC"), "All transcripts", dim_red = "umap")
# 
# plot_transcriptomic_data(umi_counts, metadata, c("MET", "HC"), "All transcripts")
# plot_transcriptomic_data(umi_counts, metadata, c("MET", "HC"), "All transcripts", dim_red = "umap")
# 
# 
# 
# 
# 
# ###############
# 
# #DE analysis for all transcriptomics data - PREOPE, MET, HC - 2024 Jan 
# 
# 
# 
# umi_counts <- read.csv("Data/RNA_all/newquant_Nov2023_umi_counts_PREOPE_MET_HC_filter90.csv", row.names = 1)
# metadata <- read.table("Data/transcriptomic_phenotype_PREOPE_MET_HC.txt", header=TRUE, sep="\t") %>%
#   filter(!is.na(data_cohort)) %>%
#   mutate(condition = case_when(Condition == "Pre-op" ~ "PREOPE",
#                                Condition == "Metastatic" ~ "MET",
#                                Condition == "Healthy Control" ~ "HC"))
# umi_counts <- umi_counts[, metadata$Sample]
# 
# dds <- DESeqDataSetFromMatrix(countData = umi_counts,
#                               colData = metadata,
#                               design = ~ condition + data_cohort)
# 
# print(dim(counts(dds)))
# 
# dds <- DESeq(dds)
# resultsNames(dds)
# 
# plot_dir_path <- "plots_RNA_all/PREOPE_MET_HC/volcano/" 
# if(!dir.exists(plot_dir_path)){
#   dir.create(plot_dir_path, recursive = TRUE)  
# }
# 
# res <- results(dds, contrast = c("condition", "PREOPE", "MET"), alpha = 0.05)
# 
# res
# print(summary(res))
# 
# result <- as.data.frame(res) %>%
#   filter(!is.na(padj)) %>%
#   rownames_to_column("rna") %>%
#   select(rna, log2FoldChange, pvalue, padj) %>%
#   dplyr::rename(Molecule = rna, adjPVal = padj, pVal = pvalue, logFC = log2FoldChange) %>%
#   arrange(logFC)
# 
# 
# 
# 
# logFC_cutoff = 0.585
# #approx 1.5 FC
# k = 10
# title = "PREOPE Vs MET"
# plot_file_name = "volcano_1_PREOPEVsMET.png"
# dir_path = plot_dir_path
# 
# diff_utils(result, title = "PREOPE Vs MET",
#            plot_file_name = "volcano_1_PREOPEVsMET.png",
#            dir_path = plot_dir_path,
#            logFC_cutoff = 0.585)
# 
# 
# 
# res <- results(dds, contrast = c("condition", "PREOPE", "HC"), alpha = 0.05)
# res
# print(summary(res))
# result <- as.data.frame(res) %>%
#   filter(!is.na(padj)) %>%
#   rownames_to_column("rna") %>%
#   select(rna, log2FoldChange, pvalue, padj) %>%
#   dplyr::rename(Molecule = rna, adjPVal = padj, pVal = pvalue, logFC = log2FoldChange) %>%
#   arrange(logFC)
# diff_utils(result, title = "PREOPE Vs HC",
#            plot_file_name = "volcano_2_PREOPEVsHC.png",
#            dir_path = plot_dir_path,
#            logFC_cutoff = 0.585)
# 
# 
# res <- results(dds, contrast = c("condition", "MET", "HC"), alpha = 0.05)
# res
# print(summary(res))
# result <- as.data.frame(res) %>%
#   filter(!is.na(padj)) %>%
#   rownames_to_column("rna") %>%
#   select(rna, log2FoldChange, pvalue, padj) %>%
#   dplyr::rename(Molecule = rna, adjPVal = padj, pVal = pvalue, logFC = log2FoldChange) %>%
#   arrange(logFC)
# diff_utils(result, title = "MET Vs HC",
#            plot_file_name = "volcano_3_METVsHC.png",
#            dir_path = plot_dir_path,
#            logFC_cutoff = 0.263)  #log2(1.2)
# 
# 
# #NOTE - do not use the below function in future - creating a common function for prot and tra in utils_diff.R
# diff_utils <- function(result, title, plot_file_name, dir_path,
#                        k = 10,
#                        logFC_cutoff) {
#   
#   
#   #arrange by significance needs to be done so as to show 
#   # the labels in the order - downregulated, not significant, upregulated in the volcano plot 
#   result <- result %>%
#     mutate(significance = case_when(logFC <= -logFC_cutoff & adjPVal <= 0.05 ~ 'Downregulated',
#                                     logFC >= logFC_cutoff & adjPVal <= 0.05 ~ 'Upregulated',
#                                     TRUE ~ 'Not significant')) %>% 
#     mutate(colour = case_when(significance == 'Downregulated' ~ '#E8495C',
#                               significance == 'Upregulated' ~ '#38ACE2',
#                               TRUE ~ 'grey')) %>%
#     arrange(significance)
#   
#   upreg <- result %>%
#     filter(significance == 'Upregulated') %>%
#     arrange(desc(logFC))
#   downreg <- result %>%
#     filter(significance == "Downregulated") %>%
#     arrange(logFC)
#   
#   if(nrow(upreg) >= k){
#     top_molecules <- upreg[1:k, "Molecule"]
#   } else{
#     top_molecules <- upreg[1:nrow(upreg), "Molecule"]
#   }
#   if(nrow(downreg) >= k){
#     top_molecules <- c(top_molecules, downreg[1:k, "Molecule"])
#   } else{
#     top_molecules <- c(top_molecules, downreg[1:nrow(downreg), "Molecule"])
#   }
#   
#   result <- result %>%
#     mutate(label = case_when(Molecule %in% top_molecules ~ Molecule,
#                              TRUE ~ NA_character_))
#   
#   # result <- result %>%
#   #   filter(!is.na(label)) %>%
#   #   arrange(desc(logFC))
#   
#   sig <- rbind(upreg, downreg) %>%
#     select(-c(significance, colour)) %>%
#     arrange(desc(logFC))
#   # print(sig)
#   # print(paste("Upregulated:", nrow(sig[sig$logFC > 0,]), sep=" "))
#   # print(paste("Downregulated:", nrow(sig[sig$logFC < 0,]), sep=" "))
#   
#   keyvals <- result$colour
#   names(keyvals) <- result$significance
#   
#   upreg_count <- sum(names(keyvals) == 'Upregulated')
#   downreg_count <- sum(names(keyvals) == 'Downregulated')
#   ns_count <- sum(names(keyvals) == 'Not significant')
#   
#   names(keyvals) <- gsub('Upregulated', paste0('Upregulated (', upreg_count, ')'), names(keyvals))
#   names(keyvals) <- gsub('Downregulated', paste0('Downregulated (', downreg_count, ')'), names(keyvals))
#   names(keyvals) <- gsub('Not significant', paste0('Not significant (', ns_count, ')'), names(keyvals))
#   
#   volcanoplot <- EnhancedVolcano(result,
#                                  lab = result$label,
#                                  x='logFC',
#                                  y='adjPVal',
#                                  ylab = bquote(~-Log[10] ~ italic(adjP)),
#                                  pCutoff = 0.05,
#                                  FCcutoff = logFC_cutoff,
#                                  colCustom = keyvals,
#                                  colAlpha = 1,
#                                  title = title,
#                                  legendPosition = 'bottom',
#                                  subtitle = '',
#                                  labSize = 3,
#                                  drawConnectors = T,
#                                  max.overlaps = Inf,
#                                  arrowheads = F)
#   plot(volcanoplot)
#   ggsave(paste0(dir_path, plot_file_name), units = "cm", width = 25)
#   
#   data_file_name <- sub(".png", ".csv", plot_file_name, fixed = TRUE)
#   write.csv(result, paste0(dir_path, data_file_name), row.names = FALSE)
# }
