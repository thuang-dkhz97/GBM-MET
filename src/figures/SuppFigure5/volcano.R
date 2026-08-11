suppressPackageStartupMessages(library(ggrepel))
library(limma)

create_volcano_plot <- function(data, phenotype, conditions, title, file_name, dir_path = "",
                                p_val_cutoff = 0.05, logFC_cutoff = 5,
                                use_p_val = FALSE,
                                k = 5,
                                colours) {
  label <- phenotype$Label[match(colnames(data), phenotype$Sample)]
  model_matrix <- model.matrix(~ 0 + label)
  colnames(model_matrix) <- sub("label", "", colnames(model_matrix), fixed = TRUE)
  contr_matrix <- makeContrasts(contrasts = paste0(conditions[2], " - ", conditions[1]),
                levels = colnames(model_matrix))
  fit <- lmFit(data, model_matrix)
  fit <- contrasts.fit(fit, contr_matrix)
  efit <- eBayes(fit)
  top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
    rownames_to_column("protein")
  result <- top.table %>%
    dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
    dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
    arrange(logFC)
  if(use_p_val){
    p_val_column <- "PVal"
  } else{
    p_val_column <- "adjPVal"
  }
  result$diffexpr <- "NO"
  result$diffexpr[result$logFC > logFC_cutoff & result[, p_val_column] < p_val_cutoff] <- "UP"
  result$diffexpr[result$logFC < -logFC_cutoff & result[, p_val_column] < p_val_cutoff] <- "DOWN"
  
  rownames(result) <- result$Molecule
  up_reg_result <- result %>%
    filter(logFC > 0 & .[[p_val_column]] < p_val_cutoff)
  if(dim(up_reg_result)[1] == 0){
    up_reg_result <- result %>%
      filter(logFC > 0)    
  }
  up_reg_indices <- order(up_reg_result$logFC,decreasing = TRUE)[1:k]
  up_reg_indices <- rownames(up_reg_result[up_reg_indices,])
  
  down_reg_result <- result %>%
    filter(logFC < 0 & .[[p_val_column]] < p_val_cutoff)
  if(dim(down_reg_result)[1] == 0){
    down_reg_result <- result %>%
      filter(logFC < 0) 
  }  
  down_reg_indices <- order(-down_reg_result$logFC, decreasing = TRUE)[1:k] 
  down_reg_indices <- rownames(down_reg_result[down_reg_indices,])
  
  result$de_mol <- NA
  result[up_reg_indices, "de_mol"] <- result[up_reg_indices, "Molecule"]
  result[down_reg_indices, "de_mol"] <- result[down_reg_indices, "Molecule"]
  
  
  # colours <- c("red", "blue", "grey")
  names(colours) <- c("UP", "DOWN", "NO")
  
  vol_plot <- ggplot(data=result, aes(x=logFC, y=-log10(!!as.symbol(p_val_column)), col=diffexpr, label=de_mol)) +
    geom_point(size=3) + 
    geom_text_repel(max.overlaps = 20, colour = "black") +
    geom_vline(xintercept=c(-logFC_cutoff, logFC_cutoff), col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(p_val_cutoff), col="black", linetype = "dashed") +
    scale_colour_manual(values = colours) +
    labs(title = title, colour = "Diff Expr",
         y=expression(~{"-"}~Log[10]~ P),
         x=expression(Log[2]~ fold ~change)) +
    theme_bw() +
    theme(legend.position="none",
          axis.text=element_text(size=25),
          axis.title=element_text(size=30),
          title=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))
  
  ggsave(file.path(dir_path, file_name), vol_plot, width=30, height=20, units="cm") 

  # save diff exp
  comp.name <- paste(rev(conditions), collapse="-vs-")
  write.xlsx(result[,1:(ncol(result)-1)], file=file.path(dir_path, paste0(comp.name, "-diff.xlsx")))
}

save.dir <- file.path("src", "figures", "SuppFigure5")

# ================== GBM V MET

gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)
colours <- c("#34a1e6", "#27d83a", "#b5b5b5")
create_volcano_plot(gbm.v.met.processed, gbm.v.met.labels, c("MET", "GBM"),
                    "GBM compared to MET", "GBMvsMET-volcano.png",
                    save.dir,
                    logFC_cutoff=log2(1.5),
                    use_p_val=T, colours=colours, k=10)

# =============== GBM V HC

gbm.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-hc-norm.csv",
                                row.names=1)
gbm.v.hc.labels <- read.table("src/process-2.6/combine-data-3/GBMVsHC_output_labels.txt", header=T)
colours <- c("#34a1e6", "#f21210", "#b5b5b5")
create_volcano_plot(gbm.v.hc.processed, gbm.v.hc.labels, c("HC", "GBM"),
                    "GBM compared to HC", "GBMvsHC-volcano.png",
                    save.dir,
                    logFC_cutoff=log2(1.5),
                    use_p_val=T,
                    colours=colours,
                    k=10)

# ============= MET V HC

met.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/met-v-hc-norm.csv",
                               row.names=1)
met.v.hc.labels <- read.table("src/process-2.6/combine-data-3/METVsHC_output_labels.txt", header=T)
colours <- c("#27d83a", "#f21210", "#b5b5b5")
create_volcano_plot(met.v.hc.processed, met.v.hc.labels, c("HC", "MET"),
                    "MET compared to HC", "METvsHC-volcano.png",
                    save.dir,
                    logFC_cutoff=log2(1.5),
                    use_p_val=T,
                    colours=colours,
                    k=10)


