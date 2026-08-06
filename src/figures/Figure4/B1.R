suppressPackageStartupMessages(library(ggplot2))


plot.boxplot <- function(data, phenotype, best_markers, conditions,
                         plot.save.dir, colours) {
  label <- phenotype$Label[match(colnames(data), phenotype$Sample)]
  
  best_markers <- best_markers %>%
    column_to_rownames("Molecule")
  
  # print(colnames(best_markers))
  # From data, extract the best markers
  # print(rownames(data))
  full.data <- data
  data <- data[rownames(data) %in% rownames(best_markers), ]
  
  # print(t(data))
  data <- t(data) %>% as.data.frame
  data$group <- label
  
  data_to_plot <- data %>%
    pivot_longer(cols=!group, names_to="biomarker", values_to="norm_expr") 
  
  biomarker.agg <- data_to_plot %>%
    group_by(biomarker) %>%
    summarise(med_expr=median(norm_expr)) %>% 
    arrange(desc(med_expr))
  
  data_to_plot <- data_to_plot %>% 
    mutate(biomarker = factor(biomarker, biomarker.agg$biomarker))
  
  
  ggplot(data_to_plot, aes(x = biomarker, 
                                y = log(norm_expr),
                                fill = group)) +
    geom_boxplot(size = 0.2, alpha = 0.99, outlier.alpha=0.5) +
    # ggtitle(paste(sub("Vs", " Vs ", sub("PREOPE", "GBM", comparison)), best_str, x_lab)) +
    ylab("Quantile normalized expression") +
    ggplot2::scale_fill_manual(name = "Condition", values = colours) +
    guides(fill = guide_legend(title = "Condition")) +
    theme(axis.text.x = element_text(size=rel(1.2), angle = 90),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.3)),
          axis.title.y = element_text(size=rel(1.3)),
          plot.title  = element_text(size=rel(1.5)))
  comp.name <- paste(rev(conditions), collapse="-vs-")
  R.utils::mkdirs(file.path(plot.save.dir, comp.name))
  ggsave(file.path(plot.save.dir, comp.name, "boxplot.png"), units="cm", width=30)
  
  
}


# ============ GBM V MET
save.dir <- file.path("src", "figures", "Figure4")

gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)

gbm.v.met.best.marker <- read.csv("src/DE_pipeline-2.8/post-RFE-select-by-AUC/select-AUC/GBM-vs-MET-final-sel.csv")

plot.boxplot(gbm.v.met.processed, gbm.v.met.labels,
             gbm.v.met.best.marker,
             c("MET", "GBM"),
             save.dir, 
             colours=c("#4cb753", "#61a5fe"))





