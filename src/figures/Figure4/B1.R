suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(umap))

plot.figs <- function(data, phenotype, cohort.labels, best_markers, conditions,
                         plot.save.dir, colours) {
  label <- phenotype$Label[match(colnames(data), phenotype$Sample)]
  cohort.matched <- cohort.labels$Cohort[match(colnames(data), cohort.labels$Sample.ID)]
  
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
                                y = norm_expr,
                                fill = group)) +
    geom_boxplot(size = 0.2, alpha = 0.99, outlier.alpha=0.5) +
    # ggtitle(paste(sub("Vs", " Vs ", sub("PREOPE", "GBM", comparison)), best_str, x_lab)) +
    ylab("Quantile normalized expression") +
    ggplot2::scale_fill_manual(name = "Condition", values = colours) +
    guides(fill = guide_legend(title = "Condition")) +
    theme_bw() +
    theme(axis.text.x = element_text(size=rel(1.2), angle = 90),
          axis.text.y = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.3)),
          axis.title.y = element_text(size=rel(1.3)),
          plot.title  = element_text(size=rel(1.5)))
  comp.name <- paste(rev(conditions), collapse="-vs-")
  R.utils::mkdirs(file.path(plot.save.dir, comp.name))
  ggsave(file.path(plot.save.dir, comp.name, "boxplot.png"), units="cm", width=30,
         height=17)
  
  # Then apply UMAP to the entire feature set
  full.result <- umap(t(full.data), n_neighbors=15)
  full.dim_red_df <- data.frame(x=full.result$layout[,1], y=full.result$layout[,2],
                                group=data$group, cohort=cohort.matched)
  ggplot(full.dim_red_df, aes(x=x, y=y, fill=group, shape=cohort)) +
    geom_point(size=8, stroke=2) +
    scale_shape_manual(name = "Data Cohort", values = c(21, 22)) +
    scale_fill_manual(name = "Condition", values = colours) +
    # guides(fill = guide_legend(override.aes = list(shape = 21, colour = colours))) +
    theme_bw() +
    theme(panel.border=element_rect(color = "gray10", fill = NA, linewidth = 1.2),
          legend.position="none",
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank())
          # axis.text=element_blank())
  ggsave(file.path(plot.save.dir, comp.name, "full-umap.png"), units="cm", width=30, height=21)
  
  # And also to the best set
  result <- umap(data %>% dplyr::select(!all_of(c("group"))),
                 n_neighbors=15)
  dim_red_df <- data.frame(x=result$layout[,1], y=result$layout[,2],
                           group=data$group, cohort=cohort.matched)
  # print(dim_red_df)
  
  ggplot(dim_red_df, aes(x=x, y=y, fill=group, shape=cohort)) +
    geom_point(size=8, stroke=2) +
    scale_shape_manual(name = "Data Cohort", values = c(21, 22)) +
    scale_fill_manual(name = "Condition", values = colours) +
    # guides(fill = guide_legend(override.aes = list(shape = 21, colour = colours))) +
    theme_bw() +
    theme(panel.border=element_rect(color = "gray10", fill = NA, linewidth = 1.2),
          legend.position="none",
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank())
          # axis.text=element_blank())
  ggsave(file.path(plot.save.dir, comp.name, "select-umap.png"), units="cm", width=30,height=21)
  
}

save.dir <- file.path("src", "figures", "Figure4")
cohort.labels <- read.xlsx("./data/class-labels.xlsx")

# ============ GBM V MET

gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)

gbm.v.met.best.marker <- read.csv("src/DE_pipeline-2.8/post-RFE-select-by-AUC/select-AUC/GBM-vs-MET-final-sel.csv")

plot.figs(gbm.v.met.processed, gbm.v.met.labels,
           cohort.labels,
           gbm.v.met.best.marker,
           c("MET", "GBM"),
           save.dir, 
           colours=c("#61a5fe", "#4cb753"))

# ============== GBM V HC

gbm.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-hc-norm.csv",
                                row.names=1)
gbm.v.hc.labels <- read.table("src/process-2.6/combine-data-3/GBMVsHC_output_labels.txt", header=T)

gbm.v.hc.best.marker <- read.csv("src/DE_pipeline-2.8/post-RFE-select-by-AUC/select-AUC/GBM-vs-HC-final-sel.csv")

plot.figs(gbm.v.hc.processed, gbm.v.hc.labels,
          cohort.labels,
          gbm.v.hc.best.marker,
          c("HC", "GBM"),
          save.dir, 
          colours=c("#63A5FE", "#FF605A"))

# ============== MET V HC

met.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/met-v-hc-norm.csv",
                               row.names=1)
met.v.hc.labels <- read.table("src/process-2.6/combine-data-3/METVsHC_output_labels.txt", header=T)

met.v.hc.best.marker <- read.csv("src/DE_pipeline-2.8/post-RFE-select-by-AUC/select-AUC/MET-vs-HC-final-sel.csv")

plot.figs(met.v.hc.processed, met.v.hc.labels,
          cohort.labels,
          met.v.hc.best.marker,
          c("HC", "MET"),
          save.dir, 
          colours=c("#4cb753", "#FF605A"))







