
library(stringr)
library(ComplexHeatmap)
library(viridis)
library(tibble)

plot_heatmap <- function(results_dir = "fem_pipeline_results",
                         best_ML_path = "",
                         save_file_path = "",
                         fsm.show,
                         fsm.add.on.names){
  # data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
  #                         sep = ',', header = TRUE)
  # fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
  #                        sep = ',', header = TRUE)
  # model_results <- read.table(paste(results_dir, "model_results_test.csv", sep = "/"),
  #                             sep = ',', header = TRUE)
  all.results.file <- list.dirs(results_dir)
  best.sel <- readRDS(best_ML_path)
  
  # For each file, open the underlying meanAUC file.
  model.results <- data.frame(matrix(data=NA, nrow=0, ncol=3))
  colnames(model.results) <- c("fsm", "model", "mean AUC")
  for (dir.name in all.results.file[2:length(all.results.file)]) {
    # Extract the fsm and number
    print(dir.name)
    
    fsm.size <- str_extract(dir.name, "(?<=num-).*")
    fsm.name <- str_extract(dir.name, "(?<=fsm-).*?(?=-feat)")
    # Open the file
    test.res <- read.csv(file.path(dir.name, 'model_results_test.csv'), header=T)
    
    fsm.name <- paste0(fsm.name, "_", fsm.size)
    # fsm.name <- paste0(fsm.name, fsm.add.on.names)
    model.results <- rbind(model.results, 
                           data.frame(fsm.name, test.res[,which(colnames(test.res) == "Model" |
                                                         colnames(test.res) == "Mean_AUC")]))
  }
  best.sel.name <- paste0(best.sel$method, "_", best.sel$sel.size)
  # best.sel.model <- best.sel$model
  print(best.sel.name)
  # Choose models to display
  model.results <- model.results[c(which(model.results$fsm.name %in% fsm.show),
                                   which(model.results$fsm.name==best.sel.name)),]

  # Now, convert to wide format
  data_to_plot <- model.results %>% 
    pivot_wider(names_from="fsm.name",
                values_from="Mean_AUC") %>%
    column_to_rownames("Model") %>%
    as.matrix
  colnames(data_to_plot) <- paste0(colnames(data_to_plot), fsm.add.on.names)
  print(data_to_plot)
  
  file_name_curr <- save_file_path
  png(file_name_curr, units = "cm", width = 20, height = 15, res = 1200)
  ht <- Heatmap(data_to_plot, name = "Mean AUC",
                col = magma(5),
                rect_gp = gpar(col = "white", lwd = 1),
                cluster_columns = T,
                # column_order = replaced_dataset_id_vec[replaced_dataset_id_vec %in% colnames(data_to_plot)],
                column_names_rot = 30,
                column_names_gp = gpar(fontsize = 10),
                row_title = "Classification Models",
                row_names_side = "left",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.4f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 10, col = "slateblue3"))
                })
  draw(ht, column_title = paste("Mean AUC results from selected common features"))
  dev.off()
  
}

R.utils::mkdirs(file.path("./src/figures/SuppFigure6/images/"))

save.file <- "./src/figures/SuppFigure6/images/GBM_VS_MET.png"
plot_heatmap("./results/second-run/GBM_VS_MET/",
             "src/pipeline-2.7/final-filtration/selected/gbm-vs-met-sels.rds",
             save.file,
             c("mrmr75_30",
               "mrmr30_30",
               "RF_RFE_28", "mrmr50_29"),
                "_GBMVsMET")



save.file <- "./src/figures/SuppFigure6/images/GBM_VS_HC.png"
plot_heatmap("./results/second-run/GBM_VS_HC/",
             "src/pipeline-2.7/final-filtration/selected/gbm-vs-hc-sels.rds",
             save.file,
             c("ranger_pos_impu_cor_29",
               "mrmr100_27",
               "RF_RFE_28", "mrmr50_29"),
             "_GBMVsHC")


save.file <- "./src/figures/SuppFigure6/images/MET_VS_HC.png"
plot_heatmap("./results/second-run/MET_VS_HC/",
             "src/pipeline-2.7/final-filtration/selected/met-vs-hc-sels.rds",
             save.file,
             c("mrmr100_28",
               "RF_RFE_27", "mrmr100_30"),
             "_METVsHC")

