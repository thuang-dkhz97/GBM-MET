library(dplyr)
library(FEMPipeline)
# source("FEMPipeline/R/pipeline.R")
# source("FEMPipeline/R/extract.R")
# source("FEMPipeline/R/fem_t_test.R")
# source("FEMPipeline/R/fem_wilcoxon_test.R")
# source("FEMPipeline/R/fem_rfrfe.R")
# source("FEMPipeline/R/fem_ga.R")
# source("FEMPipeline/R/fem_pca_transformation.R")
# source("FEMPipeline/R/fem_rf_features.R")
# source("FEMPipeline/R/fem_ranger_features.R")
# source("FEMPipeline/R/fem_mrmr_features.R")
# source("FEMPipeline/R/fem_phate_transformation.R")
# source("FEMPipeline/R/fem_umap_transformation.R")
# source("FEMPipeline/R/fem_plsda_transformation.R")
# source("FEMPipeline/R/fem_kpca_transformation.R")
# source("FEMPipeline/R/fem_mrmr_plsda_transformation.R")
# source("FEMPipeline/R/fem_ranger_plsda.R")
# source("FEMPipeline/R/fem_ranger_pos_features.R")
# source("FEMPipeline/R/fem_ipage_transformation.R")
# source("FEMPipeline/R/feature_extraction_arguments.R")
# source("FEMPipeline/R/run_fsm_and_models.R")
# source("FEMPipeline/R/preprocessing.R")
# source("FEMPipeline/R/preprocessing_norm.R")
# source("FEMPipeline/R/preprocessing_filter.R")
# source("FEMPipeline/R/run_all_models.R")
# source("FEMPipeline/R/cm_logistic_regression.R")
# source("FEMPipeline/R/cm_rf.R")
# source("FEMPipeline/R/cm_svm.R")
# source("FEMPipeline/R/compute_metrics.R")

# NEED TO CHANGE FEMS_TO_RUN = ALL ONLY, AND NEED TO RUN WITH 
# T-TEST PANEL WITH >=27, >=28, >=29, etc.
# AND REPEATS
methods <- c("all", "ga_rf", "mrmr_perc50", "mrmr30",
             "mrmr50", "mrmr75", "mrmr100", "ranger_pos_impu_cor",
             "RF_RFE", "t-test", "wilcoxontest")
# methods <- c("t-test")
for (i in 27:30) {
  # Then for each method, extract all features that were selected
  # in more than or equal to i iterations.
  for (method in methods) {
    if (method == "all" & i %in% c(27, 28, 29)) next
    dir.path <- paste0("results/second-run/GBM_VS_MET/fem_pipeline_results_combat_mod_GBM_VS_MET-fsm-", method, 
                       "-feat-num-", i)
    R.utils::mkdirs(dir.path)
    
    dparg <- list(phenotype_file_name = "data/proteomic-phenotype.txt",
                  read_count_dir_path = "src/pipeline-2.7/filtration-after-initial-selection/candidate",
                  read_count_file_name = paste0("gbm-vs-met-fsm-", method, "-feat-num-", 
                                                i, "-candidate.csv"),
                  sep = ",",
                  dataset_id = "GBM_vs_MET_combat_mod_2nd",
                  classification_criteria = "GBMVsMET",
                  classes = c("MET", "GBM"),
                  cores = 16,
                  results_dir_path = dir.path,
                  norm = "none", perform_filter = FALSE,
                  fems_to_run = c("all"))
    
    do.call(execute_pipeline, dparg)
  }
}




# dparg <- list(phenotype_file_name = "data/proteomic-phenotype.txt",
#               read_count_dir_path = "src/pipeline/filtration-after-initial-selection/candidate",
#               read_count_file_name = "gbm-vs-met-candidate.csv",
#               sep = ",",
#               dataset_id = "GBM_vs_MET_combat_mod_2nd",
#               classification_criteria = "GBMVsMET",
#               classes = c("MET", "GBM"),
#               cores = 16,
#               results_dir_path = "results/second-run/fem_pipeline_results_combat_mod_GBM_VS_MET",
#               norm = "none", perform_filter = FALSE,
#               fems_to_run = c("all",
#                               "t-test", "wilcoxontest",
#                               "ranger_pos_impu_cor",
#                               "mrmr30", "mrmr50",
#                               "mrmr75", "mrmr100",
#                               "mrmr_perc50",
#                               "RF_RFE",
#                               "ga_rf"))
# 
# do.call(execute_pipeline, dparg)
