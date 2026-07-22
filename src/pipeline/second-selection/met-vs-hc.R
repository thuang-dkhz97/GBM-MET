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



dparg <- list(phenotype_file_name = "data/proteomic-phenotype.txt",
              read_count_dir_path = "src/pipeline/filtration-after-initial-selection/candidate",
              read_count_file_name = "met-vs-hc-candidate.csv",
              sep = ",",
              dataset_id = "MET_vs_HC_combat_mod_2nd",
              classification_criteria = "METVsHC",
              classes = c("HC", "MET"),
              cores = 16,
              results_dir_path = "results/second-run/fem_pipeline_results_combat_mod_MET_VS_HC",
              norm = "none", perform_filter = FALSE,
              fems_to_run = c("all",
                              "t-test", "wilcoxontest",
                              "ranger_pos_impu_cor",
                              "mrmr30", "mrmr50",
                              "mrmr75", "mrmr100",
                              "mrmr_perc50",
                              "RF_RFE",
                              "ga_rf"))

do.call(execute_pipeline, dparg)
