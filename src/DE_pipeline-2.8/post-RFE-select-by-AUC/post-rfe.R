
get_RFE_selected_features <- function(ranked_features, 
                                      RFE_results_001_file_path,
                                      RFE_results_003_file_path,
                                      RFE_results_005_file_path,
                                      output_file_path){
  ranked_features <- ranked_features %>%
    mutate(combined_score_rank = rank(combined_score), .after = Molecule)
  ranked_features <- ranked_features %>%
    mutate(combined_score_rank = nrow(ranked_features) - combined_score_rank + 1)
  
  q75 <- quantile(ranked_features$combined_score, 0.75)
  ranked_features <- ranked_features %>% filter(combined_score >= q75)
  
  # if(filter_low_scored_features){
  #   q75 <- quantile(ranked_features$combined_score, 0.75)
  #   ranked_features <- ranked_features %>% filter(combined_score >= q75)
  # }
  
  RFE_results_001 <- read.csv(RFE_results_001_file_path)
  RFE_results_003 <- read.csv(RFE_results_003_file_path)
  RFE_results_005 <- read.csv(RFE_results_005_file_path)
  
  best_iter_001 <- get_best_iter(RFE_results_001)
  best_iter_003 <- get_best_iter(RFE_results_003)
  best_iter_005 <- get_best_iter(RFE_results_005)
  
  if(best_iter_001[2] >= best_iter_003[2]){
    if(best_iter_001[2] >= best_iter_005[2]){
      RFE_results <- RFE_results_001
      best_iter <- best_iter_001
      best_cutoff <- 0.01
      print("0.01 best")
    } else{
      RFE_results <- RFE_results_005
      best_iter <- best_iter_005
      best_cutoff <- 0.05
      print("0.05 best")
    }
  } else{
    if(best_iter_003[2] >= best_iter_005[2]){
      RFE_results <- RFE_results_003
      best_iter <- best_iter_003
      best_cutoff <- 0.03
      print("0.03 best")
    } else{
      RFE_results <- RFE_results_005
      best_iter <- best_iter_005
      print("0.05 best")
      best_cutoff <- 0.05
    }
  }
  
  RFE_results.filt <- RFE_results %>%
    filter(iter <= best_iter[1] & feature_removed)
  
  ranked_features <- ranked_features %>%
    filter(! Molecule %in% RFE_results.filt$feature)
  
  ranked_features <- ranked_features %>%
    mutate("best_cutoff" = best_cutoff, .after = Molecule) %>%
    mutate("meanAUC" = best_iter[2], .after = best_cutoff) %>%
    mutate("best_iter" = best_iter[1], .after = meanAUC)
  
  write.csv(ranked_features, output_file_path, row.names = FALSE)
}


get_best_iter <- function(RFE_results){
  best_meanAUC <- 0
  least_num_features <- RFE_results[1, "pos_logFC_feature_count"] + RFE_results[1, "non_pos_logFC_feature_count"]
  best_iter <- 0
  for(i in c(1:nrow(RFE_results))){
    meanAUC <- RFE_results[i, "MeanAUC"]
    num_features <- RFE_results[i, "pos_logFC_feature_count"] + RFE_results[i, "non_pos_logFC_feature_count"]
    if(num_features <= 25){
      if((meanAUC > best_meanAUC) | (meanAUC == best_meanAUC & num_features < least_num_features)){
        best_iter <- RFE_results[i, "iter"]
        best_meanAUC <- meanAUC
        least_num_features <- num_features
      }
    }
  }
  return(c(best_iter, best_meanAUC))
}


# ============ GBM V MET ========================
save.dir <- file.path("src", "DE_pipeline-2.8", "post-RFE-select-by-AUC", "select-AUC")
R.utils::mkdirs(save.dir)
gbm.v.met.ranking <- read.csv("./src/DE_pipeline-2.8/compute-feature-ranks/feature-ranks/GBM-vs-MET-de.csv")
gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)

gbm.v.met.001.path <- "src/DE_pipeline-2.8/RFE/rfe-results/GBM-vs-MET-rfe-AUC_decrease-0.01-pos_logFC_prop-0.25.csv"
gbm.v.met.003.path <- "src/DE_pipeline-2.8/RFE/rfe-results/GBM-vs-MET-rfe-AUC_decrease-0.03-pos_logFC_prop-0.25.csv"
gbm.v.met.005.path <- "src/DE_pipeline-2.8/RFE/rfe-results/GBM-vs-MET-rfe-AUC_decrease-0.05-pos_logFC_prop-0.25.csv"

get_RFE_selected_features(gbm.v.met.ranking, 
                          gbm.v.met.001.path,
                          gbm.v.met.003.path,
                          gbm.v.met.005.path,
                          file.path(save.dir, "GBM-vs-MET-final-sel.csv"))

# =============== GBM V HC

gbm.v.hc.ranking <- read.csv("./src/DE_pipeline-2.8/compute-feature-ranks/feature-ranks/GBM-vs-HC-de.csv")
gbm.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-hc-norm.csv",
                                row.names=1)
gbm.v.hc.labels <- read.table("src/process-2.6/combine-data-3/GBMVsHC_output_labels.txt", header=T)

gbm.v.hc.001.path <- "src/DE_pipeline-2.8/RFE/rfe-results/GBM-vs-HC-rfe-AUC_decrease-0.01-pos_logFC_prop-0.25.csv"
gbm.v.hc.003.path <- "src/DE_pipeline-2.8/RFE/rfe-results/GBM-vs-HC-rfe-AUC_decrease-0.03-pos_logFC_prop-0.25.csv"
gbm.v.hc.005.path <- "src/DE_pipeline-2.8/RFE/rfe-results/GBM-vs-HC-rfe-AUC_decrease-0.05-pos_logFC_prop-0.25.csv"

get_RFE_selected_features(gbm.v.hc.ranking, 
                          gbm.v.hc.001.path,
                          gbm.v.hc.003.path,
                          gbm.v.hc.005.path,
                          file.path(save.dir, "GBM-vs-HC-final-sel.csv"))

# =============== MET V HC

met.v.hc.ranking <- read.csv("./src/DE_pipeline-2.8/compute-feature-ranks/feature-ranks/MET-vs-HC-de.csv")
met.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/met-v-hc-norm.csv",
                               row.names=1)
met.v.hc.labels <- read.table("src/process-2.6/combine-data-3/METVsHC_output_labels.txt", header=T)

met.v.hc.001.path <- "src/DE_pipeline-2.8/RFE/rfe-results/MET-vs-HC-rfe-AUC_decrease-0.01-pos_logFC_prop-0.25.csv"
met.v.hc.003.path <- "src/DE_pipeline-2.8/RFE/rfe-results/MET-vs-HC-rfe-AUC_decrease-0.03-pos_logFC_prop-0.25.csv"
met.v.hc.005.path <- "src/DE_pipeline-2.8/RFE/rfe-results/MET-vs-HC-rfe-AUC_decrease-0.05-pos_logFC_prop-0.25.csv"

get_RFE_selected_features(met.v.hc.ranking, 
                          met.v.hc.001.path,
                          met.v.hc.003.path,
                          met.v.hc.005.path,
                          file.path(save.dir, "MET-vs-HC-final-sel.csv"))
