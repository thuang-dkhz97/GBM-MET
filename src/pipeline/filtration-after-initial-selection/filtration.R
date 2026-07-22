
library(dplyr)
library(scales)

# NEED TO REDO THIS TO MAKE FILES FOR COUNTS OF 27, 28, 29, 30
obtain.candidate <- function(feature.mat, perform.mat) {
  # Remove metrics related to logistic regression, as there is likely overfitting.
  perform.mat <- perform.mat[perform.mat$Model != "Simple logistic regression",]
 
  # First, determine which feature selection had the highest performance
  # Group by FSM
  grouped.perform <- perform.mat %>% group_by(FSM)
  
  # Take the mean of a column
  grouped.perform.auc.mean <- grouped.perform %>%
    summarise(mean_auc = mean(Mean_AUC))
  print(grouped.perform.auc.mean)
  # Obtain fsm with highest mean auc not including "all" and "ga_rf
  grouped.perform.auc.mean <- grouped.perform.auc.mean[grouped.perform.auc.mean$FSM != "all" &
                                                         grouped.perform.auc.mean$FSM != "ga_rf",]
  opt.fsm <- grouped.perform.auc.mean$FSM[which(grouped.perform.auc.mean$mean_auc == max(grouped.perform.auc.mean$mean_auc))]
  print(opt.fsm)
  # Obtain the counts for the best AUC
  best.fsm <- feature.mat[which(feature.mat$FSM %in% opt.fsm),-c(1, 2)]
  # count how many times each feature was selected and choose those with
  # selection >=27
  # If the optimal fsm is ga_rf, scale it so that the highest selected
  # feature is 30 - ga_rf tends to select under 27. - Not used
  total.counts <- colSums(best.fsm)
  # if (opt.fsm == "ga_rf") {
  #   total.counts <- rescale(total.counts, to=c(min(total.counts), 30))
  # }
  
  candidate <- colnames(best.fsm)[total.counts >= 27]
  print(candidate)
  candidate
}

gbm.vs.met.perform.mat <- read.csv("./results/initial-run/fem_pipeline_results_combat_mod_GBM_VS_MET/model_results_test.csv")
gbm.vs.met.feature.mat <- read.csv("./results/initial-run/fem_pipeline_results_combat_mod_GBM_VS_MET/GBM_vs_MET_combat_mod_GBMVsMET_features.csv")
gbm.vs.met.candidate <- obtain.candidate(gbm.vs.met.feature.mat, gbm.vs.met.perform.mat)

gbm.vs.hc.perform.mat <- read.csv("./results/initial-run/fem_pipeline_results_combat_mod_GBM_Vs_HC/model_results_test.csv")
gbm.vs.hc.feature.mat <- read.csv("./results/initial-run/fem_pipeline_results_combat_mod_GBM_Vs_HC/GBM_vs_HC_combat_mod_GBMVsHC_features.csv")
gbm.vs.hc.candidate <- obtain.candidate(gbm.vs.hc.feature.mat, gbm.vs.hc.perform.mat)
# Use the candidates to change the counts and save

met.vs.hc.perform.mat <- read.csv("./results/initial-run/fem_pipeline_results_combat_mod_MET_Vs_HC/model_results_test.csv")
met.vs.hc.feature.mat <- read.csv("./results/initial-run/fem_pipeline_results_combat_mod_MET_Vs_HC/MET_vs_HC_combat_mod_METVsHC_features.csv")
met.vs.hc.candidate <- obtain.candidate(met.vs.hc.feature.mat, met.vs.hc.perform.mat)

# Use the candidates to change the counts and save
gbm.vs.met.og.counts <- read.csv("src/process/combine-data-3/gbm-v-met-norm.csv")
gbm.vs.met.candidate.counts <- gbm.vs.met.og.counts[gbm.vs.met.og.counts[,1] %in% gbm.vs.met.candidate,]
save.dir <- "./src/pipeline/filtration-after-initial-selection/candidate"
R.utils::mkdirs(save.dir)
write.csv(gbm.vs.met.candidate.counts, file.path(save.dir, "gbm-vs-met-candidate.csv"))
# Use the candidates to change the counts and save
gbm.vs.hc.og.counts <- read.csv("src/process/combine-data-3/gbm-v-hc-norm.csv")
gbm.vs.hc.candidate.counts <- gbm.vs.hc.og.counts[gbm.vs.hc.og.counts[,1] %in% gbm.vs.hc.candidate,]
write.csv(gbm.vs.hc.candidate.counts, file.path(save.dir, "gbm-vs-hc-candidate.csv"))
# Use the candidates to change the counts and save
met.vs.hc.og.counts <- read.csv("src/process/combine-data-3/met-v-hc-norm.csv")
met.vs.hc.candidate.counts <- met.vs.hc.og.counts[met.vs.hc.og.counts[,1] %in% met.vs.hc.candidate,]
write.csv(met.vs.hc.candidate.counts, file.path(save.dir, "met-vs-hc-candidate.csv"))

