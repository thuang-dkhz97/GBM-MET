library(edgeR)
library(tidyverse)

source('src/DE_pipeline-2.8/cm_rf.R')

de.ml.biomarker.discovery.compute.feature.rank <- function(data, phenotype,
                                                           conditions, pval_cutoff,
                                                           best_features_ml,
                                                           result_file_path) {
  
  # Sort the labels.
  label <- phenotype$Label[match(colnames(data), phenotype$Sample)]

  # Extract features selected from 2.7
  best_features <- best_features_ml$features
  
  # Filter data
  model_matrix <- model.matrix(~ 0 + label)
  colnames(model_matrix) <- sub("label", "", colnames(model_matrix), fixed = TRUE)
  keep <- edgeR::filterByExpr(data, group = model_matrix)
  data <- data[keep, ]
  
  set.seed(1000)
  k_val <- 5
  times_val <- 6
  k_prod_times <- k_val * times_val
  train_index <- caret::createMultiFolds(y = label, k = k_val, times = times_val)
  
  for (i in c(1:k_prod_times)) {
    
    # i <- 1
    label.train <- label[train_index[[i]]]
    label.test <- label[-train_index[[i]]]
    # print(head(label.train$Sample))
    # print(head(label.test$Sample))
    
    # print(summary(factor(label.train$PREOPE_MET_HC)))
    # print(summary(factor(label.test$PREOPE_MET_HC)))
    
    data.train <- data[,train_index[[i]]]
    data.test <- data[,-train_index[[i]]]
    
    result.de <- perform_de(data.train, label.train, conditions) %>%
      mutate("iter" = i) %>%
      mutate("significant" = case_when((PVal < pval_cutoff) ~ 1,
                                       TRUE ~ 0))
    
    ##########################################
    # label.train <- label.train
    # label.test <- label.test
    # data.train <- data[, label.train$Sample]
    # data.test <- data[, label.test$Sample]
    
    result_list <- rf_model(as.data.frame(t(data.train)), label.train,
                            as.data.frame(t(data.test)), label.test,
                            conditions, classifier_feature_imp = TRUE)
    feature_importance <- result_list[[2]]
    feature_importance <- feature_importance %>%
      arrange(desc(MeanDecreaseGini)) %>%
      dplyr::rename(c("Molecule" = "feature"))
    
    result <- result.de %>%
      inner_join(feature_importance)
    ##########################################
    
    #obtain features selected by ranger to count number of times feature is chosen - to be used as stability score
    
    ranger_model <- ranger::ranger(x = as.data.frame(t(data.train)),
                                   y = factor(label.train), importance = "impurity_corrected")
    features <- which(ranger_model$variable.importance > 0)
    
    feature_chosen <- data.frame(Molecule = rownames(data.train))
    feature_chosen <- feature_chosen %>%
      mutate(ranger_chosen = ifelse(Molecule %in% names(features), 1, 0))
    
    result <- result %>%
      inner_join(feature_chosen)
    
    if(i == 1){
      result_all <- result
    } else{
      result_all <- rbind(result_all, result)
    }
    
  }
  
  result.agg <- result_all %>%
    group_by(Molecule) %>%
    summarize(sig_count = sum(significant),
              fischers_combined_p = -2*sum(log(PVal)),
              mean_logFC = mean(logFC),
              MeanDecreaseGini_mean = mean(MeanDecreaseGini),
              stability = sum(ranger_chosen)) %>%
    mutate(identified_biomarker = case_when(Molecule %in% best_features ~ 1,
                                            TRUE ~ 0))
  
  
  #obtain ranking and combined_score
  result.agg <- result.agg %>%
    mutate(rank.fischers_combined_p = rank(fischers_combined_p)) %>%
    mutate(rank.abs_mean_logFC = rank(abs(mean_logFC))) %>%
    mutate(rank.MeanDecreaseGini_mean = rank(MeanDecreaseGini_mean)) %>%
    mutate(rank.stability = rank(stability)) %>%
    mutate(score.identified_biomarker = 2 * nrow(result.agg) * identified_biomarker) #twice the weightage for selected biomarkers
  
  result.agg <- result.agg %>%
    mutate(combined_score = rank.fischers_combined_p + rank.abs_mean_logFC +
             rank.MeanDecreaseGini_mean + rank.stability +
             score.identified_biomarker,
           .after = Molecule) %>%
    arrange(desc(combined_score))
  
  R.utils::mkdirs(result_file_path)
  res.file.name <- paste(rev(conditions), collapse="-vs-")
  write.csv(result.agg, file.path(result_file_path, paste0(res.file.name, "-de.csv")), row.names = FALSE)
}

perform_de <- function(data, label, conditions){
  model_matrix <- model.matrix(~ 0 + label)
  colnames(model_matrix) <- sub("label", "", colnames(model_matrix), fixed = TRUE)

  contr_matrix <- makeContrasts(contrasts = paste0(conditions[2], " - ", conditions[1]),
                                levels = colnames(model_matrix))
  
  fit <- lmFit(log2(data), model_matrix)
  # print(head(coef(fit)))
  fit <- contrasts.fit(fit, contr_matrix)
  # print(head(coef(fit)))
  efit <- eBayes(fit)
  top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
    rownames_to_column("molecule")
  result <- top.table %>%
    dplyr::select(molecule, logFC, P.Value, adj.P.Val) %>%
    dplyr::rename(Molecule = molecule, adjPVal = adj.P.Val, PVal = P.Value) %>%
    arrange(logFC)
  
  return(result)
}


save.dir <- file.path("src", "DE_pipeline-2.8", "feature-ranks")

# Conduct combined feature selection with DE and data file path
gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
# TEMPORARILY REMOVE SB12.1 - TO BE UPDATED WHEN SUSANNAH RESPONDS TO ME
gbm.v.met.processed <- gbm.v.met.processed %>% dplyr::select(-all_of("SB12.1"))
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)
gbm.v.met.ml.pipeline.sels <- readRDS("src/pipeline-2.7/final-filtration/selected/gbm-vs-met-sels.rds")
de.ml.biomarker.discovery.compute.feature.rank(gbm.v.met.processed,
                                               gbm.v.met.labels,
                                               c("MET", "GBM"),
                                               0.05,
                                               gbm.v.met.ml.pipeline.sels,
                                               file.path(save.dir))

# ========= GBM V HC

gbm.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-hc-norm.csv",
                               row.names=1)
gbm.v.hc.processed <- gbm.v.hc.processed %>% dplyr::select(-all_of("SB12.1"))
gbm.v.hc.labels <- read.table("src/process-2.6/combine-data-3/GBMVsHC_output_labels.txt", header=T)
gbm.v.hc.ml.pipeline.sels <- readRDS("src/pipeline-2.7/final-filtration/selected/gbm-vs-hc-sels.rds")
de.ml.biomarker.discovery.compute.feature.rank(gbm.v.hc.processed,
                                               gbm.v.hc.labels,
                                               c("HC", "GBM"),
                                               0.05,
                                               gbm.v.hc.ml.pipeline.sels,
                                               file.path(save.dir))

# ============ MET V HC

met.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/met-v-hc-norm.csv",
                               row.names=1)
met.v.hc.labels <- read.table("src/process-2.6/combine-data-3/METVsHC_output_labels.txt", header=T)
met.v.hc.ml.pipeline.sels <- readRDS("src/pipeline-2.7/final-filtration/selected/met-vs-hc-sels.rds")
de.ml.biomarker.discovery.compute.feature.rank(met.v.hc.processed,
                                               met.v.hc.labels,
                                               c("HC", "MET"),
                                               0.05,
                                               met.v.hc.ml.pipeline.sels,
                                               file.path(save.dir))

