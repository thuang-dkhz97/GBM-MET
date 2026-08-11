
suppressPackageStartupMessages(library(dplyr))

compute_mean_AUC <- function(data, label, conditions){
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
    
    data.train <- data[, train_index[[i]]]
    data.test <- data[, -train_index[[i]]]
    
    ##########################################
    # label.train <- label.train %>%
    #   filter(!is.na(Label))
    # label.test <- label.test %>%
    #   filter(!is.na(Label))
    # data.train <- data[, train_index[[i]]]
    # data.test <- data[, -train_index[[i]]]
    
    result <- rf_model(as.data.frame(t(data.train)), label.train,
                       as.data.frame(t(data.test)), label.test,
                       conditions)
    result.test <- result %>%
      filter(Type == "test")
    
    pr <- ROCR::prediction(result.test$Pred_prob, result.test$TrueLabel, label.ordering = conditions)
    auc.test <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
    
    if(i == 1){
      result_all <- auc.test
    } else{
      result_all <- c(result_all, auc.test)
    }
  }
  mean_AUC <- mean(result_all)
  return (mean_AUC)
}

conduct.rfe <- function(data, phenotype, conditions,
                        rankings, results_file_path,
                        meanAUC_decrease_cutoff=0.01,
                        pos_logFC_min_proportion = 0.25) {
  label <- phenotype$Label[match(colnames(data), phenotype$Sample)]
  
  ranked_features <- rankings %>% 
    column_to_rownames("Molecule") %>%
    arrange(combined_score)
  
  # Store performance of every protein used
  meanAUC <- compute_mean_AUC(data, label, conditions)
  feature_set <- rownames(ranked_features)
  feature_set.pos_logFC <- rownames(ranked_features %>%               #NOTE : this list of features is ranked from lowest score to highest
                                      filter(mean_logFC > 0))
  feature_set.non_pos_logFC <- setdiff(feature_set, feature_set.pos_logFC)
  
  results <- data.frame(iter = 0, MeanAUC = meanAUC, 
                        pos_logFC_feature_count = length(feature_set.pos_logFC),
                        non_pos_logFC_feature_count = length(feature_set.non_pos_logFC),
                        feature = NA, feature_removed = NA)
  
  # Extract the top 25%
  q75 <- quantile(ranked_features$combined_score, 0.75)
  ranked_features <- ranked_features %>% filter(combined_score >= q75)
  
  data <- data[rownames(ranked_features),]
  meanAUC <- compute_mean_AUC(data, label, conditions)
  feature_set <- rownames(ranked_features)
  feature_set.pos_logFC <- rownames(ranked_features %>%               #NOTE : this list of features is ranked from lowest score to highest
                                      filter(mean_logFC > 0))
  feature_set.non_pos_logFC <- setdiff(feature_set, feature_set.pos_logFC)    
  
  results <- rbind(results,
                   data.frame(iter = 0, MeanAUC = meanAUC, 
                              pos_logFC_feature_count = length(feature_set.pos_logFC),
                              non_pos_logFC_feature_count = length(feature_set.non_pos_logFC),
                              feature = NA, feature_removed = NA))
  
  for(i in c(1:nrow(ranked_features))){
    # i <- 1
    if(i %% 50 == 1){
      print(paste("Iter", i, "out of", nrow(ranked_features)))
    }
    
    feature_considered <- rownames(ranked_features)[i]
    
    if(feature_considered %in% feature_set.pos_logFC){
      # print("Next")
      # print(length(feature_set.pos_logFC))
      # print(length(feature_set.non_pos_logFC))
      # print((length(feature_set.pos_logFC)-1)/(length(feature_set.pos_logFC)+length(feature_set.non_pos_logFC)-1))
      if(length(feature_set.pos_logFC) > 1 & (length(feature_set.pos_logFC)-1)/(length(feature_set.pos_logFC)+length(feature_set.non_pos_logFC)-1) > pos_logFC_min_proportion){
        
        current_iter_results <- data.frame(matrix(ncol = 6, nrow = 0,
                                                  dimnames = list(c(), 
                                                                  c('iter', 'MeanAUC', 
                                                                    'pos_logFC_feature_count', 'non_pos_logFC_feature_count',
                                                                    'feature', 'feature_removed'))))  
        while((length(feature_set.pos_logFC)-1)/(length(feature_set.pos_logFC)+length(feature_set.non_pos_logFC)-1) > pos_logFC_min_proportion){
          
          #select pos logFC feature with least score since previously a lower scored pos logfc 
          #               feature might not have been removed due to pos_logFC_min_proportion
          feature_considered.pos_logFC <- feature_set.pos_logFC[1]
          
          data_updated <- data[rownames(data) != feature_considered.pos_logFC, ]
          meanAUC_updated <- compute_mean_AUC(data_updated, label, conditions)
          
          feature_set.pos_logFC.updated <- setdiff(feature_set.pos_logFC, feature_considered.pos_logFC)
          
          if(meanAUC - meanAUC_updated < meanAUC_decrease_cutoff){ 
            feature_set.pos_logFC <- feature_set.pos_logFC.updated
            data <- data_updated
            meanAUC <- meanAUC_updated
            is_feature_removed <- TRUE
          } else{
            is_feature_removed <- FALSE
          }
          
          current_iter_results <- rbind(current_iter_results,
                                        data.frame(iter = i, MeanAUC = meanAUC, 
                                                   pos_logFC_feature_count = length(feature_set.pos_logFC),
                                                   non_pos_logFC_feature_count = length(feature_set.non_pos_logFC),
                                                   feature = feature_considered.pos_logFC, feature_removed = is_feature_removed))
          
          #while looping through pos_logFC features if we reach till current considered feature, then come out of the loop 
          #also come out of loop if current pos_logFC feature can't be removed due to meanAUC large decrease. 
          #       We assume that any higher scored features also can't be removed, and so exit the loop
          if(feature_considered.pos_logFC == feature_considered | !is_feature_removed){ 
            break
          }
        }
      } else{
        feature_considered <- feature_set.pos_logFC[1]
        is_feature_removed <- FALSE
        current_iter_results <- data.frame(iter = i, MeanAUC = meanAUC, 
                                           pos_logFC_feature_count = length(feature_set.pos_logFC),
                                           non_pos_logFC_feature_count = length(feature_set.non_pos_logFC),
                                           feature = feature_considered, feature_removed = is_feature_removed)
      }
    } else{
      data_updated <- data[rownames(data) != feature_considered, ]
      meanAUC_updated <- compute_mean_AUC(data_updated, label, conditions)
      feature_set.non_pos_logFC.updated <- setdiff(feature_set.non_pos_logFC, feature_considered)
      if(meanAUC - meanAUC_updated < meanAUC_decrease_cutoff){ 
        feature_set.non_pos_logFC <- feature_set.non_pos_logFC.updated   
        data <- data_updated
        meanAUC <- meanAUC_updated
        is_feature_removed <- TRUE
      } else{
        is_feature_removed <- FALSE
      }
      current_iter_results <- data.frame(iter = i, MeanAUC = meanAUC, 
                                         pos_logFC_feature_count = length(feature_set.pos_logFC),
                                         non_pos_logFC_feature_count = length(feature_set.non_pos_logFC),
                                         feature = feature_considered, feature_removed = is_feature_removed)
    }
    #create dataframe of current mean auc, number of features, current_feature_considered, is_feature_removed, iter number
    results <- rbind(results,
                     current_iter_results)
  }
  
  
  print(results)
  R.utils::mkdirs(results_file_path)
  res.file.name <- paste(rev(conditions), collapse="-vs-")
  write.csv(results, file.path(results_file_path, paste0(res.file.name, "-rfe-AUC_decrease-", meanAUC_decrease_cutoff,
                                                         "-pos_logFC_prop-",
                                                        pos_logFC_min_proportion, ".csv")), row.names = FALSE)
}


# ============ GBM V MET ========================
save.dir <- file.path("src", "DE_pipeline-2.8", "RFE", "rfe-results")
gbm.v.met.ranking <- read.csv("./src/DE_pipeline-2.8/compute-feature-ranks/feature-ranks/GBM-vs-MET-de.csv")
gbm.v.met.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-met-norm.csv",
                                row.names=1)
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)

# conduct.rfe(data=gbm.v.met.processed,
#             phenotype=gbm.v.met.labels,
#             conditions=c("MET", "GBM"),
#             rankings=gbm.v.met.ranking,
#             file.path(save.dir),
#             0.01)

conduct.rfe(data=gbm.v.met.processed,
            phenotype=gbm.v.met.labels,
            conditions=c("MET", "GBM"),
            rankings=gbm.v.met.ranking,
            file.path(save.dir),
            0.03)

conduct.rfe(data=gbm.v.met.processed,
            phenotype=gbm.v.met.labels,
            conditions=c("MET", "GBM"),
            rankings=gbm.v.met.ranking,
            file.path(save.dir),
            0.05)

# ================== HC v GBM

gbm.v.hc.ranking <- read.csv("./src/DE_pipeline-2.8/compute-feature-ranks/feature-ranks/GBM-vs-HC-de.csv")
gbm.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/gbm-v-hc-norm.csv",
                                row.names=1)
gbm.v.hc.labels <- read.table("src/process-2.6/combine-data-3/GBMVsHC_output_labels.txt", header=T)

conduct.rfe(data=gbm.v.hc.processed,
            phenotype=gbm.v.hc.labels,
            conditions=c("HC", "GBM"),
            rankings=gbm.v.hc.ranking,
            file.path(save.dir),
            0.01)

conduct.rfe(data=gbm.v.hc.processed,
            phenotype=gbm.v.hc.labels,
            conditions=c("HC", "GBM"),
            rankings=gbm.v.hc.ranking,
            file.path(save.dir),
            0.03)

conduct.rfe(data=gbm.v.hc.processed,
            phenotype=gbm.v.hc.labels,
            conditions=c("HC", "GBM"),
            rankings=gbm.v.hc.ranking,
            file.path(save.dir),
            0.05)

# ==================== MET V HC

met.v.hc.ranking <- read.csv("./src/DE_pipeline-2.8/compute-feature-ranks/feature-ranks/MET-vs-HC-de.csv")
met.v.hc.processed <- read.csv("src/process-2.6/combine-data-3/met-v-hc-norm.csv",
                               row.names=1)
met.v.hc.labels <- read.table("src/process-2.6/combine-data-3/METVsHC_output_labels.txt", header=T)

conduct.rfe(data=met.v.hc.processed,
            phenotype=met.v.hc.labels,
            conditions=c("HC", "MET"),
            rankings=met.v.hc.ranking,
            file.path(save.dir),
            0.01)

conduct.rfe(data=met.v.hc.processed,
            phenotype=met.v.hc.labels,
            conditions=c("HC", "MET"),
            rankings=met.v.hc.ranking,
            file.path(save.dir),
            0.03)

conduct.rfe(data=met.v.hc.processed,
            phenotype=met.v.hc.labels,
            conditions=c("HC", "MET"),
            rankings=met.v.hc.ranking,
            file.path(save.dir),
            0.05)





