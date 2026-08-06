# label.train <- output_labels.train
# label.test <- output_labels.test
# random_seed <- 1000

rf_model <- function(data.train, label.train, data.test, label.test, 
                     classes, random_seed = 1000, classifier_feature_imp = FALSE, ...){
  
  model_name <- "Random Forest"
  
  feature_imp <- data.frame(matrix(nrow = 0, ncol = 2,
                                   dimnames = list(NULL, c("feature", "meanDecreaseGini"))))
  
  try({
    set.seed(random_seed)
    data.sanity_check <- data.train - colMeans(data.train)
    if(sum(data.sanity_check != 0) != 0){
      #ensure that atleast one column is not a constant vector
      #does this by checking if there is atleast one non-zero entry after subtracting column mean
      
      #all columns constant causes the below line to run forever
      model <- randomForest::randomForest(x = data.train, y = factor(label.train, levels = classes))
      
      if(classifier_feature_imp){
        feature_imp <- data.frame(randomForest::importance(model))
        feature_imp <- cbind(feature = rownames(feature_imp),
                             data.frame(feature_imp, row.names = NULL))
      }
      
      best_acc <- -1
      best_cut_off <- 0.5
      for(cut_off in c(0.5, seq(0.2, 0.8, 0.01))){
        # print(cut_off)
        pred_prob.train <- predict(model, data.train, type = "prob")
        pred_prob.train <- data.frame(pred_prob.train)[classes[2]]
        pred.train <- ifelse(pred_prob.train > best_cut_off, classes[2], classes[1])
        acc <- mean(pred.train == label.train)
        if(acc > best_acc){
          best_acc <- acc
          best_cut_off <- cut_off
        }
      }
      
      
      pred_prob.train <- predict(model, data.train, type = "prob")
      pred_prob.train <- data.frame(pred_prob.train)[classes[2]]
      pred.train <- ifelse(pred_prob.train > best_cut_off, classes[2], classes[1])
      
      
      pred_prob <- predict(model, data.test, type="prob")
      pred_prob <- data.frame(pred_prob)[classes[2]]
      pred <- ifelse(pred_prob > best_cut_off, classes[2], classes[1])
      
      result_df.train <- data.frame("TrueLabel" = label.train,
                                    "Pred_prob" = pred_prob.train[,1],
                                    "PredictedLabel" = pred.train[,1],
                                    "Type" = "train",
                                    "cutoff" = best_cut_off)
      
      result_df.test <- data.frame("TrueLabel" = label.test,
                                   "Pred_prob" = pred_prob[,1],
                                   "PredictedLabel" = pred[,1],
                                   "Type" = "test",
                                   "cutoff" = best_cut_off)
      
      result_df <- rbind(result_df.train %>%
                           rownames_to_column("sample"), 
                         result_df.test %>%
                           rownames_to_column("sample"))
      
      result_df <- result_df %>%
        mutate(Pred_prob = as.double(Pred_prob))
      
    } else{
      print("data to RF : all fields constant!")
      print(dim(data.train))
      
      result_df <- data.frame("sample" = NA, 
                              "TrueLabel" = NA,
                              "Pred_prob" = NA,
                              "PredictedLabel" = NA,
                              "Type" = NA,
                              "cutoff" = NA)
    }
    
  })

  result_df <- result_df %>%
    mutate(model = model_name)
  
  if(classifier_feature_imp){
    return (list(result_df, feature_imp))   
  } else{
    return (result_df)     
  }
  
}
