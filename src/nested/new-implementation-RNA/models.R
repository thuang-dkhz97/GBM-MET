suppressPackageStartupMessages(library(caret))

# Train the model
decide.model <- function(x, y, model.name, negative.class, positive.class) {
  if (model.name == "ridge") {
    fold.id <- drop(caret::createFolds(y, k = 10, list = FALSE))
    x <- as.matrix(x)
    model <- cv.glmnet(x,
                       y,
                       alpha=0,
                       family=binomial,
                       foldid=fold.id)
    train.preds <- c(predict(model, x, type="response", s="lambda.1se"))
  } else if (model.name == "rf") {
    tc <- trainControl(method='LGOCV',
                       number=100,
                       p=0.9,
                       allowParallel = F,
                       search='random',
                       savePredictions='final',
                       classProbs=T,
                       summaryFunction=twoClassSummary)
    model <- train(x=x,
                   y=y,
                   trControl=tc,
                   method="rf",
                   metric="ROC")
    train.preds <- predict(model, x, type='prob')[,2]
  }
  roc_obj <- invisible(roc(y, train.preds, levels=c(negative.class, positive.class)))
  best_coords <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "specificity", "sensitivity"))
  
  list("model"=model,
       "threshold"=best_coords$threshold[1])
}

# Obtain the metrics on x and y
obtain.metrics <- function(x, y, model, threshold, model.name, negative.class, positive.class) {
  if (model.name == "ridge") {
    x <- as.matrix(x)
    val.preds <- c(predict(model, x, type="response", s="lambda.1se"))
  } else if (model.name == "rf") {
    val.preds <- c(predict(model, x, type='prob')[,2])
  }
  
  # Calculate auc, etc.
  pred.class <- ifelse(val.preds >= threshold, positive.class, negative.class)
  acc.val <- MLmetrics::Accuracy(y_true =y, y_pred=pred.class)
  sens.val <- MLmetrics::Sensitivity(y_true =y, y_pred=pred.class, positive=positive.class)
  if (is.na(sens.val) & positive.class %in% y) {
    sens.val <- 0
  }
  spec.val <- MLmetrics::Specificity(y_true =y, y_pred=pred.class, positive=positive.class)
  if (is.na(spec.val) & negative.class %in% y) {
    spec.val <- 0
  }
  auc.val <- MLmetrics::AUC(y_pred=val.preds, y_true=ifelse(y == positive.class, 1, 0))
  list("auc"=auc.val,
       "sens"=sens.val,
       "spec"=spec.val,
       "acc"=acc.val)
}

