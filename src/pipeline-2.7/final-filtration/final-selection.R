library(stringr)
# In this file, loop through the data and get the selection of the features
# that have been selected in this ML pipeline.

decide.ml.panel <- function(source.folder = "results/second-run/GBM_VS_MET") {
  
  internal.folders <- list.dirs(path=source.folder)[-c(1)]
  
  auc.data <- data.frame(matrix(data=NA, ncol=5, nrow=0))
  
  # For each folder, load the meanAUC's obtained by each panel.
  for (folder in internal.folders) {
    # Open folder with metrics
    metric.file <- file.path(folder, "model_results_test.csv")
    metric.data <- read.csv(metric.file)
    # Extract mean AUC.
    save.met.data <- metric.data[, c("FSM", "Model", "Mean_AUC")]
    # Extract FSM.
    fsm.method <- str_extract(folder, paste0("(?<=", "fsm-", ").*?(?=", "-feat", ")"))
    save.met.data$FSM <- fsm.method
    # Extract the threshold selection size of the FSM.
    fsm.panel.size <- str_sub(folder, -2, -1)
    save.met.data$sel.size <- fsm.panel.size
    
    # Extract number of features used.
    inner.folder.list <- list.files(folder, full.names=T)
    sel.data <- read.csv(inner.folder.list[which(grepl("_features", inner.folder.list))])
    num.features <- length(colnames(sel.data)) - 2
    save.met.data$num.features <- num.features
    
    auc.data <- rbind(auc.data, save.met.data)
  }
  
  auc.data

  # Identify the best performing models and extract the number of features
  # used.
  max.auc.data <- auc.data[which(auc.data$Mean_AUC == max(auc.data$Mean_AUC)),]
  # Extract the fsm with the smallest number of features.
  max.auc.data.min.feature <- max.auc.data[which(max.auc.data$num.features == min(max.auc.data$num.features)),]
  
  max.auc.data.min.feature
  
}

# ========== GBM V MET
data <- decide.ml.panel()
print(data)
d <- data[1,]
# Take the first model, and record/save which features led to the highest performance
counts <- read.csv(file.path("src/pipeline-2.7/filtration-after-initial-selection/candidate",
                             paste0("gbm-vs-met-fsm-", d$FSM, "-feat-num-",
                                    d$sel.size, "-candidate.csv")),
                   row.names=1)
feature.sel <- rownames(counts)

gbm.v.met.save <- list("method"=d$FSM,
                       "sel.size"=d$sel.size,
                       "features"=feature.sel)
print(gbm.v.met.save)
# Save
save.dir <- file.path("src/pipeline-2.7/final-filtration/selected")
R.utils::mkdirs(save.dir)

saveRDS(gbm.v.met.save, file.path(save.dir, "gbm-vs-met-sels.rds"))

# ========== GBM V HC
data <- decide.ml.panel("results/second-run/GBM_VS_HC")
print(data)
d <- data[1,]
# Take the first model, and record/save which features led to the highest performance
counts <- read.csv(file.path("src/pipeline-2.7/filtration-after-initial-selection/candidate",
                             paste0("gbm-vs-hc-fsm-", d$FSM, "-feat-num-",
                                    d$sel.size, "-candidate.csv")),
                   row.names=1)
feature.sel <- rownames(counts)

gbm.v.hc.save <- list("method"=d$FSM,
                      "sel.size"=d$sel.size,
                      "features"=feature.sel)
print(gbm.v.hc.save)
# Save
save.dir <- file.path("src/pipeline-2.7/final-filtration/selected")
R.utils::mkdirs(save.dir)

saveRDS(gbm.v.hc.save, file.path(save.dir, "gbm-vs-hc-sels.rds"))

# ========== MET V HC
data <- decide.ml.panel("results/second-run/MET_VS_HC")
print(data)
d <- data[1,]
# Take the first model, and record/save which features led to the highest performance
counts <- read.csv(file.path("src/pipeline-2.7/filtration-after-initial-selection/candidate",
                             paste0("met-vs-hc-fsm-", d$FSM, "-feat-num-",
                                    d$sel.size, "-candidate.csv")),
                   row.names=1)
feature.sel <- rownames(counts)

met.v.hc.save <- list("method"=d$FSM,
                      "sel.size"=d$sel.size,
                      "features"=feature.sel)
print(met.v.hc.save)
# Save
save.dir <- file.path("src/pipeline-2.7/final-filtration/selected")
R.utils::mkdirs(save.dir)

saveRDS(met.v.hc.save, file.path(save.dir, "met-vs-hc-sels.rds"))

