# This file is for 10-fold cross validation
# without information leakage.
############################
# Load required libraries
############################
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(mRMRe))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(optparse))

source("src/nested/new-implementation-RNA/de.R")
source("src/nested/new-implementation-RNA/selections.R")
source("src/nested/new-implementation-RNA/helper.R")
source("src/nested/new-implementation-RNA/models.R")

# source("src/novel-discovery/helper.R")
# source("src/novel-discovery/obtain-data.R")
# source("src/novel-discovery/load-data.R")

option_list <- list( 
  make_option(c("-s", "--seed"), type="integer", default=42,
              help="Seed to initialze RNG with [default %default]"),
  # make_option(c("-M", "--max-iter"), type="integer", default=100,
  #             help="max iterations"),
  # make_option(c("-k", "--k"), type="integer", default=5,
  #             help="K-folds"),
  make_option(c("-c", "--comparison"), type="character", default="gbm.v.met",
              help="Which comparison to look at"),
  make_option(c("-P", "--sel-prop"), type="numeric", default=0.01,
              help="Selection proportion"),
  make_option(c("-k", "--top-k"), type="numeric", default=NULL,
              help="Select the top K? Can be used in conjunction with sel-prop"),
  make_option(c("-B", "--bootstraps"), type="numeric", default=50,
              help="# boots"),
  make_option(c("-S", "--selection"), type="character", default="lasso",
              help="The feature selection to conduct"),
  make_option(c("-m", "--model"), type="character", default="ridge",
              help="Model to use"),
  make_option(c("-f", "--fold-change"), type="numeric", default=1.5,
              help="The fold change to use.")
)

############################
# Parse arguments
############################

opt <- parse_args(OptionParser(option_list=option_list))

seed <- opt$seed
set.seed(seed)

# max.iter <- opt$`max-iter`

# k <- opt$k

sel.prop <- opt$`sel-prop`

boots <- opt$bootstraps

selection <- opt$selection

model.to.use <- opt$model

comparison <- opt$comparison

top.k <- opt$`top-k`

fold.change <- opt$`fold-change`

# folder.name <- paste0("model-", model.to.use, "_feature-sel-",
#                       selection,
#                       "_max-iter-", max.iter,
#                       "_k-", k, "_boots-", boots,
#                       "_sel-prop-", sel.prop, "_comparison-",
#                       comparison)
folder.name <- paste0("fc-", fold.change, "_model-", model.to.use, 
                      "_feature-sel-",
                      selection,
                      "_boots-", boots,
                      "_sel-prop-", sel.prop, "_sel-k-",
                      top.k, "_comparison-",
                      comparison, ".rds")

############################
# Load in the data
############################

fold.splits <- readRDS(file.path("./src/nested/preprocess-RNA/processed-fold-data/",
                                 comparison, "processed-fold-data.rds"))

# ================== 
# Get total # of cores
# ==================

get.n.cores <- function() {
  slurm <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  if (!is.na(slurm)) slurm else max(1L, detectCores() - 1L)
}
n.cores <- get.n.cores()

############################
# Initiate LSOCV
############################

# count - the counts for the patients
# sample_info - the corresponding meta data for each patient.
run.fold <- function(train.x, train.y, test.x, test.y, negative.class,
                     positive.class) {
  # Create a dataframe to save our results. It'll save the Sample ID,
  # parent code, true class label, prediction score.
  pred.df <- data.frame(matrix(data=NA, nrow=0, ncol=4))
  colnames(pred.df) <- c("sample.id", "parent", "true", "pred")
  
  acc.c <- c()
  spec.c <- c()
  auc.c <- c()
  sens.c <- c()
  run.i <- 0
  
  sel.features <- list()
  # print(dim(train.x))
  
  # Run DE
  de.res.i <- compute.de(train.x, train.y, c(negative.class, positive.class), fold.change)
  
  train.x <- train.x[,de.res.i,drop=F]
  test.x <- test.x[,de.res.i,drop=F]
  
  # Run feature selection
  sel.res <- decide.selection(selection, boots, train.x, train.y, n.cores, positive.class)
  sel.res <- as.data.frame(do.call(cbind, sel.res))
  # Maybe I should change below to obtain top 25 features that are at least above sel.prop
  sel.props <- sel.res %>% rowMeans
  # print(sum(sel.props >= sel.prop))
  if (is.null(top.k)) {
    sel.res.i <- which(sel.props >= sel.prop)
  } else {
    sel.res.i <- which(sel.props >= sel.prop & rank(-sel.props, ties.method="min") <= top.k)
  }
  print(sel.res.i)
  train.x <- train.x[,sel.res.i,drop=F]
  test.x <- test.x[,sel.res.i,drop=F]
  
  model.data <- decide.model(train.x, train.y, model.to.use, negative.class, positive.class)
  model <- model.data$model
  threshold <- model.data$threshold
  
  # Obtain metrics
  metric.list <- obtain.metrics(test.x, test.y, 
                                model, threshold, model.to.use,
                                negative.class, positive.class)
  print(test.y)
  print(metric.list)
  list(metric=metric.list,
       features=colnames(train.x))  
}

seeds <- sample(.Machine$integer.max, length(fold.splits))

# Get the location of the directory to save in/load from
out.dir <- file.path('results', "nested", paste0(comparison, ".RNA"),
                     paste0("model-", model.to.use, "_fs-", selection))
R.utils::mkdirs(out.dir)

auc.c <- c()
acc.c <- c()
sens.c <- c()
spec.c <- c()
features.sel <- list()

for (i in 1:length(fold.splits)) {
  ret.list <- run.fold(fold.splits[[i]]$train.x, factorise(fold.splits[[i]]$train.y, comparison),
           fold.splits[[i]]$test.x, factorise(fold.splits[[i]]$test.y, comparison),
           get.negative.class(comparison), get.positive.class(comparison))
  metric.list <- ret.list$metric
  auc.c <- c(auc.c, metric.list$auc)
  acc.c <- c(acc.c, metric.list$acc)
  sens.c <- c(sens.c, metric.list$sens)
  spec.c <- c(spec.c, metric.list$spec)
  features.sel[[i]] <- ret.list$features
}

save.list <- list(features=features.sel,
                  auc=auc.c,
                  acc=acc.c,
                  sens=sens.c,
                  spec=spec.c)

saveRDS(save.list, file.path(out.dir, folder.name))
