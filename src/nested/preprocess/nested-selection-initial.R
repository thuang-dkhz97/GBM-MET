
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(missForest))
suppressPackageStartupMessages(library(missForestPredict))
suppressPackageStartupMessages(library(openxlsx))
library(factoextra)

fold.splits <- readRDS("./src/nested/preprocess/fold-data/fold-data.rds")

# Use the data below to determine which samples are from
# hcbb and sbtb.
hcbb.filter <- read.csv("./src/process-2.6/combine-data-1/hcbb-combine-data-1.csv",
                        row.names=1, check.names=F)
hcbb.samples <- colnames(hcbb.filter)
sbtb.filter <- read.csv("./src/process-2.6/combine-data-1/sbtb-combine-data-1.csv",
                        row.names=1, check.names=F)
sbtb.samples <- colnames(sbtb.filter)

processed.fold.data <- list()


process.data.na <- function(train.sbtb.x, test.sbtb.x, train.hcbb.x, test.hcbb.x) {
  # For each training set, check which features are completely missing.
  # This may happen due to random split in data.
  sbtb.missing.feat.i <- which(colSums(is.na(train.sbtb.x)) == nrow(train.sbtb.x) |
                                 colSums(is.na(train.sbtb.x)) == nrow(train.sbtb.x) - 1)
  hcbb.missing.feat.i <- which(colSums(is.na(train.hcbb.x)) == nrow(train.hcbb.x) |
                                 colSums(is.na(train.hcbb.x)) == nrow(train.hcbb.x) - 1)
  # Set all of these features to 0.
  train.sbtb.x[,union(sbtb.missing.feat.i, hcbb.missing.feat.i)] <- 0
  test.sbtb.x[,union(sbtb.missing.feat.i, hcbb.missing.feat.i)] <- 0
  train.hcbb.x[,union(sbtb.missing.feat.i, hcbb.missing.feat.i)] <- 0
  test.hcbb.x[,union(sbtb.missing.feat.i, hcbb.missing.feat.i)] <- 0
  list(train.sbtb.x, test.sbtb.x, train.hcbb.x, test.hcbb.x)
}

apply.missforest <- function(train.x, test.x) {
  miss.model <- missForest(as.data.frame(train.x), verbose=F)
  train.x <- miss.model$ximp
  # Then compute missforest for the corresponding test sets.
  test.x <- missForestPredict(miss.model, newdata=as.data.frame(test.x))
  list(train.x, test.x)
}

apply.quantile <- function(train.x, test.x) {
  train.m <- t(as.matrix(train.x))             
  test.m  <- t(as.matrix(test.x))
  target  <- preprocessCore::normalize.quantiles.determine.target(train.m)
  train.n <- preprocessCore::normalize.quantiles.use.target(train.m, target)
  test.n  <- preprocessCore::normalize.quantiles.use.target(test.m,  target)
  dimnames(train.n) <- dimnames(train.m)
  dimnames(test.n)  <- dimnames(test.m)
  list(as.data.frame(t(train.n)), as.data.frame(t(test.n)))
  # norm_data <- preprocessCore::normalize.quantiles(as.matrix(t(train.x)))
  # norm_data <- data.frame(norm_data, row.names = colnames(train.x))
  # colnames(norm_data) <- rownames(train.x)
  # train.x <- norm_data
  # 
  # norm_data <- preprocessCore::normalize.quantiles(as.matrix(t(test.x)))
  # norm_data <- data.frame(norm_data, row.names = colnames(test.x))
  # colnames(norm_data) <- rownames(test.x)
  # test.x <- norm_data
  # 
  # train.x <- as.data.frame(t(as.matrix(train.x)))
  # test.x <- as.data.frame(t(as.matrix(test.x)))
  # list(train.x, test.x)
}

perform.combat.correct <- function(x.train, x.test) {
  x.train <- as.data.frame(as.matrix(x.train))
  train.samples <- rownames(x.train)
  x.test <- as.data.frame(as.matrix(x.test))
  test.samples <- rownames(x.test)
  
  sample.data <- read.xlsx("./data/class-labels.xlsx")
  # Extract the gbm cohorts
  gbm.cohorts.i <- which(sample.data$Class == "GBM")
  gbm.cohorts.meta <- sample.data[gbm.cohorts.i,]
  
  apply.combat.gbm <- function(v.data) {
    # Need to apply COMBAT to the gbm cohort.
    # Now, extract the gbm cohorts from our data
    cross.i <- match(gbm.cohorts.meta$Sample.ID, rownames(v.data))
    
    # Here, I need a which, since I want to extract the positions
    # which don't have na's to get the matching index.
    gbm.cohorts.data <- v.data[cross.i[which(!is.na(cross.i))],]
    # Create the batch list from this
    gbm.batch.label <- rep(1, nrow(gbm.cohorts.data))
    # Set the second batch to cohort B.
    # We don't have a which below, since I only want to extract those
    # rows with matches
    gbm.batch.label[gbm.cohorts.meta$Cohort[!is.na(cross.i)] == "B"] <- 2
    # Apply combat
    combat.data <- suppressMessages(ComBat(t(gbm.cohorts.data), gbm.batch.label, ref.batch=1))
    v.data[cross.i[which(!is.na(cross.i))],] <- t(combat.data)
    v.data
  }
  
  # Below is to apply combat across a train and test set.
  apply.combat <- function(x.train, x.test) {
    combine.set <- rbind(x.train, x.test)
    batch.label <- rep(1, nrow(x.train) + nrow(x.test))
    batch.label[(nrow(x.train)+1):(nrow(x.train) + nrow(x.test))] <- 2
    combat.data <- suppressMessages(t(ComBat(t(combine.set), batch.label, ref.batch=1)))
    
    combat.data[(nrow(x.train)+1):(nrow(x.train) + nrow(x.test)),]
  }
  
  # First, correct GBM within the training set
  train.cross.i <- match(gbm.cohorts.meta$Sample.ID, train.samples)
  if (any(gbm.cohorts.meta$Cohort[!is.na(train.cross.i)] == "B")) {
    x.train <- apply.combat.gbm(x.train)
  }
  
  # Then correct the full test set against the train set.
  # x.test <- apply.combat(x.train, x.test)
  
  # =============== Below has been commented out, the
  # =============== reviewer wants realistic performance.
  # Then correct GBM Cohort 2 using the training set GBMs
  test.samples <- rownames(x.test)
  
  cross.i <- match(gbm.cohorts.meta$Sample.ID, test.samples)
  # Extract samples which are cohort 2.
  x.test.gbm.cohort.data <- x.test[cross.i[which(!is.na(cross.i))],]
  # Extract those in cohort B.
  if (any(gbm.cohorts.meta$Cohort[!is.na(cross.i)] == "B")) {
    # Only correct the cohort B's against the train set gbm cohorts, if any.
    x.test.gbm.cohort.b.data <- x.test.gbm.cohort.data[gbm.cohorts.meta$Cohort[!is.na(cross.i)] == "B",]
    # Extract the training set GBM
    # ================= TO DO =================
    train.cross.i <- match(gbm.cohorts.meta$Sample.ID, train.samples)
    x.train.gbm <- x.train[train.cross.i[which(!is.na(train.cross.i))],,drop=F]
    
    # Obtain the index of the gbm and cohort B's.
    x.test.gbm.cohort.data[gbm.cohorts.meta$Cohort[!is.na(cross.i)] == "B",] <- apply.combat(x.train.gbm, x.test.gbm.cohort.b.data)
    x.test[cross.i[which(!is.na(cross.i))], ] <- x.test.gbm.cohort.data
  }
  
  list(x.train, x.test)
}


parse.folds <- function(fold.splits, comparison="gbm.v.met") {
  R.utils::mkdirs(file.path("./src/nested/preprocess/processed-fold-data", comparison))
  if (comparison != "met.v.hc") {
    for (i in 1:length(fold.splits[[comparison]]$fold.i)) {
      # Perform missforest for sbtb and hcbb separately
      # Then perform quantile norm
      # Then apply BE correction on the GBM classes only.
      train.fold.i <- fold.splits[[comparison]]$fold.i[[i]]
      
      # Extract hcbb and sbtb data
      train.fold.data <- t(fold.splits[[comparison]]$data[,train.fold.i,drop=F])
      train.fold.sbtb.data <- train.fold.data[rownames(train.fold.data) %in% sbtb.samples,,drop=F]
      train.fold.hcbb.data <- train.fold.data[rownames(train.fold.data) %in% hcbb.samples,,drop=F]
      
      test.fold.data <- t(fold.splits[[comparison]]$data[,-train.fold.i,drop=F])
      test.fold.sbtb.data <- test.fold.data[rownames(test.fold.data) %in% sbtb.samples,,drop=F]
      test.fold.hcbb.data <- test.fold.data[rownames(test.fold.data) %in% hcbb.samples,,drop=F]
      
      processed.data <- process.data.na(train.fold.sbtb.data, test.fold.sbtb.data,
                                        train.fold.hcbb.data, test.fold.hcbb.data)
      train.fold.sbtb.data <- processed.data[[1]]
      test.fold.sbtb.data <- processed.data[[2]]
      train.fold.hcbb.data <- processed.data[[3]]
      test.fold.hcbb.data <- processed.data[[4]]

      # Combine train data
      sbtb.data <- apply.missforest(train.fold.sbtb.data, test.fold.sbtb.data)
      train.fold.sbtb.data <- sbtb.data[[1]]
      test.fold.sbtb.data <- sbtb.data[[2]]
      hcbb.data <- apply.missforest(train.fold.hcbb.data, test.fold.hcbb.data)
      train.fold.hcbb.data <- hcbb.data[[1]]
      test.fold.hcbb.data <- hcbb.data[[2]]

      train.fold.data <- rbind(train.fold.sbtb.data, train.fold.hcbb.data)
      # Combine test data
      test.fold.data <- rbind(test.fold.sbtb.data, test.fold.hcbb.data)

      # ========== MissForest complete
      # ========= Quantile norm
      quantile.norm <- apply.quantile(train.fold.data, test.fold.data)

      # Combine train data
      train.fold.data <- quantile.norm[[1]]
      test.fold.data <- quantile.norm[[2]]
      # =====================

      # ========= ComBat ==================
      combat.data <- perform.combat.correct(train.fold.data, test.fold.data)
      train.fold.data <- combat.data[[1]]
      test.fold.data <- combat.data[[2]]

      sample.data <- read.xlsx("./data/class-labels.xlsx")
      train.match.samples.i <- match(rownames(train.fold.data), sample.data$Sample.ID)
      test.match.samples.i <- match(rownames(test.fold.data), sample.data$Sample.ID)
      # Save the processed data.
      train.y <- sample.data$Class[train.match.samples.i]
      test.y <- sample.data$Class[test.match.samples.i]
      processed.fold.data[[i]] <- list(train.x=train.fold.data, train.y=train.y,
                                       test.x=test.fold.data, test.y=test.y)

      R.utils::mkdirs(file.path("./src/nested/preprocess/processed-fold-data", comparison, i))
      write.csv(t(train.fold.data), file.path("./src/nested/preprocess/processed-fold-data", comparison, i, "train.csv"))
      write.csv(t(test.fold.data), file.path("./src/nested/preprocess/processed-fold-data", comparison, i, "test.csv"))
      
      # load the proteomic-phenotype-txt and only store
      # for the samples in this fold
      proteomic.phenotype <- read.table("data/proteomic-phenotype.txt", header=T, sep="\t")
      proteomic.phenotype <- proteomic.phenotype[proteomic.phenotype$Sample.ID %in% rownames(train.fold.data),]
      
      print(nrow(proteomic.phenotype))
      
      write.table(proteomic.phenotype, file.path("./src/nested/preprocess/processed-fold-data", comparison, i, "proteomic-phenotype.txt"),
                  sep="\t", row.names=FALSE)
      
    }
  } else {
    for (i in 1:length(fold.splits[[comparison]]$fold.i)) {
      train.fold.i <- fold.splits[[comparison]]$fold.i[[i]]
      
      # Extract hcbb and sbtb data
      train.fold.data <- t(fold.splits[[comparison]]$data[,train.fold.i,drop=F])
      test.fold.data <- t(fold.splits[[comparison]]$data[,-train.fold.i,drop=F])
      
      train.missing.feat.i <- which(colSums(is.na(train.fold.data)) == nrow(train.fold.data) |
                                      colSums(is.na(train.fold.data)) == nrow(train.fold.data) - 1)
      # Set all of these features to NULL.
      train.fold.data[,train.missing.feat.i] <- 0
      test.fold.data[,train.missing.feat.i] <- 0
      
      imputed.data <- apply.missforest(train.fold.data, test.fold.data)
      train.fold.data <- imputed.data[[1]]
      test.fold.data <- imputed.data[[2]]
      
      quantile.data <- apply.quantile(train.fold.data, test.fold.data)
      train.fold.data <- quantile.data[[1]]
      test.fold.data <- quantile.data[[2]]
      
      sample.data <- read.xlsx("./data/class-labels.xlsx")
      train.match.samples.i <- match(rownames(train.fold.data), sample.data$Sample.ID)
      test.match.samples.i <- match(rownames(test.fold.data), sample.data$Sample.ID)
      # Save the processed data.
      train.y <- sample.data$Class[train.match.samples.i]
      test.y <- sample.data$Class[test.match.samples.i]
      processed.fold.data[[i]] <- list(train.x=train.fold.data, train.y=train.y,
                                       test.x=test.fold.data, test.y=test.y)
      
      R.utils::mkdirs(file.path("./src/nested/preprocess/processed-fold-data", comparison, i))
      write.csv(t(train.fold.data), file.path("./src/nested/preprocess/processed-fold-data", comparison, i, "train.csv"))
      write.csv(t(test.fold.data), file.path("./src/nested/preprocess/processed-fold-data", comparison, i, "test.csv"))
    
      # load the proteomic-phenotype-txt and only store
      # for the samples in this fold
      proteomic.phenotype <- read.table("data/proteomic-phenotype.txt", header=T, sep="\t")
      proteomic.phenotype <- proteomic.phenotype[proteomic.phenotype$Sample.ID %in% rownames(train.fold.data),]
      
      write.table(proteomic.phenotype, file.path("./src/nested/preprocess/processed-fold-data", comparison, i, "proteomic-phenotype.txt"),
                  sep="\t", row.names=FALSE)  
    }
  }
  saveRDS(processed.fold.data, file.path("./src/nested/preprocess/processed-fold-data", comparison, "processed-fold-data.rds"))
}

parse.folds(fold.splits, "gbm.v.met")
parse.folds(fold.splits, "gbm.v.hc")
parse.folds(fold.splits, "met.v.hc")



