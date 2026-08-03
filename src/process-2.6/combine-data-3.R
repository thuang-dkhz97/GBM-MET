library(dplyr)
library(sva)
library(openxlsx)

# Then generate 3 separate data files for comparing
# GBM vs MET, GBM vs HC, MET vs HC. 
# Conduct quantile norm and if required, BE correction
# COMBAT is applied only to the GBM class.

# ================== Load data =======================

hcbb.filter <- read.csv("./src/process-2.6/combine-data-2/hcbb-combine-data-2.csv",
                        row.names=1, check.names=F)
sbtb.filter <- read.csv("./src/process-2.6/combine-data-2/sbtb-combine-data-2.csv",
                        row.names=1, check.names=F)
combine.filter <- cbind(hcbb.filter, sbtb.filter)

# ============== Separate data based on comparisons

sample.data <- read.xlsx("./data/class-labels.xlsx")
obtain.class.data <- function(data, classes) {
  sample.ids <- sample.data$Sample.ID[sample.data$Class %in% classes]
  list("data"=data[,colnames(data) %in% sample.ids],
       "class"=sample.data$Class[which(sample.data$Sample.ID %in% sample.ids)])
}

gbm.v.met.data <- obtain.class.data(combine.filter, c("GBM", "MET"))
gbm.v.met <- gbm.v.met.data$data
gbm.v.hc.data <- obtain.class.data(combine.filter, c("GBM", "HC"))
gbm.v.hc <- gbm.v.hc.data$data
met.v.hc.data <- obtain.class.data(combine.filter, c("MET", "HC"))
met.v.hc <- met.v.hc.data$data

# ================ Quantile norm

perform.quantile.norm <- function(data) {
  data.train.rank <- apply(data, 2, rank, ties.method="average")
  data.train.sorted <- data.frame(apply(data, 2, sort))
  data.train.mean <- apply(data.train.sorted, 1, mean)
  
  index_to_mean <- function(index, data_mean){
    #index can be int or int+0.5
    #if int+0.5, take average of the numbers in those positions
    int.result <- data_mean[index]
    index.int <- floor(index)
    #some of the values in point5.result might be NA
    #but they won't be chosen
    point5.result <- (data_mean[index.int] + data_mean[index.int+1])/2
    point5.indices <- index%%1 != 0
    result <- int.result
    result[point5.indices] <- point5.result[point5.indices]
    return (result)
  }
  data.train.norm <- apply(data.train.rank, 2, index_to_mean, data_mean = data.train.mean)
  rownames(data.train.norm) <- rownames(data)
  data.train.norm
}

gbm.v.met.norm <- perform.quantile.norm(gbm.v.met)
gbm.v.hc.norm <- perform.quantile.norm(gbm.v.hc)
met.v.hc.norm <- perform.quantile.norm(met.v.hc)

# ================= COMBAT (for gbm.v.met.norm and gbm.v.hc.norm)

apply.combat <- function(v.data) {
  # Need to apply COMBAT to the gbm cohort.
  # Extract the gbm cohorts
  gbm.cohorts.i <- which(sample.data$Class == "GBM")
  gbm.cohorts.meta <- sample.data[gbm.cohorts.i,]
  
  # Now, extract the gbm cohorts from our data
  cross.i <- match(gbm.cohorts.meta$Sample.ID, colnames(v.data))
  gbm.cohorts.data <- v.data[,cross.i]
  # Create the batch list from this
  gbm.batch.label <- rep(1, ncol(gbm.cohorts.data))
  # Set the second batch to cohort B.
  gbm.batch.label[gbm.cohorts.meta$Cohort == "B"] <- 2
  # Apply combat
  combat.data <- ComBat(gbm.cohorts.data, gbm.batch.label, ref.batch=1)
  v.data[,cross.i] <- combat.data
  v.data
}

gbm.v.met.norm <- apply.combat(gbm.v.met.norm)
gbm.v.hc.norm <- apply.combat(gbm.v.hc.norm)

# ============= Save

if (!dir.exists("./src/process-2.6/combine-data-3")) {
  mkdirs("./src/process-2.6/combine-data-3")
}

write.csv(gbm.v.met.norm, "./src/process-2.6/combine-data-3/gbm-v-met-norm.csv")
write.csv(gbm.v.hc.norm, "./src/process-2.6/combine-data-3/gbm-v-hc-norm.csv")
write.csv(met.v.hc.norm, "./src/process-2.6/combine-data-3/met-v-hc-norm.csv")
