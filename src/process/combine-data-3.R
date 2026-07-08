library(dplyr)
library(openxlsx)

# Then generate 3 separate data files for comparing
# GBM vs MET, GBM vs HC, MET vs HC. 
# Conduct quantile norm and if required, BE correction

# ================== Load data =======================

hcbb.filter <- read.csv("./src/process/combine-data-2/hcbb-combine-data-2.csv",
                        row.names=1, check.names=F)
sbtb.filter <- read.csv("./src/process/combine-data-2/sbtb-combine-data-2.csv",
                        row.names=1, check.names=F)
combine.filter <- cbind(hcbb.filter, sbtb.filter)

# ============== Separate data based on comparisons

sample.data <- read.xlsx("./data/class-labels.xlsx")
obtain.class.data <- function(data, classes) {
  sample.ids <- sample.data$Sample.ID[sample.data$Class %in% classes]
  data[,colnames(data) %in% sample.ids]
}

gbm.v.met <- obtain.class.data(combine.filter, c("GBM", "MET"))
gbm.v.hc <- obtain.class.data(combine.filter, c("GBM", "HC"))
met.v.hc <- obtain.class.data(combine.filter, c("MET", "HC"))

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

# ================= COMBAT (for gbm.v.met.norm)

# UP TO HERE


# ============= Save

if (!dir.exists("./src/process/combine-data-3")) {
  mkdirs("./src/process/combine-data-3")
}

# write.csv(hcbb.filter, "./src/process/combine-data-3/hcbb-combine-data-3.csv")
# write.csv(sbtb.filter, "./src/process/combine-data-3/sbtb-combine-data-3.csv")
