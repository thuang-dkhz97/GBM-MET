suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))


hcbb.filter <- read.csv("./src/process-2.6/combine-data-1/hcbb-combine-data-1.csv",
                        row.names=1, check.names=F)
sbtb.filter <- read.csv("./src/process-2.6/combine-data-1/sbtb-combine-data-1.csv",
                        row.names=1, check.names=F)

stopifnot(identical(rownames(hcbb.filter), rownames(sbtb.filter)))
stopifnot(!anyDuplicated(colnames(combine.filter)))

# We will conduct imputation in the nested fold later
# Will be followed by quantile normalisation
# And then BE correction --> will need to know which samples came from
# sbtb and hcbb.

# First, combine the two datasets.
combine.filter <- cbind(hcbb.filter, sbtb.filter)

# ============== Separate data based on comparisons

sample.data <- read.xlsx("./data/class-labels.xlsx")
obtain.class.data <- function(data, classes) {
  sample.ids <- sample.data$Sample.ID[sample.data$Class %in% classes]
  keep <- colnames(data)[colnames(data) %in% sample.ids]
  list("data"=data[,keep, drop=F],
       "class"=sample.data$Class[match(keep, sample.data$Sample.ID)])
}

gbm.v.met.data <- obtain.class.data(combine.filter, c("GBM", "MET"))
gbm.v.met <- gbm.v.met.data$data
gbm.v.hc.data <- obtain.class.data(combine.filter, c("GBM", "HC"))
gbm.v.hc <- gbm.v.hc.data$data
met.v.hc.data <- obtain.class.data(combine.filter, c("MET", "HC"))
met.v.hc <- met.v.hc.data$data


create.fold.split <- function(data, phenotype, classification_criteria, k=5, sample.total=30) {
  # Function to create the stratified folds.
  set.seed(1000)
  filter <- expression(TRUE)
  extracted_samples <- phenotype %>% subset(eval(filter))
  extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
  extracted_samples$Sample <- factor(extracted_samples$Sample)
  filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
  colnames(filtered_samples_output_labels) <- c("Sample", "Label")
  
  # Sort filtered_samples_output_labels based on data
  filtered_samples_output_labels <- filtered_samples_output_labels[match(colnames(data), filtered_samples_output_labels$Sample),]
  
  # filtered_samples_output_labels should be aligned to the columns of data.
  iter <- sample.total/k
  train.fold.i <- createMultiFolds(y=filtered_samples_output_labels$Label, k=k, times=iter)
  train.fold.i
}

# Create fold splits.
labels.table <- read.table("data/proteomic-phenotype.txt", header=TRUE, sep="\t")
gbm.v.met.fold.i <- create.fold.split(gbm.v.met, labels.table, "GBMVsMET")
gbm.v.hc.fold.i <- create.fold.split(gbm.v.hc, labels.table, "GBMVsHC")
met.v.hc.fold.i <- create.fold.split(met.v.hc, labels.table, "METVsHC")

# Save
save.data <- list("gbm.v.met"=list("fold.i"=gbm.v.met.fold.i,
                                   "data"=gbm.v.met),
                  "gbm.v.hc"=list("fold.i"=gbm.v.hc.fold.i,
                                  "data"=gbm.v.hc),
                  "met.v.hc"=list("fold.i"=met.v.hc.fold.i,
                                  "data"=met.v.hc))

R.utils::mkdirs("./src/nested/fold-data")
saveRDS(save.data, "./src/nested/fold-data/fold-data.rds")

R.utils::mkdirs("./src/nested/combine-data-2")
write.csv(gbm.v.met, "./src/nested/combine-data-2/gbm-v-met.csv")
write.csv(gbm.v.hc, "./src/nested/combine-data-2/gbm-v-hc.csv")
write.csv(met.v.hc, "./src/nested/combine-data-2/met-v-hc.csv")



