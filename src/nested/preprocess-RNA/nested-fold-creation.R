suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tibble))

hcbb.data <- read.xlsx("./data/umi_counts_initial_cohort.xlsx")
hcbb.data <- column_to_rownames(hcbb.data, var="X1")
sbtb.data <- read.xlsx("./data/umi_counts_validation_cohort.xlsx")
sbtb.data <- column_to_rownames(sbtb.data, var="X1")
# Change the names in data.2, add an "S"
colnames(sbtb.data) <- paste0("S", colnames(sbtb.data))

sample.data <- read.xlsx("./data/class-labels.xlsx")

# Obtain rows that are in both.
common.genes <- intersect(rownames(hcbb.data), rownames(sbtb.data))

# Extract the filtered data
hcbb.filter <- hcbb.data[common.genes,colnames(hcbb.data) %in% sample.data$Sample.ID]
sbtb.filter <- sbtb.data[common.genes,colnames(sbtb.data) %in% sample.data$Sample.ID]

stopifnot(identical(rownames(hcbb.filter), rownames(sbtb.filter)))
stopifnot(!anyDuplicated(colnames(combine.filter)))

# We will conduct imputation in the nested fold later
# Will be followed by quantile normalisationZ
# And then BE correction --> will need to know which samples came from
# sbtb and hcbb.

# First, combine the two datasets.
combine.filter <- cbind(hcbb.filter, sbtb.filter)
# nrow(combine.filter)
# ncol(combine.filter)
# setdiff(sample.data$Sample.ID, colnames(combine.filter))
stopifnot(!anyDuplicated(colnames(combine.filter)))

# ============== Separate data based on comparisons

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

R.utils::mkdirs("./src/nested/preprocess-RNA/fold-data")
saveRDS(save.data, "./src/nested/preprocess-RNA/fold-data/fold-data.rds")

# R.utils::mkdirs("./src/nested/preprocess-RNA/combine-data-2")
# write.csv(gbm.v.met, "./src/nested/preprocess-RNA/combine-data-2/gbm-v-met.csv")
# write.csv(gbm.v.hc, "./src/nested/preprocess-RNA/combine-data-2/gbm-v-hc.csv")
# write.csv(met.v.hc, "./src/nested/preprocess-RNA/combine-data-2/met-v-hc.csv")



