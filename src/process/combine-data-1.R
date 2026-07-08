suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(openxlsx))

# The purpose of this file is to first process the protein data
# separately - i.e. change the column names, and 
# remove features with >50% missingness. 
# 2 new files will be generated for HCBB and SBTB.

# ================== HANDLE THE COLUMN NAMES ================

hcbb <- read.table("data/HCBB.gg_matrix.tsv", header=T, sep="\t", check.names=F, row.names=1)
# Change the name of the columns
colnames(hcbb) <- gsub("R:\\PRJ-LargeDataAnalyses_20\\PEV_HBB\\MZMLHBB\\20250825_Eclipse_SH_", 
                       "", colnames(hcbb), fixed=T)
colnames(hcbb) <- gsub("R:\\PRJ-LargeDataAnalyses_20\\PEV_HBB\\MZMLHBB\\20250915_Eclipse_SH_", 
                       "", colnames(hcbb), fixed=T)
colnames(hcbb) <- gsub(".mzML", "",
                       colnames(hcbb), fixed=T)
# Remove trailing values after "_"
colnames(hcbb) <- strsplit(colnames(hcbb), split="_") %>%
  lapply(\(x) x[[1]][1]) %>%
  unlist

sbtb <- read.table("data/SBTB2022.gg_matrix.tsv", header=T, sep="\t", check.names=F, row.names=1)
colnames(sbtb) <- gsub("R:\\PRJ-LargeDataAnalyses_20\\PEV_SBTB2022\\MZMLSBTB\\20220714_Eclipse-OT_SH_",
                       "", colnames(sbtb), fixed=T)
colnames(sbtb) <- gsub("R:\\PRJ-LargeDataAnalyses_20\\PEV_SBTB2022\\MZMLSBTB\\20220719_Eclipse-OT_SH_",
                       "", colnames(sbtb), fixed=T)
colnames(sbtb) <- gsub(".mzML", "",
                       colnames(sbtb), fixed=T)
colnames(sbtb) <- gsub("_200ng|_300ng", "",
                       colnames(sbtb))

# ======================================================

# ============= Extract the samples that we need ===========

sample.sheet <- read.xlsx("./data/class-labels.xlsx")
sbtb <- sbtb %>% dplyr::select(any_of(sample.sheet$Sample.ID))
hcbb <- hcbb %>% dplyr::select(any_of(sample.sheet$Sample.ID))

# ============= Remove features with >50% missingness ===========
# This filtration was done by combining the dataset.
# Note that at this stage, we have only combined the datasets
# to check which features have a lot of missingness. We will
# conduct imputation in each dataset separately to reduce the
# bias from batch effects.

combined.data.uncorrected <- cbind(hcbb, sbtb) 

# Remove features with more than 50% missingness
more.than.0.5.i <- which(rowMeans(is.na(combined.data.uncorrected)) > 0.5)
combined.data.uncorrected <- combined.data.uncorrected[-more.than.0.5.i,]

kept.features <- rownames(combined.data.uncorrected)

# ============ Extract the features that have been kept =============

hcbb.filter <- hcbb[rownames(hcbb) %in% kept.features,]
sbtb.filter <- sbtb[rownames(sbtb) %in% kept.features,]

# ============ Save output ====================

if (!dir.exists("./src/process/combine-data-1")) {
  mkdirs("./src/process/combine-data-1")
}

write.csv(hcbb.filter, file="./src/process/combine-data-1/hcbb-combine-data-1.csv")
write.csv(sbtb.filter, file="./src/process/combine-data-1/sbtb-combine-data-1.csv")
