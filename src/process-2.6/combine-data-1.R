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
ncol(hcbb)
# Remove all instances of SB12 in hcbb
hcbb <- hcbb[,!grepl("SB12", colnames(hcbb))]

sbtb <- read.table("data/SBTB2022.gg_matrix.tsv", header=T, sep="\t", check.names=F, row.names=1)
colnames(sbtb) <- gsub("R:\\PRJ-LargeDataAnalyses_20\\PEV_SBTB2022\\MZMLSBTB\\20220714_Eclipse-OT_SH_",
                       "", colnames(sbtb), fixed=T)
colnames(sbtb) <- gsub("R:\\PRJ-LargeDataAnalyses_20\\PEV_SBTB2022\\MZMLSBTB\\20220719_Eclipse-OT_SH_",
                       "", colnames(sbtb), fixed=T)
colnames(sbtb) <- gsub(".mzML", "",
                       colnames(sbtb), fixed=T)
colnames(sbtb) <- gsub("_200ng|_300ng", "",
                       colnames(sbtb))
colnames(sbtb)
colnames(sbtb) <- gsub("-", "",
                       colnames(sbtb), fixed=T)
# Remove all instances of SB12, except for SB12_01. This is the
# one recommended for use by Susannah.
sbtb <- sbtb[,!(grepl("SB12", colnames(sbtb)) & (colnames(sbtb) != "SB12_01"))]
# Rename SB12_01 to SB12.
colnames(sbtb)[colnames(sbtb) == "SB12_01"] <- "SB12"

# ======================================================

# ============= Extract the samples that we need ===========

sample.sheet <- read.xlsx("./data/class-labels.xlsx")
sbtb <- sbtb %>% dplyr::select(any_of(sample.sheet$Sample.ID))
hcbb <- hcbb %>% dplyr::select(any_of(sample.sheet$Sample.ID))

sbtb <- log2(sbtb)
hcbb <- log2(hcbb)

# ============= Remove features with >50% missingness ===========
# This filtration was done by combining the dataset.
# Note that at this stage, we have only combined the datasets
# to check which features have a lot of missingness. We will
# conduct imputation in each dataset separately to reduce the
# bias from batch effects.

combined.data.uncorrected <- cbind(hcbb, sbtb) 

# ====================== Side quest: Figure 3
# Remove features with more than 50% missingness in each group separately.
# For Fig 3
gbm.a <- hcbb[,match(sample.sheet$Sample.ID[sample.sheet$Class == "GBM" & sample.sheet$Cohort == "A"],
                    colnames(hcbb))]
met <- hcbb[,match(sample.sheet$Sample.ID[sample.sheet$Class == "MET"], 
                   colnames(hcbb))]
hc <- hcbb[,match(sample.sheet$Sample.ID[sample.sheet$Class == "HC"],
                  colnames(hcbb))]
gbm.b <- sbtb

# Extract the features with less than 50% missingness
extract.non.missing.features <- function(data) {
  more.than.0.5.i <- which(rowMeans(is.na(data)) > 0.5)
  data[-more.than.0.5.i,]
}

gbm.a <- extract.non.missing.features(gbm.a)
met <- extract.non.missing.features(met)
hc <- extract.non.missing.features(hc)
gbm.b <- extract.non.missing.features(gbm.b)
if (!dir.exists("./src/process-2.6/combine-data-1/fig3")) {
  mkdirs("./src/process-2.6/combine-data-1/fig3")
}

# Save each file.
write.csv(gbm.a, file="./src/process-2.6/combine-data-1/fig3/gbm-a.csv")
write.csv(met, file="./src/process-2.6/combine-data-1/fig3/met.csv")
write.csv(hc, file="./src/process-2.6/combine-data-1/fig3/hc.csv")
write.csv(gbm.b, file="./src/process-2.6/combine-data-1/fig3/gbm-b.csv")

# =================================================

more.than.0.5.i <- which(rowMeans(is.na(combined.data.uncorrected)) > 0.5)
combined.data.uncorrected <- combined.data.uncorrected[-more.than.0.5.i,]

kept.features <- rownames(combined.data.uncorrected)

# ============ Extract the features that have been kept =============

hcbb.filter <- hcbb[rownames(hcbb) %in% kept.features,]
sbtb.filter <- sbtb[rownames(sbtb) %in% kept.features,]

# ============ Save output ====================

if (!dir.exists("./src/process-2.6/combine-data-1")) {
  mkdirs("./src/process-2.6/combine-data-1")
}

write.csv(hcbb.filter, file="./src/process-2.6/combine-data-1/hcbb-combine-data-1.csv")
write.csv(sbtb.filter, file="./src/process-2.6/combine-data-1/sbtb-combine-data-1.csv")
