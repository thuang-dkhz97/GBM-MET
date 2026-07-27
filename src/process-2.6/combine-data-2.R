library(dplyr)
library(missForest)

# Apply MissForest to the datasets from the previous step.
# Then save into a new file.


# ================== Load data =======================

hcbb.filter <- read.csv("./src/process/combine-data-1/hcbb-combine-data-1.csv",
                        row.names=1, check.names=F)
sbtb.filter <- read.csv("./src/process/combine-data-1/sbtb-combine-data-1.csv",
                        row.names=1, check.names=F)

# ============ Perform missforest imputation ==================
# Perform imputation on each dataset SEPARATELY, so that
# we don't have bias seeping in from batch effects.

hcbb.filter <- missForest(hcbb.filter)$ximp
sbtb.filter <- missForest(sbtb.filter)$ximp

# ============= Save

if (!dir.exists("./src/process/combine-data-2")) {
  mkdirs("./src/process/combine-data-2")
}

write.csv(hcbb.filter, "./src/process/combine-data-2/hcbb-combine-data-2.csv")
write.csv(sbtb.filter, "./src/process/combine-data-2/sbtb-combine-data-2.csv")
