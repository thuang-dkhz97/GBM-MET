library(openxlsx)
library(dplyr)
library(umap)
library(ggplot2)
library(sva)

plot.figs <- function(data, phenotype, cohort.labels, conditions,
                      plot.save.dir, extra.name, colours) {
  label <- phenotype$Label[match(colnames(data), phenotype$Sample)]
  cohort.matched <- cohort.labels$Cohort[match(colnames(data), cohort.labels$Sample.ID)]
  print(cohort.matched)
  print(colnames(data))
  print(label)
  # best_markers <- best_markers %>%
  #   column_to_rownames("Molecule")
  
  # print(colnames(best_markers))
  # From data, extract the best markers
  # print(rownames(data))
  # full.data <- data
  # data <- data[rownames(data) %in% rownames(best_markers), ]
  
  # print(t(data))
  data <- t(data) %>% as.data.frame
  data$group <- label
  
  # And also to the best set
  result <- umap(data %>% dplyr::select(!all_of(c("group"))),
                 n_neighbors=15)
  dim_red_df <- data.frame(x=result$layout[,1], y=result$layout[,2],
                           group=data$group, cohort=cohort.matched)
  # print(dim_red_df)
  
  ggplot(dim_red_df, aes(x=x, y=y, fill=group, shape=cohort)) +
    geom_point(size=3, stroke=1) +
    scale_shape_manual(name = "Data Cohort", values = c(21, 22),
                       labels=c("Cohort 1", "Cohort 2")) +
    scale_fill_manual(name = "Condition", values = colours) +
    guides(fill = guide_legend(override.aes = list(shape = 21, colour = colours))) +
    theme_bw() +
    theme(legend.text=element_text(size=18),
          legend.title=element_text(size=20),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))+
    labs(y="UMAP2", x="UMAP1")
  comp.name <- paste(rev(conditions), collapse="-vs-")
  ggsave(file.path(plot.save.dir, paste0(comp.name, "-", extra.name, "-umap.png")), units="cm", width=35, height=20)
  
}

# Extract the before batch effect correction

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

# cohort.labels <- read.xlsx("./data/class-labels.xlsx")
save.dir <- file.path("src", "figures", "SuppFigure3", "figures")
R.utils::mkdirs(save.dir)
# ========= No BE
gbm.v.met.labels <- read.table("src/process-2.6/combine-data-3/GBMVsMET_output_labels.txt", header=T)
plot.figs(gbm.v.met, gbm.v.met.labels, sample.data, c("MET", "GBM"),
          save.dir, "no-be", colours=c("#9f21ee", "#00fefe"))

gbm.v.hc.labels <- read.table("src/process-2.6/combine-data-3/GBMVsHC_output_labels.txt", header=T)
plot.figs(gbm.v.hc, gbm.v.hc.labels, sample.data, c("HC", "GBM"),
          save.dir, "no-be", colours=c("#09fe05", "#00fefe"))

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


# ========= Apply BE
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

plot.figs(gbm.v.met.norm, gbm.v.met.labels, sample.data, c("MET", "GBM"),
          save.dir, "combat-be", colours=c("#00fefe", "#9f21ee"))

plot.figs(gbm.v.hc.norm, gbm.v.hc.labels, sample.data, c("HC", "GBM"),
          save.dir, "combat-be", colours=c("#00fefe", "#09fe05"))

