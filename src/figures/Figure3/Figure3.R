library(VennDiagram)
library(gridExtra)
library(grid)

gbm.a <- read.csv("./src/process/combine-data-1/fig3/gbm-a.csv")
gbm.b <- read.csv("./src/process/combine-data-1/fig3/gbm-b.csv")
hc <- read.csv("./src/process/combine-data-1/fig3/hc.csv")
met <- read.csv("./src/process/combine-data-1/fig3/met.csv")

gbm.a.expr <- gbm.a[,1]
gbm.b.expr <- gbm.b[,1]
hc.expr <- hc[,1]
met.expr <- met[,1]

venn.plot <- venn.diagram(
  x=list("GBM Cohort A"=gbm.a.expr,
         "GBM Cohort B"=gbm.b.expr,
         "BMET"=met.expr,
         "HC"=hc.expr),
  category.names = c("", "", "", ""),
  # filename="src/figures/figure-3.png",
  filename=NULL,
  fill=c("#ddffff", "#fbeefa", "#c9eedb", "#ffc7c7"),
  cat.pos=c(180, 180, -45, 45),
  cat.dist=c(0.32, 0.32, 0.12, 0.12),
  cat.box=c(T, T, T, T),
  cat.box.col="black",
  cat.box.lwd=1,
  disable.logging=T,
  cex=1.3
)

boxed_label <- function(
    x, y, label,
    fill = "white",
    width = NULL,
    height = unit(0.08, "npc"),
    default.units = "npc",
    ...
) {
  
  if (is.null(width)) {
    if (label != "HC") {
      width <- grobWidth(textGrob(label, gp = gpar(...))) +
        unit(6, "mm")
    } else {
      width <- grobWidth(textGrob(label, gp = gpar(...))) +
        unit(10, "mm")
    }
  }
  
  grid.rect(
    x = unit(x, default.units),
    y = unit(y, default.units),
    width = width,
    height = height,
    gp = gpar(fill = fill)
  )

  grid.text(
    label,
    x = unit(x, default.units),
    y = unit(y, default.units),
    gp = gpar(...)
  )
}
png(filename = "src/figures/Figure3/figure-3.png", width = 1200, height = 1200, res=275)
grid.newpage()
grid.draw(venn.plot)

boxed_label(
  x = 0.15,
  y = 0.83,
  label = "BMET",
  fill = "#c9eedb",
  fontface="bold"
)

boxed_label(
  x = 0.83,
  y = 0.83,
  label = "HC",
  fill = "#ffc7c7",
  fontface="bold"
)

boxed_label(
  x = 0.17,
  y = 0.2,
  label = "GBM Cohort A",
  fill = "#ddffff",
  fontface="bold"
)

boxed_label(
  x = 0.83,
  y = 0.2,
  label = "GBM Cohort B",
  fill = "#fbeefa",
  fontface="bold"
)

dev.off()



