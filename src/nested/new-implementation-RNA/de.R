library(edgeR)
library(tibble)

compute.de <- function(x, y, conditions, fc) {
  model_matrix <- model.matrix(~ 0 + y)
  colnames(model_matrix) <- sub("y", "", colnames(model_matrix), fixed = TRUE)
  
  contr_matrix <- makeContrasts(contrasts = paste0(conditions[2], " - ", conditions[1]),
                                levels = colnames(model_matrix))

  fit <- lmFit(t(x), model_matrix)
  fit <- contrasts.fit(fit, contr_matrix)
  efit <- eBayes(fit)
  top.table <- topTable(efit, n = Inf, sort.by="none") %>%
    rownames_to_column("molecule")
  result <- top.table %>%
    dplyr::select(molecule, logFC, P.Value, adj.P.Val) %>%
    dplyr::rename(Molecule = molecule, adjPVal = adj.P.Val, PVal = P.Value)
  # Select where logFC >= log2(1.2) and PVal <= 0.05
  which(abs(result$logFC) >= log2(fc) & result$PVal <= 0.05)
}
