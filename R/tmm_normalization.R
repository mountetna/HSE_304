library(edgeR)

# normalization for expression counts of a genes by samples matrix
tmm_norm <- function(counts_matrix) {
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <-
    cpm(
      dge,
      normalized.lib.sizes = TRUE,
      log = TRUE,
      prior.count = 0.25
    )
  logCPM
}