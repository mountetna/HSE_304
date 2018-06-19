library(edgeR)
library(limma)
library(methods)
library(statmod)

#e.g. of differential gene expression of 2 groups

counts_matrix <- matrix() # genes by samples raw counts matrix
group1 <- list() # list of column names in group 1
group2 <- list() # list of column names in group 2

groups <- factor(
  c(
    rep("group1", length(group1)),
    rep("group2", length(group2))
  ),
  levels = c("group1", "group2"))

# design matrix for differential gene expression
design <- model.matrix(~groups)
colnames(design) <- levels(groups)

reordered_counts < counts_matrix[, c(group1, group2)]
dge <- DGEList(counts = reordered_counts)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# table of genes sorted by P.Value
tt <- topTable(fit, number = Inf, coef = ncol(design), sort.by = "p")
