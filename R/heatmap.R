library(gplots)
library(RColorBrewer)

# standard color ramp for gene expression - green = low expresssion, red = high expression
green_to_red_color_ramp <- colorRampPalette(c("green", "black", "red"))(80)

# input is matrix of genes as rows and samples as columns
# used for z-scoreing the gene values e.g. log(TPM's + .001) or normalized logCPM
scale_gene_matrix <- function(mat) {
  t(scale(t(mat)))
}

# this is useful to remove genes with all 0's across the samples
remove_genes_with_no_exp <- function(mat) {
  rows <- apply(mat, 1, function(row) sum(row != 0) != 0)
  hm_matrix <- hm_matrix[rows, ]
}

# our favorite cluster method ward.D2
# by default the heatmap.2 distance method is euclidean
hclust.ward = function(d) hclust(d,method="ward.D2")


# color breaks for heatmap
# z scored values should have max red > 1.95 and max green < 1.95
breaks <- (-40:40 / 20)
breaks[81] = 10
breaks[1] = -10

# e.g. heatmap
heatmap.2(
  mat, # this is where the genes by samples matrix goes
  trace = "none",
  key = TRUE, # display the color key
  breaks = breaks,
  hclust = hclust.ward,
  col = green_to_red_color_ramp,
  srtCol = 45,
  margins = c(9, 8)#,
  #ColSideColors, - argument to put in color labels for columns c("blue", "yellow", "blue")
  #RowSideColors - argument to put in color labels for rows c("blue", "yellow", "blue")
)

# color annotation labels
legend("topright",
       legend=c("high", "low", "", "label 1", "label 2"), #labels
       fill=c("red", "green", "white", "purple", "blue"), #colors
       border=FALSE,
       bty="n",
       y.intersp = 0.7, cex=0.7)

# e.g. getting the reordered columns or rows from the heatmap
orig_mat <- matrix() #some gene by sample matrix
hm <- heatmap.2(orig_mat)
reordered_rows <- rownames(orig_mat)[hm$rowInd]
reordered_cols <- rownames(orig_mat)[hm$colInd]
