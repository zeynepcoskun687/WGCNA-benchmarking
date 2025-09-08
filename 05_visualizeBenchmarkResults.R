#############################################################
# script: visualizeBenchmarkResults.R
# Author: Zeynep Coskun
# Project: WGCNA Benchmarking
# Purpose:
#   - Compare hub gene overlap across networks
#   - Quantify and visualize similarity of hub sets
#   - Align module structures to a reference network
#   - Explore cross-network eigengene correlations
#   - Summarize hub gene frequency across methods
# Outputs:
#   - Similarity heatmaps
#   - Dendrogram with aligned module colors
#   - Barplot of hub gene frequencies

#############################################################

# -----------------------------------------------------------
# 1. Hub gene overlap across methods
# -----------------------------------------------------------
# Goal:
#   - Identify whether the same genes appear as top hubs
#     across different methods and modules.
#   - Build a gene × Method_Module matrix of kME values.

combinedHubData$MethodModule = paste(combinedHubData$Method,
                                     combinedHubData$Module,
                                     sep = "_")

all_genes = unique(combinedHubData$Gene)
allMethodModules = unique(combinedHubData$MethodModule)

kme_matrix = matrix(
  NA,
  nrow = length(all_genes),
  ncol = length(allMethodModules),
  dimnames = list(all_genes, allMethodModules)
)

# Fill matrix with kME values
for (i in seq_len(nrow(combinedHubData))) {
  gene = combinedHubData$Gene[i]
  methodMod = combinedHubData$MethodModule[i]
  row_idx = match(gene, all_genes)
  col_idx = match(methodMod, allMethodModules)
  kme_matrix[row_idx, col_idx] = combinedHubData$kME[i]
}


# -----------------------------------------------------------
# 2. Jaccard similarity of hub sets
# -----------------------------------------------------------
# Goal:
#   - Quantify overlap in hub gene lists between all
#     Method_Module combinations.
#   - Use Jaccard index (intersection / union).
#   - Visualize as a labeled heatmap.

hubLists = split(combinedHubData$Gene, combinedHubData$MethodModule)
methods = names(hubLists)

jaccardMat = matrix(
  0,
  nrow = length(methods),
  ncol = length(methods),
  dimnames = list(methods, methods)
)

for (i in seq_along(methods)) {
  for (j in seq_along(methods)) {
    intersectSize = length(intersect(hubLists[[i]], hubLists[[j]]))
    unionSize = length(union(hubLists[[i]], hubLists[[j]]))
    jaccardMat[i, j] = intersectSize / unionSize
  }
}

# Heatmap of Jaccard similarity
labeledHeatmap(
  Matrix = jaccardMat,
  xLabels = colnames(jaccardMat),
  yLabels = rownames(jaccardMat),
  colors = colorRampPalette(c("white", "lightblue", "blue"))(50),
  textMatrix = signif(jaccardMat, 2),
  colorLabels = FALSE,
  setStdMargins = FALSE,
  cex.lab.x = 0.6,
  cex.lab.y = 0.6,
  cex.text = 0.5,
  main = "Jaccard Similarity of Hub Gene Sets"
)


# -----------------------------------------------------------
# 3. Align module structures across methods
# -----------------------------------------------------------
# Goal:
#   - Use signed_pearson as the reference network.
#   - Align module colors from other networks using
#     matchLabels().
#   - Plot a single dendrogram with aligned color bars.

refMethod = "signed_pearson"

# Reconstruct full signed_pearson network without block splitting
net_full = blockwiseModules(
  datExpr0,
  power = 16,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  maxBlockSize = 6000,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

# Extract key components
moduleLabels1 = net_full$colors
moduleColors1 = labels2colors(net_full$colors)
geneTree1 = net_full$dendrograms[[1]]
MEs1 = net_full$MEs
refColors = moduleColors

# Initialize matched color list
allMethods = c("signed_pearson", "signed_bicor",
               "signed_hybrid_pearson", "signed_hybrid_bicor",
               "unsigned_pearson")

#create matched colors list to store the color assignments to each gene 
#in each network:

matchedColorsList = list("Reference (signed_pearson)" = matrix(refColors, ncol = 1))

for (method in allMethods) {
  if (method == refMethod) next
  load(paste0("~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis/BenchmarkData-2-", method, ".RData"))
  matched = matchLabels(moduleColors, refColors)
  matchedColorsList[[method]] = matrix(matched, ncol = 1)
}

names(matchedColorsList) = c("Reference (signed_pearson)",
                             "signed_bicor",
                             "signed_hybrid_pearson",
                             "signed_hybrid_bicor",
                             "unsigned_pearson")

# Plot dendrogram with aligned module colors
plotDendroAndColors(
  dendro = geneTree1,
  colors = do.call(cbind, matchedColorsList),
  groupLabels = names(matchedColorsList),
  main = "Module Comparison Aligned to signed_pearson Dendrogram",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

# Save alignment results
save(matchedColorsList, geneTree1, MEs1, refColors,
     file = "~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis/BenchmarkData-4-moduleComparison.RData")


# -----------------------------------------------------------
# 4. Cross-network eigengene correlations
# -----------------------------------------------------------
# Goal:
#according to the results, we look at cross-network module correspondences 
#(e.g., turquoise ↔ pink, blue ↔ red) and to ask whether the eigengene correlations
#between these modules indicate that they might be biologically similar 
#(or even candidates for merging):

# Get color vectors from matched colors list for the modules to compare
colors_ref = matchedColorsList[["Reference (signed_pearson)"]]
colors_bicor = matchedColorsList[["signed_bicor"]]

# Compute module eigengenes
ME_ref = moduleEigengenes(datExpr0, colors_ref)
ME_bicor = moduleEigengenes(datExpr0, colors_bicor)

# Select modules of interest
ME_ref_turquoise = ME_ref[["eigengenes"]][["MEturquoise"]]
ME_bicor_pink = ME_bicor[["eigengenes"]][["MEpink"]]
ME_ref_blue = ME_ref[["eigengenes"]][["MEblue"]]
ME_bicor_red = ME_bicor[["eigengenes"]][["MEred"]]

# Compute correlations
cor_value_pink = cor(ME_ref_turquoise, ME_bicor_pink, use = "p")
p_value_pink   = corPvalueStudent(cor_value_pink, nrow(datExpr0))

cor_value_red = cor(ME_ref_blue, ME_bicor_red, use = "p")
p_value_red   = corPvalueStudent(cor_value_red, nrow(datExpr0))

# Example interpretation:
# turquoise (signed_pearson) vs. pink (signed_bicor)
#   Correlation = 0.175, p = 0.5177 → not significant
# blue (signed_pearson) vs. red (signed_bicor)
#   Correlation = 0.274, p = 0.305 → not significant

#The eigengene correlation between the turquoise module in the signed_pearson
#network and the pink module in the signed_bicor network was low (r = 0.175) and 
#not statistically significant (p = 0.5177), suggesting no meaningful co-expression 
#similarity between these modules across samples.
#We therefore cannot conclude that these modules represent the same underlying
#expression pattern, despite visual similarity in dendrogram-based module color mapping.

# -----------------------------------------------------------
# 5. Hub gene frequency analysis
# -----------------------------------------------------------
# Goal:
#   - Identify which genes are most frequently selected
#     as hub genes across all networks.
#   - Visualize top recurrent hub genes.

HubGeneFrequency = as.data.frame(table(combinedHubData$Gene))
colnames(HubGeneFrequency) = c("Gene", "Frequency")

# Select top 52 most frequent hub genes
topGenes = HubGeneFrequency[order(-HubGeneFrequency$Frequency), ][1:52, ]

# Plot frequency barplot
pdf("HubGeneFrequencyPlot.pdf", width = 10, height = 12)
barplot(
  topGenes$Frequency,
  names.arg = topGenes$Gene,
  las = 2,
  col = "skyblue",
  main = "Top Most Frequent Hub Genes",
  ylab = "Frequency",
  cex.names = 0.6
)
dev.off()

# Save frequency data
save(HubGeneFrequency, topGenes,
     file = "~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis/BenchmarkData-5-hubGeneFrequency.RData")

#############################################################
# End of script
#############################################################