############################################################
# Script: dataInput_preprocess.R
# Author: Zeynep Coskun
# Project: WGCNA Benchmarking
# Purpose: Clean and preprocess expression data for WGCNA
# Output: BenchmarkData-0-dataInput.RData
# Based on: WGCNA tutorials (Langfelder & Horvath), adapted.
############################################################

#-----------------------------------------------------------
# 1. Setup
#-----------------------------------------------------------

# Change working directory to project folder
setwd("~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis") 

# Load required package
library(WGCNA)
options(stringsAsFactors = FALSE)

#-----------------------------------------------------------
# 2. Load expression data set
#-----------------------------------------------------------

# Expression data: rows = genes, columns = samples
# This data set does not contain extra metadata to remove
load("/Users/zeynepmacbook/Downloads/CraftedFor_WGCNA.RData")

# Inspect dimensions and sample names
dim(bCorrected)
names(bCorrected)

#-----------------------------------------------------------
# 3. Data formatting
#-----------------------------------------------------------

# Transpose the data so that:
#   rows = samples
#   columns = genes
datExpr0 = as.data.frame(t(bCorrected))

#-----------------------------------------------------------
# 4. Quality control
#-----------------------------------------------------------

# Check for missing values, genes with zero variance, or bad samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK   # Returns TRUE if dataset passes QC

#-----------------------------------------------------------
# 5. Sample clustering (to detect outliers)
#-----------------------------------------------------------

# Calculate Euclidean distance between samples and perform hierarchical clustering
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the clustering dendrogram
par(cex = 0.6)              # shrink text size
par(mar = c(0,4,2,0))       # set margins
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Note: No clear outliers are visible (sample size is already small)

#-----------------------------------------------------------
# 6. Save cleaned dataset
#-----------------------------------------------------------

save(datExpr0, file = "BenchmarkData-0-dataInput.RData")
