#############################################################
# Script: softThresholding_benchmark.R
# Author: Zeynep Coskun
# Project: WGCNA Benchmarking
# Goal: Identify optimal soft-thresholding power for 
#       network construction across multiple 
#       correlation functions and network types.
# Based on: WGCNA tutorials (Langfelder & Horvath), adapted 
#           and modularized into reusable code.
#############################################################

# ------------------------------
# 1. Setup
# ------------------------------
# Start a fresh R session before running.
# Set working directory to your project folder:
setwd("~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis")

# Load WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)

# Load pre-processed expression dataset (gene expression matrix)
# File was generated in the previous data cleaning step.
lnames = load(file = "BenchmarkData-0-dataInput.RData")
lnames  

# ------------------------------
# 2. Define reusable function
# ------------------------------
# This function benchmarks soft-thresholding power across:
#   - Multiple correlation functions: Pearson (standard),
#    Biweight midcorrelation (bicor, robust to outliers),
#    Kendall, Spearman
#   - Multiple network types: signed, unsigned, signed hybrid

# Output:
#   - Summary table of optimal parameters
#   - Results list with detailed pickSoftThreshold output


benchmarkSoftThreshold = function(exprData, 
                                  powers = c(1:10, seq(12, 20, 2)),
                                  networkTypes = c("signed", "unsigned", "signed hybrid"),
                                  corFns = list(
                                    pearson = cor,
                                    bicor = bicor,
                                    kendall = function(x, y, ...) cor(x, y, method = "kendall", ...),
                                    spearman = function(x, y, ...) cor(x, y, method = "spearman", ...)
                                  ),
                                  r2_threshold = 0.8,
                                  connectivity_threshold = 10) {
  
  results = list() # Empty list to store results
  summaryTable = data.frame() # store summary across all runs
  
  for (net in networkTypes) {
    for (corName in names(corFns)) {
      corFunc = corFns[[corName]]
      label = paste(net, corName, sep = "_") # e.g., "signed_pearson"
      message("Running for: ", label)
      
      # Run WGCNA soft-thresholding for all parameter combinations
      sft = pickSoftThreshold(exprData,
                               powerVector = powers,
                               networkType = net,
                               corFnc = corFunc,
                               verbose = 5)
      results[[label]] = sft # Store results in list
      
      # Evaluate results and pick best parameters
      # Criteria:
      # - Prefer R² ≥ 0.8 (good scale-free fit)
      # - Mean connectivity ≥ 10
      # - If multiple, choose the lowest power that satisfies conditions
      # - Otherwise, take the power with the max R²
      fitIndices = sft$fitIndices
      signedR2 = -sign(fitIndices[, 3]) * fitIndices[, 2]
      
      # Filter by thresholds
      validRows = which(signedR2 >= r2_threshold & 
                           fitIndices[, 5] >= connectivity_threshold)
      
      if (length(validRows) > 0) {
        # Choose lowest power among valid rows
        bestRowIndex = validRows[which.min(fitIndices[validRows, 1])]
      } else {
        # If no valid row, choose power with max R²
        bestRowIndex = which.max(signedR2)
      }
      
      bestRow = fitIndices[bestRowIndex, ]
      # Append to summary table
      summaryTable = rbind(summaryTable, data.frame(
        Setting = label,
        NetworkType = net,
        Correlation = corName,
        Power = bestRow[1],
        R2 = bestRow[2],
        Slope = bestRow[3],
        MeanConnectivity = bestRow[5]
      ))
    }
  }
  
  return(list(summary = summaryTable, results = results))
}

# ------------------------------
# 3. Apply function to dataset
# ------------------------------
benchmarkResults = benchmarkSoftThreshold(
  exprData = datExpr0)

# Inspect summary
print(benchmarkResults$summary)

# 4. Generate diagnostic plots
# ------------------------------
# (a) Scale independence plots: R² vs power
pdf("scale_independence_plots.pdf", width = 18, height = 16)
par(mfrow = c(3,4))  # multiple panels per page
for (name in names(benchmarkResults$results)) {
  sft = benchmarkResults$results[[name]]
  plot(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R²",
       type = "b", 
       main = paste("Scale independence:", name))
  abline(h = 0.90, col = "red")  # threshold for scale-free fit
}
dev.off()

# (b) Mean connectivity plots
pdf("mean_connectivity_plots.pdf", width = 18, height = 16)
par(mfrow = c(3,4))
for (name in names(benchmarkResults$results)) {
  sft = benchmarkResults$results[[name]]
  plot(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "b",
       main = paste("Mean connectivity:", name))
}
dev.off()

# Save objects for next step of pipeline
save(datExpr0, corFns, networkTypes, benchmarkResults$summary, 
     file = "BenchmarkData-1-dataInput.RData")

#############################################################
# End of script
#############################################################
