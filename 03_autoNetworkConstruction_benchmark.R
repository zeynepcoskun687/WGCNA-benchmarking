############################################################
# script: autoNetworkConstruction_benchmark.R
# Author: Zeynep Coskun
# Project: WGCNA Benchmarking
# Description:
#   - Constructs networks with different correlation + network type combinations
#   - Detects gene co-expression modules for each method
#   - Saves results for comparison and visualization
# Based on:  WGCNA tutorials (Langfelder & Horvath), adapted.
############################################################

# -----------------------------------------------------------
# 1. Setup
# -----------------------------------------------------------
# Start a fresh R session before running.
# Set working directory to your project folder:
setwd("~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis")

#load the saved objects from the previous script 
#which are the data expression matrix and the summary table of soft thresholding 
#powers for different method combinations

load("BenchmarkData-1-dataInput.RData")

# Load WGCNA package
library(WGCNA)
options(stringsAsFactors = FALSE)


# -----------------------------------------------------------
# 2. Define benchmarking methods
# -----------------------------------------------------------
# Each method specifies:
#   - network type: signed, signed hybrid, or unsigned
#   - correlation type: Pearson or biweight midcorrelation (bicor)
#   - soft-thresholding power chosen from previous analysis
methods = list(
  "signed_hybrid_bicor"   = list(networkType = "signed hybrid", corType = "bicor",   power = 9),
  "signed_hybrid_pearson" = list(networkType = "signed hybrid", corType = "pearson", power = 8),
  "signed_bicor"          = list(networkType = "signed",        corType = "bicor",   power = 20),
  "signed_pearson"        = list(networkType = "signed",        corType = "pearson", power = 16),
  "unsigned_pearson"      = list(networkType = "unsigned",      corType = "pearson", power = 14)
)


# -----------------------------------------------------------
# 3. Construct networks and detect modules
# ------------------------------------------------------
for (name in names(methods)) {
  print(paste("Running:", name))
  method = methods[[name]]
  
  net = blockwiseModules(
    datExpr0,
    power = method$power,
    networkType = method$networkType,
    corType = method$corType,
    TOMType = ifelse(method$networkType == "unsigned", "unsigned", "signed"), 
    # Note: TOMType cannot be "signed hybrid" → must be "signed" in this case
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = paste0("Benchmark_TOM_", name),
    verbose = 3
  )
  
  # Save each network object as Rdata files
  save(net, file = paste0("net_", name, ".RData"))
}
# Free memory
gc()


# -----------------------------------------------------------
# 4. Summarize and visualize detected modules
# -----------------------------------------------------------
methods = c(
  "signed_bicor",
  "signed_pearson",
  "signed_hybrid_bicor",
  "signed_hybrid_pearson",
  "unsigned_pearson"
)

for (method in methods) {
  print(paste("Processing:", method))
  
  # Load the saved network object
  load(paste0("net_", method, ".RData"))
  
  # Summarize module sizes
  print(table(net$colors))
  
  # Convert numeric module labels → color names
  mergedColors = labels2colors(net$colors)
  
  # Plot gene dendrogram with module colors
  plotDendroAndColors(
    net$dendrograms[[1]],
    mergedColors[net$blockGenes[[1]]],
    paste("Module colors:", method),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  
  # Extract module data
  moduleLabels = net$colors
  moduleColors = mergedColors
  MEs = net$MEs
  geneTree = net$dendrograms[[1]]
  
  # Save results for further use
  save(MEs, moduleLabels, moduleColors, geneTree,
       file = paste0("~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis/BenchmarkData-2-", method,".RData"))
}

#############################################################
# End of script
#############################################################