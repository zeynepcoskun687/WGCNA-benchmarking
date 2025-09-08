#############################################################
# script: extractHubGenes_benchmark.R
# Author: Zeynep Coskun
# Project: WGCNA Benchmarking
# Purpose: Identify and compare top hub genes per module for multiple WGCNA networks
#generated with different parameters
# Outputs:
#   - A list of data frames per method
#   - A combined data frame of all methods
#Based on: WGCNA tutorials (Langfelder & Horvath), adapted 
#           and modularized into reusable code.
#############################################################


# ------------------------------
# 1. Setup
# ------------------------------
# Start a fresh R session before running
# Set working directory to your project folder
setwd("~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis")

# Load required packages
library(WGCNA)
options(stringsAsFactors = FALSE)

# Load the original expression data
load(file = "BenchmarkData-1-dataInput.RData")  # contains datExpr0, corFns,
#networkTypes, summaryTable

# ------------------------------
# 2. Define reusable function
# ------------------------------
# This function extracts top hub genes per module for multiple WGCNA networks
# Inputs:
#   - exprData: expression matrix
#   - methods: character vector of method names corresponding to saved .RData files
#   - dataDir: path where module .RData files are stored
#   - topN: number of top hub genes to extract per module
# Outputs:
#   - hubGeneTables: list of data frames (one per method)
#   - combinedHubData: single data frame containing all methods

extractHubGenes = function(exprData,
                            methods = c("signed_pearson", "signed_bicor", 
                                        "signed_hybrid_pearson", "signed_hybrid_bicor", 
                                        "unsigned_pearson"),
                            dataDir = "data/",
                            topN = 10) {
  
  hubGeneTables = list()  # store results per method
  
  for (method in methods) {
    cat("Processing:", method, "\n")
    
    # Load module data (moduleLabels, moduleColors, MEs) with the name it was saved
    load(file.path(dataDir, paste0("BenchmarkData-2-", method, ".RData")))
    
    # Compute module eigengenes
    MEs0 = moduleEigengenes(exprData, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    modNames = substring(names(MEs), 3)
    
    # Determine correlation function based on method
    if (grepl("bicor", method)) {
      corFnc = bicor  
    } 
    else {
      corFnc = cor 
    }
    
    # Identify if network is unsigned
    isUnsigned = grepl("unsigned", method)
    
    # Compute module membership (kME)
    MM = as.data.frame(corFnc(exprData, MEs, use = "p"))
    if (nrow(MM) == ncol(exprData)) {
      rownames(MM) = colnames(exprData)
    } else {
      stop("Dimension mismatch between MM rows and gene names!")
    }
    if (isUnsigned) {
      MM = abs(MM) #Unsigned networks use absolute correlation
    }
    
    names(MM)= paste("MM", modNames, sep="")
    rownames(MM) = colnames(exprData)
    
    # Extract top N hub genes per module
    topHubTable = data.frame()
    for (mod in modNames) {
      mm_col = paste0("MM", mod)
      genesInModule = which(moduleColors == mod)
      kme_values = MM[genesInModule, mm_col]
      gene_names = rownames(MM)[genesInModule]
      top_idx = order(kme_values, decreasing = TRUE)[1:min(topN, length(kme_values))]
      top_genes = gene_names[top_idx]
      top_kmes = kme_values[top_idx]
      df = data.frame(Method = method, Module = mod, Gene = top_genes, 
                      kME = round(top_kmes,4), stringsAsFactors = FALSE)
      topHubTable = rbind(topHubTable, df)
    }
    
    hubGeneTables[[method]] = topHubTable
  }
  
  combinedHubData = do.call(rbind, hubGeneTables)
  
  return(list(hubGeneTables = hubGeneTables,
              combinedHubData = combinedHubData))
}

# ------------------------------
# 3. Apply function to dataset
# ------------------------------
hubResults = extractHubGenes(
  exprData = datExpr0,
  dataDir = "~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis",  # relative path to saved module data
  topN = 10           # number of hub genes per module
)

# ------------------------------
# 4. Inspect results
# ------------------------------
# View combined hub genes
hubResults$combinedHubData

# Access hub genes from a specific method
hubResults$hubGeneTables$signed_pearson

# ------------------------------
# 5. Save results for downstream analysis
# ------------------------------
save(hubResults, file = "~/Desktop/WGCNA_GeneNetwork/benchmarkAnalysis/BenchmarkData-3-hubGenes.RData")

#############################################################
# End of script
#############################################################