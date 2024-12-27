
setwd("C:/Users/uabic/Desktop")

# Load the necessary libraries
library(SingleCellExperiment)
library(BASiCS)
library(readr)

# Load the data
data <- read.csv("IDHWT_log_normalized.csv")

# Extract gene names and sample data
genes <- data$Gene
counts <- as.matrix(data[, -1])
rownames(counts) <- genes

# Create placeholder batch information with two batches
set.seed(123)  # For reproducibility
n_samples <- ncol(counts)
batch_info <- data.frame(BatchInfo = sample(c(1, 2), n_samples, replace = TRUE))

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = batch_info)

# Check the SingleCellExperiment object
print(sce)

# Use the created SingleCellExperiment object for BASiCS
Data <- sce

# Run BASiCS MCMC
Chain <- BASiCS_MCMC(Data = Data, 
                     N = 100, 
                     Thin = 5, 
                     Burn = 50, 
                     Regression = FALSE, 
                     WithSpikes = FALSE, 
                     PrintProgress = FALSE)

# Summarize the MCMC chain
ChainSummary <- Summary(Chain)

# Inspect the structure of ChainSummary
str(ChainSummary)

# Visualize the summaries
par(mfrow = c(2, 2))
plot(ChainSummary, Param = "mu", main = "All genes", log = "y")
plot(ChainSummary, Param = "mu", Genes = 1:10, main = "First 10 genes")
plot(ChainSummary, Param = "delta", main = "All genes")
plot(ChainSummary, Param = "delta", Genes = c(2, 5, 10, 50), main = "5 custom genes")

# Check if 'phi' is a valid parameter before plotting
if ("phi" %in% colnames(ChainSummary$Table)) {
  par(mfrow = c(1, 2))
  plot(ChainSummary, Param = "phi")
  plot(ChainSummary, Param = "s", Cells = 1:5)
} else {
  message("'phi' is not a valid parameter in ChainSummary.")
}

# Scatterplot of posterior estimates for gene-specific parameters
par(mfrow = c(1, 2))
plot(ChainSummary, Param = "mu", Param2 = "delta", log = "x", SmoothPlot = FALSE)
plot(ChainSummary, Param = "mu", Param2 = "delta", log = "x", SmoothPlot = TRUE)

# Detect highly variable genes (HVG) and lowly variable genes (LVG)
HVG <- BASiCS_DetectHVG(Chain, VarThreshold = 0.6, Plot = TRUE)
tryCatch({
  LVG <- BASiCS_DetectLVG(Chain, VarThreshold = 0.2, Plot = TRUE)
}, error = function(e) {
  message("EFDR calibration failed for LVG detection. Using default ProbThreshold.")
  LVG <- BASiCS_DetectLVG(Chain, VarThreshold = 0.2, Plot = TRUE, ProbThreshold = 0.5)
})

# Inspect the structure of HVG and LVG objects
str(HVG)
str(LVG)

# Access the results of these tests
if (inherits(HVG, "BASiCS_ResultVG")) {
  HVG_results <- as.data.frame(HVG@Table)
  head(HVG_results)
}

if (inherits(LVG, "BASiCS_ResultVG")) {
  LVG_results <- as.data.frame(LVG@Table)
  head(LVG_results)
}


# Save the results to CSV files
write.csv(as.data.frame(HVG_results), file = "HVG_results.csv", row.names = FALSE)
write.csv(as.data.frame(LVG_results), file = "LVG_results.csv", row.names = FALSE)

# Summarize the MCMC chain
ChainSummary <- Summary(Chain)

# Plot the summaries and add points for HVG and LVG
par(mfrow = c(2, 2))
plot(ChainSummary, Param = "mu", Param2 = "delta", log = "xy")
if (exists("HVG_results")) {
  with(HVG_results[HVG_results$HVG == TRUE,], points(Mu, Delta, col = 'red'))
}
if (exists("LVG_results")) {
  with(LVG_results[LVG_results$LVG == TRUE,], points(Mu, Delta, col = 'blue'))
}
