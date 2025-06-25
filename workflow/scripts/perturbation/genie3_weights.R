suppressPackageStartupMessages({
  library(GENIE3)
  library(rhdf5)
})

set.seed(110203) # For reproducibility of results

sn <- get("snakemake", envir = .GlobalEnv)

# Access parameters directly from config
genie3_config <- sn@config[["genie3"]]

opt <- list(
  input = sn@input[["expr"]],
  output_weights = sn@output[["weights"]],
  tree_method = genie3_config[["tree_method"]],
  K = genie3_config[["K"]],
  n_trees = as.integer(genie3_config[["n_trees"]]),
  verbose = genie3_config[["verbose"]]
)

# Read expression matrix from TSV
cat("Reading expression matrix from TSV file\n")
exprMatr <- as.matrix(read.table(
  opt$input,
  header=TRUE,
  row.names=1,
  sep="\t",
  check.names=FALSE
))

# Run GENIE3
cat("Calculating weight matrix\n")
weightMat <- GENIE3(exprMatr, 
                    nCores=sn@threads, 
                    treeMethod=opt$tree_method,
                    K=opt$K,
                    nTrees=opt$n_trees,
                    verbose=opt$verbose)
cat("Weight matrix is calculated\n")

# Output weight matrix to TSV
cat("Writing weight matrix to TSV file\n")
write.table(
  weightMat,
  file=opt$output_weights,
  sep="\t",
  row.names=TRUE,
  quote=FALSE
)
cat("Weight matrix TSV written\n")
