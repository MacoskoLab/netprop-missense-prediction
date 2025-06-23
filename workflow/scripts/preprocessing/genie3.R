suppressPackageStartupMessages({
  library(GENIE3)
})

set.seed(123) # For reproducibility of results

sn <- get("snakemake", envir = .GlobalEnv)
opt <- list(
  input = sn@input[["expr"]],
  output_links = sn@output[["links"]],
  output_weights = sn@output[["weights"]],
  tree_method = sn@params[["tree_method"]],
  K = sn@params[["K"]],
  n_trees = as.integer(sn@params[["n_trees"]])
)

# Read expression matrix
exprMatr <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, sep="\t", check.names=FALSE))

# Run GENIE3
cat("Calculating weight matrix\n")
weightMat <- GENIE3(exprMatr, 
                    nCores=sn@threads, 
                    treeMethod=opt$tree_method,
                    K=opt$K,
                    nTrees=opt$n_trees,
                    verbose=TRUE)
cat("Weight matrix is calculated\n")

# Output weight matrix
cat("Writing weight matrix\n")
write.table(weightMat,
            file=opt$output_weights,
            sep="\t",
            row.names=TRUE,
            quote=FALSE)
cat("Weight matrix is written\n")

# Obtain ranked list of regulatory links
cat("Getting link list\n")
linkList <- getLinkList(weightMat)

# Write output
cat("Writing link list\n")
write.table(linkList,
            file=opt$output_links,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)
cat("Link list is written\n")
