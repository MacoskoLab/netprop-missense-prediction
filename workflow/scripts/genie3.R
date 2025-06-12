set.seed(42) # For reproducibility

suppressPackageStartupMessages({
  library(GENIE3)
})

sn <- get("snakemake", envir = .GlobalEnv)
opt <- list(
  input = sn@input[["expr"]],
  output = sn@output[["links"]],
  ntrees = as.integer(sn@params[["n_trees"]])
)

# Read expression matrix
exprMatr <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, sep="\t", check.names=FALSE))

# Run GENIE3
weightMat <- GENIE3(exprMatr, nCores=sn@threads, nTrees=opt$ntrees, verbose=TRUE)

# Obtain ranked list of regulatory links
linkList <- getLinkList(weightMat)

# Write output
write.table(linkList,
            file=opt$output,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)