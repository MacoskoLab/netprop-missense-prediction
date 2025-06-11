# GENIE3 wrapper script
# Usage: Rscript genie3.r -i expr_matrix.tsv -o output_links.tsv [-r regulators.txt] [--ntrees N] [--ncores C]

suppressPackageStartupMessages({
  library(GENIE3)
})

sn <- get("snakemake", envir = .GlobalEnv)
opt <- list(
  input = sn@input[["expr"]],
  output = sn@output[["links"]],
  ntrees = if (!is.null(sn@params[["n_trees"]])) as.integer(sn@params[["n_trees"]]) else 1000,
  min_size = if (!is.null(sn@params[["min_size"]])) as.integer(sn@params[["min_size"]]) else NULL,
  max_depth = if (!is.null(sn@params[["max_depth"]])) as.integer(sn@params[["max_depth"]]) else NULL
)

# Read expression matrix
exprMatr <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, sep="\t", check.names=FALSE))

# Run GENIE3
weightMat <- GENIE3(exprMatr,
                    nTrees=opt$ntrees,
                    min_size=opt$min_size,
                    max_depth=opt$max_depth)

# Obtain ranked list of regulatory links
linkList <- getLinkList(weightMat)

# Write output
write.table(linkList,
            file=opt$output,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)
