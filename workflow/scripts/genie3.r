# GENIE3 wrapper script
# Usage: Rscript genie3.r -i expr_matrix.tsv -o output_links.tsv [-r regulators.txt] [--ntrees N] [--ncores C]

suppressPackageStartupMessages({
  library(optparse)
  library(GENIE3)
})

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to expression matrix (TSV) with genes as rows, samples as columns", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", help="Path to output regulatory links (TSV)", metavar="FILE"),
  make_option(c("-r", "--regulators"), type="character", default=NULL, help="Optional: file listing candidate regulators (one gene per line)", metavar="FILE"),
  make_option(c("-n", "--ntrees"), type="integer", default=1000, help="Number of trees per ensemble [default = %default]", metavar="INT"),
  make_option(c("-c", "--ncores"), type="integer", default=1, help="Number of cores to use [default = %default]", metavar="INT")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  quit(status = 1)
}

# Read expression matrix
exprMatr <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, sep="\t", check.names=FALSE))

# Load candidate regulators if provided
if (!is.null(opt$regulators)) {
  regs <- readLines(opt$regulators)
  weightMat <- GENIE3(exprMatr,
                      regulators=regs,
                      nTrees=opt$ntrees,
                      nCores=opt$ncores)
} else {
  weightMat <- GENIE3(exprMatr,
                      nTrees=opt$ntrees,
                      nCores=opt$ncores)
}

# Obtain ranked list of regulatory links
linkList <- getLinkList(weightMat)

# Write output
write.table(linkList,
            file=opt$output,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)
