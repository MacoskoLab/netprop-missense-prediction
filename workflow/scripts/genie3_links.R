suppressPackageStartupMessages({
  library(GENIE3)
})

sn <- get("snakemake", envir = .GlobalEnv)
opt <- list(
  input_weights = sn@input[["weights"]],
  output_links = sn@output[["links"]]
)

# Read weight matrix
cat("Reading weight matrix\n")
weightMat <- as.matrix(read.table(opt$input_weights, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
cat("Weight matrix loaded\n")

# Obtain ranked list of regulatory links
cat("Getting link list\n")
linkList <- getLinkList(weightMat)
cat("Link list is obtained\n")

# Write output
cat("Writing link list\n")
write.table(linkList,
            file=opt$output_links,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)
cat("Link list is written\n")
