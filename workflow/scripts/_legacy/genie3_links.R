suppressPackageStartupMessages({
  library(GENIE3)
})

sn <- get("snakemake", envir = .GlobalEnv)
opt <- list(
  input_weights = sn@input[["weights"]],
  output_links = sn@output[["links"]],
  top_n_links = sn@params[["top_n_links"]],
  threshold = sn@params[["threshold"]]
)

# Check for mutually exclusive parameters
if (!is.null(opt$top_n_links) && !is.null(opt$threshold)) {
  stop("Error: top_n_links and threshold are mutually exclusive. Please set only one of them in the config.")
}

# Read weight matrix
cat("Reading weight matrix\n")
weightMat <- as.matrix(read.table(opt$input_weights, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
cat("Weight matrix loaded\n")

# Obtain ranked list of regulatory links
cat("Getting link list\n")

if (!is.null(opt$top_n_links) && opt$top_n_links != "null") {
  cat(paste("Using top_n_links =", opt$top_n_links, "\n"))
  linkList <- getLinkList(weightMat, reportMax = as.numeric(opt$top_n_links))
} else if (!is.null(opt$threshold) && opt$threshold != "null") {
  cat(paste("Using threshold =", opt$threshold, "\n"))
  linkList <- getLinkList(weightMat, threshold = as.numeric(opt$threshold))
} else {
  cat("No filtering parameters specified, using default getLinkList behavior\n")
  linkList <- getLinkList(weightMat)
}

cat("Link list is obtained\n")

# Write output
cat("Writing link list\n")
write.table(linkList,
            file=opt$output_links,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)
cat("Link list is written\n")
