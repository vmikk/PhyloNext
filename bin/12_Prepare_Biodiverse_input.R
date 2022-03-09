#!/usr/bin/Rscript

## Prepare input data for Biodiverse

############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the directory with filtered data"),
  make_option(c("-p", "--phytree"), action="store", default=NA, type='character', help="Phylogenetic tree (pre-computed in Newick format or get it from API)"),
  make_option(c("-g", "--taxgroup"), action="store", default="All life", type='character', help="Taxonomy group in OpenTree"),

  make_option(c("-t", "--threads"), action="store", default=2L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$input)){
  stop("Input is not specified.\n")
}
if(is.na(opt$output)){
  stop("Output directory is not specified.\n")
}


## Assign variables
INPUT <- opt$input
PHYTREE <- opt$phytree
TAXGROUP <- opt$taxgroup
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input directory with filtered occurrences: ", INPUT, "\n", sep=""))
cat(paste("Pre-computed phylogenetic tree: ", PHYTREE, "\n", sep=""))
cat(paste("Taxonomy group in OpenTree: ", TAXGROUP, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output directory: ", OUTPUT, "\n", sep=""))

cat("\n")


## Create output directory
dir.create(path = OUTPUT, showWarnings = F, recursive = TRUE)
