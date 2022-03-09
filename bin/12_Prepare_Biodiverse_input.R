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


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("rotl")     # Interface to the 'Open Tree of Life' API
# load_pckg("rgbif")
load_pckg("ape")

cat("\n")


## Set CPU thread pool
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table

## Start local cluster
if(CPUTHREADS > 1){
  load_pckg("doFuture")
  registerDoFuture()
  plan(multicore, workers = CPUTHREADS)
  options(future.globals.maxSize = 1e10)

  parall <- TRUE
} else {
  parall <- FALSE
}


############################################## Main pipeline

## Load phylogenetic tree
if(!is.na(PHYTREE)){
  
  ## Function to extract file extension
  file_extension <- function(x){
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }

  if(file_extension(PHYTREE) %in% c("nwk", "newick", "Newick", "tre", "tree")){

    ## Load Newick tree
    cat("Loading pre-computed phylogenetic tree\n")
    TREE <- ape::read.tree(PHYTREE)
  
  } else if(PHYTREE %in% "API") {
    
    cat("Fetching phylogenetic tree from Open Tree API\n")
    # ////////////////////////////////////////  -- TO DO

  } else {

    stop("Unknown tree option is specified. Please check `--phytree` argument.\n")

  }
}



## Load and combine filtered data
cat("Loading filtered occurrences\n")
fls <- list.files(
  path = INPUT,
  pattern = ".RData$",
  full.names = T, recursive = T)

datt <- alply(.data = fls, .margins = 1,
  .fun = function(z){ readRDS(z) },
  .progress = "none", .parallel = parall)

cat("..Merging filtered occurrences from different specieskeys\n")
datt <- rbindlist(datt)

## Check data
if(!nrow(datt) > 1){
  stop("No species occurrences found. Please check the `--input` argument.\n")
}

## Data summary
cat("..Total number of records: ", nrow(datt), "\n")
cat("..Number of unique specieskeys: ", length(unique(datt$specieskey)), "\n")
cat("..Number of unique gridcells: ", length(unique(datt$H3)), "\n")


