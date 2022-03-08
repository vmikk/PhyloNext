#!/usr/bin/Rscript


## Additional filtering of GBIF occurrences and spatial binning



############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the directory with pre-filtered GBIF occurrences in Parquet format"),
  make_option(c("-s", "--specieskey"), action="store", default=NA, type='character', help="GBIF species ID (optional)"),

  ## Additional filters
  make_option(c("-l", "--terrestrial"), action="store", default=NA, type='character', help="Remove non-terrestrial occurrences, provide land polygon in sf-format"),

  ## DBSCAN options
  make_option(c("-d", "--dbscan"), action="store", default=FALSE, type='logical', help="Remove spatial outliers with density-based clustering"),
  make_option(c("-e", "--epsilon"), action="store", default=700, type='double', help="DBSCAN parameter epsilon, km"),
  make_option(c("-p", "--minpts"), action="store", default=3, type='double', help="DBSCAN min number of points"),

  ## Spatial aggregation
  make_option(c("-r", "--resolution"), action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),

  make_option(c("-t", "--threads"), action="store", default=2L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output directory")
  )


## Validation of the required argiments
if(is.na(opt$input)){
  cat("Input is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output directory is not specified.\n", file=stderr())
  stop()
}


## Assign variables
INPUT <- opt$input
SPECIESKEY <- opt$specieskey
TERRESTRIAL <- opt$terrestrial

DBSCAN <- opt$dbscan
DBSCAN_EPS <- as.numeric(opt$epsilon)
DBSCAN_PTS <- as.numeric(opt$minpts)

RESOLUTION <- as.integer(opt$resolution)
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input occurrences: ", INPUT, "\n", sep=""))
cat(paste("GBIF specieskey: ", SPECIESKEY, "\n", sep=""))
cat(paste("Terrestrial data: ", TERRESTRIAL, "\n", sep=""))

cat(paste("Perform DBSCAN: ", DBSCAN, "\n", sep=""))
if(DBSCAN == TRUE){
  cat(paste("DBSCAN parameter epsilon: ", DBSCAN_EPS, "\n", sep=""))
  cat(paste("DBSCAN min number of points: ", DBSCAN_PTS, "\n", sep=""))
}


cat(paste("Spatial resolution: ", RESOLUTION, "\n", sep=""))

cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output directory: ", OUTPUT, "\n", sep=""))



############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("arrow")
load_pckg("data.table")
load_pckg("dplyr")
load_pckg("h3")

cat("\n")

## Set CPU thread pool
cat("Number of available CPU threads: ", cpu_count(), "\n")
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
set_cpu_count(CPUTHREADS)           # for libarrow
setDTthreads(threads = CPUTHREADS)  # for data.table


## Load land mask
if(!is.na(TERRESTRIAL)){
  cat("Loading land mask\n")
  TERRESTRIAL <- readRDS(TERRESTRIAL)
}




## Open dataset
cat("Loading Parquet data\n")
ds <- arrow::open_dataset(INPUT)

## Select species and collect the data
cat("Collecting data\n")
if(!is.na(SPECIESKEY)){
  datt <- ds %>%
    filter(specieskey %in% SPECIESKEY) %>% 
    collect()
} else {
  datt <- ds %>%
    collect()
}

## Convert to data.table
setDT(datt)

