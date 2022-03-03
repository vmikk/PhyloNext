#!/usr/bin/Rscript

## Filter GBIF occurrence data and perform spatial binning


## Input data could be located in the cloud.
## E.g., to query GBIF AWS snapshot use: 
##  --input "s3://gbif-open-data-eu-central-1/occurrence/2021-11-01/occurrence.parquet"

## Usage:
# ./Filter_and_aggregate.R \
#    --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet" \
#    --family "Fabaceae" \
#    --country "AU" \
#    --resolution "4" \
#    --output "Fabaceae_in_AU"


############################################## Parse input parameters

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the directory with GBIF snapshot of occurrence records in Parquet format"),
  
  ## Taxonomy filters
  make_option("--phylum", action="store", default=NA, type='character', help="Comma-separated list of phyla to select"),
  make_option("--class", action="store", default=NA, type='character', help="Comma-separated list of classes to select"),
  make_option("--order", action="store", default=NA, type='character', help="Comma-separated list of orders to select"),
  make_option("--family", action="store", default=NA, type='character', help="Comma-separated list of families to select"),
  
  ## Spatial filters
  make_option("--country", action="store", default=NA, type='character', help="Comma-separated list of country codes (e.g., AU,CA)"),
  make_option("--latmin", action="store", default=NA, type='double', help="Minimum latitude"),
  make_option("--latmax", action="store", default=NA, type='double', help="Maximum latitude"),
  make_option("--lonmin", action="store", default=NA, type='double', help="Minimum longitude"),
  make_option("--lonmax", action="store", default=NA, type='double', help="Maximum longitude"),

  ## Additional filters
  make_option("--minyear"), action="store", default=1945, type='integer', help="Minimum year of occurrence (default, 1945)"),
  make_option("--terrestrial", action="store", default=FALSE, type='logical', help="Remove non-terrestrial occurrences"),
  make_option("--noextinct", action="store", default=NA, type='character', help="Remove extinct species (provide a file with extinct specieskeys)"),

  make_option(c("-r", "--resolution"), action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output prefix")
)
opt <- parse_args(OptionParser(option_list=option_list))



## Validation of the required argiments
if(is.na(opt$input)){
  cat("Input is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output prefix is not specified.\n", file=stderr())
  stop()
}

## Assign variables
INPUT <- opt$input

PHYLUM <- opt$phylum
CLASS <- opt$class
ORDER <- opt$order
FAMILY <- opt$family

COUNTRY <- opt$country
LATMIN <- opt$latmin
LATMAX <- opt$latmax
LONMIN <- opt$lonmin
LONMAX <- opt$lonmax

MINYEAR <- opt$minyear
EXTINCT <- opt$noextinct
RESOLUTION <- opt$resolution
CPUTHREADS <- opt$threads
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input occurrences: ", INPUT, "\n", sep=""))

cat(paste("Selected phyla: ",    PHYLUM, "\n", sep = ""))
cat(paste("Selected classes: ",  CLASS, "\n", sep = ""))
cat(paste("Selected orders: ",   ORDER, "\n", sep = ""))
cat(paste("Selected families: ", FAMILY, "\n", sep = ""))

cat(paste("Country codes: ",     COUNTRY, "\n", sep = ""))
cat(paste("Minimum latitude: ",  LATMIN, "\n", sep = ""))
cat(paste("Maximum latitude: ",  LATMAX, "\n", sep = ""))
cat(paste("Minimum longitude: ", LONMIN, "\n", sep = ""))
cat(paste("Maximum longitude: ", LONMAX, "\n", sep = ""))

cat(paste("Minimum year of occurrence: ", MINYEAR, "\n", sep=""))
cat(paste("List of extict species: ", EXTINCT, "\n", sep=""))
cat(paste("Spatial resolution: ", RESOLUTION, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output file prefix: ", OUTPUT, "\n", sep=""))
cat("\n")


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

############################################## Main pipeline


## Open dataset
ds <- arrow::open_dataset(INPUT)

## General filtering pipeline
## Based on scipts by John Waller
## https://data-blog.gbif.org/post/gbif-filtering-guide/
dsf <- ds %>%
  select(-mediatype,-issue) %>%
  filter(!is.na(species)) %>%
  filter(taxonrank %in% c("SPECIES", "SUBSPECIES", "VARIETY", "FORM")) %>%
  filter(occurrencestatus == "PRESENT") %>%
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%
  filter(!establishmentmeans %in% c("MANAGED", "INTRODUCED", "INVASIVE", "NATURALISED")) %>%
  filter(!is.na(decimallongitude)) %>% 
  filter(!is.na(decimallatitude)) %>% 
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>%
  filter(decimallatitude != decimallongitude) %>%
  filter(year >= MINYEAR) %>% 
  filter(coordinateprecision < 0.1 | is.na(coordinateprecision)) %>% 
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301, 3036, 999, 9999))



## Taxonomy filters
if(is.na(PHYLUM)){
  PHYLUM <- strsplit(x = PHYLUM, split = ",")[[1]]  # split multiple records
  dsf <- dsf %>% filter(phylum %in% PHYLUM)
}

if(is.na(CLASS)){
  CLASS <- strsplit(x = CLASS, split = ",")[[1]]
  dsf <- dsf %>% filter(class %in% CLASS)
}

if(is.na(ORDER)){
  ORDER <- strsplit(x = ORDER, split = ",")[[1]]
  dsf <- dsf %>% filter(order %in% ORDER)
}

if(is.na(FAMILY)){
  FAMILY <- strsplit(x = FAMILY, split = ",")[[1]]
  dsf <- dsf %>% filter(family %in% FAMILY)
}


## Spatial filters
if(is.na(COUNTRY)){
  COUNTRY <- strsplit(x = COUNTRY, split = ",")[[1]]
  dsf <- dsf %>% filter(countrycode %in% COUNTRY)
}

if(is.na(LATMIN)){
  dsf <- dsf %>% filter(decimallatitude >= LATMIN)
}
if(is.na(LATMAX)){
  dsf <- dsf %>% filter(decimallatitude <= LATMAX)
}
if(is.na(LONMIN)){
  dsf <- dsf %>% filter(decimallongitude >= LONMIN)
}
if(is.na(LONMAX)){
  dsf <- dsf %>% filter(decimallongitude <= LONMAX)
}

