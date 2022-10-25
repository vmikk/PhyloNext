#!/usr/bin/env Rscript

## Filter GBIF occurrence data

cat("Filtering GBIF species occurrence data\n")

## Input data could be located in the cloud.
## E.g., to query GBIF AWS snapshot use: 
##  --input "s3://gbif-open-data-eu-central-1/occurrence/2022-01-01/occurrence.parquet"

## Usage:
# ./10_Filter_occurrences.R \
#    --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet" \
#    --family "Fabaceae" \
#    --country "AU" \
#    --latmin -55.3228175 \
#    --latmax -9.0882278 \
#    --lonmin 72.2460938 \
#    --lonmax 168.2249543 \
#    --minyear 1945 \
#    --roundcoords 2 \
#    --threads 10 \
#    --noccurrences 30 \
#    --output "Fabaceae_in_AU"

## Note about rounding of coordinates (`--roundcoords 2`):
## A value in decimal degrees to an accuracy of 2 decimal places is accurate to ~1.11 km at the equator

# TO DO:
# - add country name validation
# - check if the input folder exists (except of 's3')

############################################## Parse input parameters

## Check time
start_time <- Sys.time()


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
  make_option("--genus", action="store", default=NA, type='character', help="Comma-separated list of genera to select"),
  make_option("--specieskeys", action="store", default=NA, type='character', help="File with user-supplied GBIF specieskeys"),
  
  ## Spatial filters
  make_option("--country", action="store", default=NA, type='character', help="Comma-separated list of country codes (e.g., AU,CA)"),
  make_option("--latmin", action="store", default=NA, type='double', help="Minimum latitude"),
  make_option("--latmax", action="store", default=NA, type='double', help="Maximum latitude"),
  make_option("--lonmin", action="store", default=NA, type='double', help="Minimum longitude"),
  make_option("--lonmax", action="store", default=NA, type='double', help="Maximum longitude"),

  ## Additional filters
  make_option("--minyear", action="store", default=1945, type='integer', help="Minimum year of occurrence (default, 1945)"),
  make_option("--noextinct", action="store", default=NA, type='character', help="Remove extinct species (provide a file with extinct specieskeys)"),
  make_option("--excludehuman", action="store", default=TRUE, type='logical', help="Exclude human records (genus Homo)"),

  make_option(c("--roundcoords"), action="store", default=2L, type='integer', help="Round spatial coordinates to the N decimal places, to reduce the dataset size (default, 2). To disable, set to a negative value."),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-n", "--noccurrences"), action="store", default=30, type='double', help="Occurrence threshold (used for parquet partitioning)"),
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

## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% "NA"){ x <- NA }
  return(x)
}

## Assign variables
INPUT <- opt$input

PHYLUM <- to_na( opt$phylum )
CLASS  <- to_na( opt$class )
ORDER  <- to_na( opt$order )
FAMILY <- to_na( opt$family )
GENUS  <- to_na( opt$genus )
SPECIESKEYS  <- to_na( opt$specieskeys )

COUNTRY <- to_na( opt$country )
LATMIN <- as.numeric( to_na(opt$latmin) )
LATMAX <- as.numeric( to_na(opt$latmax) )
LONMIN <- as.numeric( to_na(opt$lonmin) )
LONMAX <- as.numeric( to_na(opt$lonmax) )

MINYEAR <- as.numeric(to_na( opt$minyear) )
EXTINCT <- to_na( opt$noextinct)
EXCLUDEHUMAN <- as.logical( opt$excludehuman )
ROUNDCOORDS <- as.numeric( opt$roundcoords )
CPUTHREADS   <- as.numeric(opt$threads)
OCCURRENCES  <- as.numeric(opt$noccurrences)
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input occurrences: ", INPUT, "\n", sep=""))

cat(paste("Selected phyla: ",    PHYLUM, "\n", sep = ""))
cat(paste("Selected classes: ",  CLASS,  "\n", sep = ""))
cat(paste("Selected orders: ",   ORDER,  "\n", sep = ""))
cat(paste("Selected families: ", FAMILY, "\n", sep = ""))
cat(paste("Selected genera: ",   GENUS,  "\n", sep = ""))
cat(paste("File with GBIF specieskeys: ", SPECIESKEYS,  "\n", sep = ""))

cat(paste("Country codes: ",     COUNTRY, "\n", sep = ""))
cat(paste("Minimum latitude: ",  LATMIN,  "\n", sep = ""))
cat(paste("Maximum latitude: ",  LATMAX,  "\n", sep = ""))
cat(paste("Minimum longitude: ", LONMIN,  "\n", sep = ""))
cat(paste("Maximum longitude: ", LONMAX,  "\n", sep = ""))

cat(paste("Minimum year of occurrence: ", MINYEAR, "\n", sep=""))
cat(paste("List of extict species: ", EXTINCT, "\n", sep=""))
cat(paste("Exclusion of human records: ", EXCLUDEHUMAN, "\n", sep=""))
cat(paste("Round coordinates: ", ROUNDCOORDS, "\n", sep=""))

cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Occurrence threshold for parquet partitioning: ", OCCURRENCES, "\n", sep=""))
cat(paste("Output prefix: ", OUTPUT, "\n", sep=""))
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

cat("\n")

## Set CPU thread pool
cat("Number of available CPU threads: ", cpu_count(), "\n")
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
set_cpu_count(CPUTHREADS)           # for libarrow
setDTthreads(threads = CPUTHREADS)  # for data.table


## Load extinct species list
if(!is.na(EXTINCT)){
  EXTINCT <- fread(file = EXTINCT, sep = "\t")
  colnames(EXTINCT) <- "specieskey"
  cat("Extinct species list loaded. Number of records: ", nrow(EXTINCT), "\n")
}

## Load user-supplied specieskeys
if(!is.na(SPECIESKEYS)){
  SPECIESKEYS <- fread(file = SPECIESKEYS, sep = "\t")
  colnames(SPECIESKEYS) <- "specieskey"
  SPECIESKEYS <- unique(na.omit(SPECIESKEYS))
  cat("Specieskey list loaded. Number of records: ", nrow(SPECIESKEYS), "\n")
}


############################################## Main pipeline


## Open dataset
cat("Loading Parquet data\n")
ds <- arrow::open_dataset(INPUT)

## General filtering pipeline
## Based on scipts by John Waller
## https://data-blog.gbif.org/post/gbif-filtering-guide/
cat("General data filteing:\n")
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
  filter(coordinateprecision < 0.1 | is.na(coordinateprecision)) %>% 
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301, 3036, 999, 9999))

## Year
if(!is.na(MINYEAR)){
  cat("..Filtering by collection date\n")
  dsf <- dsf %>% filter(year >= MINYEAR)
}

## Taxonomy filters
if(!is.na(PHYLUM)){
  cat("..Filtering by Phylum\n")
  PHYLUM <- strsplit(x = PHYLUM, split = ",")[[1]]  # split multiple records
  dsf <- dsf %>% filter(phylum %in% PHYLUM)
}

if(!is.na(CLASS)){
  cat("..Filtering by Class\n")
  CLASS <- strsplit(x = CLASS, split = ",")[[1]]
  dsf <- dsf %>% filter(class %in% CLASS)
}

if(!is.na(ORDER)){
  cat("..Filtering by Order\n")
  ORDER <- strsplit(x = ORDER, split = ",")[[1]]
  dsf <- dsf %>% filter(order %in% ORDER)
}

if(!is.na(FAMILY)){
  cat("..Filtering by Family\n")
  FAMILY <- strsplit(x = FAMILY, split = ",")[[1]]
  dsf <- dsf %>% filter(family %in% FAMILY)
}

if(!is.na(GENUS)){
  cat("..Filtering by Genus\n")
  GENUS <- strsplit(x = GENUS, split = ",")[[1]]
  dsf <- dsf %>% filter(genus %in% GENUS)
}

## Custom specieskeys
if(!is.na(SPECIESKEYS[[1]][1])){
  cat("..Filtering by specieskeys\n")
  dsf <- dsf %>% filter(specieskey %in% SPECIESKEYS$specieskey)
}


## Spatial filters
if(!is.na(COUNTRY)){
  cat("..Filtering by Country\n")
  COUNTRY <- strsplit(x = COUNTRY, split = ",")[[1]]
  dsf <- dsf %>% filter(countrycode %in% COUNTRY)
}

if(!is.na(LATMIN)){
  cat("..Filtering by min latitude\n")
  dsf <- dsf %>% filter(decimallatitude >= LATMIN)
}
if(!is.na(LATMAX)){
  cat("..Filtering by max latitude\n")
  dsf <- dsf %>% filter(decimallatitude <= LATMAX)
}
if(!is.na(LONMIN)){
  cat("..Filtering by min longitude\n")
  dsf <- dsf %>% filter(decimallongitude >= LONMIN)
}
if(!is.na(LONMAX)){
  cat("..Filtering by max longitude\n")
  dsf <- dsf %>% filter(decimallongitude <= LONMAX)
}


## Remove extinct species
if(!is.na(EXTINCT[[1]][1])){
  cat("..Filtering out extinct species\n")
  dsf <- dsf %>% filter(!specieskey %in% EXTINCT$specieskey)
}

## Removal of human records
if(EXCLUDEHUMAN == TRUE){
  cat("..Excluding human records (genus Homo)\n")
  dsf <- dsf %>% filter(!genus %in% "Homo")
}


## Round coordiantes, to reduce the dataset size
cat("..Rounding coordinates\n")
if(ROUNDCOORDS >= 0){
  cat("..Rounding coordinates\n")
  dsf <- dsf %>%
    mutate(
      decimallatitude  = round(decimallatitude,  ROUNDCOORDS),
      decimallongitude = round(decimallongitude, ROUNDCOORDS))
}

## Select columns and remove duplicated records
cat("..Column selection\n")
dsf <- dsf %>%
  select(specieskey, species, decimallongitude, decimallatitude) %>%
  distinct()



## Count number of records by species
cat("Counting number of records per species\n")
sp_counts <- dsf %>%
  count(specieskey) %>% 
  collect() %>% 
  mutate(Partition = case_when(n <= OCCURRENCES ~ "low", 
                              n > OCCURRENCES ~ "high"))

smr <- table(sp_counts$Partition)
if("low" %in% names(smr)){
  cat("..Number of species with low number of occurrences: ", smr[["low"]], "\n")
}
if("high" %in% names(smr)){
  cat("..Number of species with high number of occurrences: ", smr[["high"]], "\n")
}

## Export species counts
fwrite(x = sp_counts,
  file = "SpeciesCounts.txt",
  sep = "\t")


## Group by species  -- there would be too many partitions (default arrow max = 1024)
# cat("Groupping by species\n")
# dsf <- dsf %>%
#   group_by(specieskey)


## Add partition ID to the dataset
dsf <- dsf %>% 
  left_join(sp_counts[, c("specieskey", "Partition")]) %>%
  group_by(Partition)


## Export species occurrence, partition by species
cat("Exporting filtered occurrence data\n")

# outdir <- paste0(OUTPUT, ".parquet")
outdir <- OUTPUT
if(! outdir %in% "."){ dir.create(outdir, showWarnings = FALSE) }

dir_high <- paste0(outdir, "/Partition=high")
dir_low <- paste0(outdir, "/Partition=low")
dir.create(dir_high, showWarnings = FALSE)
dir.create(dir_low, showWarnings = FALSE)

write_dataset(
  dataset = dsf,
  path = outdir,
  format = "parquet",
  basename_template = "occ_{i}")




##################### Session info

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
