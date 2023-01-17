#!/usr/bin/env Rscript

## Derived datasets are citable records of GBIF-mediated occurrence data

cat("Preparation of derived dataset\n")
cat("Script name: 16_Derived_dataset.R\n")

## Usage example:
# ./16_Derived_dataset.R \
#   --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet" \
#    --family "Fabaceae" \
#    --country "AU" \
#    --latmin -55.3228175 \
#    --latmax -9.0882278 \
#    --lonmin 72.2460938 \
#    --lonmax 168.2249543 \
#    --minyear 1945 \
#    --roundcoords 2 \
#    --terrestrial        "pipeline_data/Land_Buffered_025_dgr.RData" \
#    --rmcountrycentroids "pipeline_data/CC_CountryCentroids_buf_1000m.RData" \
#    --rmcountrycapitals  "pipeline_data/CC_Capitals_buf_10000m.RData" \
#    --rminstitutions     "pipeline_data/CC_Institutions_buf_100m.RData" \
#    --rmurban            "pipeline_data/CC_Urban.RData" \
#    --speciestree        "02.Biodiverse_input/Trimmed_occurrences.csv" \
#    --threads 10 \
#    --output "Dataset_DOIs.txt"

## TO DO:
# - run `rgbif::gbif_citation` in parallel


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(

  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the directory with GBIF snapshot of occurrence records in Parquet format"),
  make_option(c("-o", "--output"), action="store", default="Dataset_DOIs.txt", type='character', help="Output file name"),

  ## Taxonomy filters
  make_option("--phylum", action="store", default=NA, type='character', help="Comma-separated list of phyla to select"),
  make_option("--class", action="store", default=NA, type='character', help="Comma-separated list of classes to select"),
  make_option("--order", action="store", default=NA, type='character', help="Comma-separated list of orders to select"),
  make_option("--family", action="store", default=NA, type='character', help="Comma-separated list of families to select"),
  make_option("--genus", action="store", default=NA, type='character', help="Comma-separated list of genera to select"),
  make_option("--specieskeys", action="store", default=NA, type='character', help="File with user-supplied GBIF specieskeys"),
  make_option(c("-s", "--speciestree"), action="store", default=NA, type='character', help="File with specieskeys found in the phylogenetic tree (output of `12_Prepare_Biodiverse_input.R`)"),
  
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
  make_option("--basisofrecordinclude",  action="store", default=NA, type='character', help="Basis of record to include from the data"),
  make_option("--basisofrecordexclude", action="store", default="FOSSIL_SPECIMEN,LIVING_SPECIMEN", type='character', help="Basis of record to exclude from the data"),

  ## Coordinate precision and uncertainty filters
  make_option("--coordprecision", action="store", default=0.1, type='double', help="Coordinate precision threshold (max allowed value)"),
  make_option("--coorduncertainty", action="store", default=10000, type='double', help="Maximum allowed coordinate uncertainty, meters"),
  make_option("--coorduncertaintyexclude", action="store", default="301,3036,999,9999", type='character', help="Black-listed values of coordinate uncertainty"),

  ## Spatial filters (shapefile-based)
  make_option(c("-l", "--terrestrial"), action="store", default=NA, type='character', help="Remove non-terrestrial occurrences, provide land polygon in sf-format"),
  make_option(c("-c", "--rmcountrycentroids"), action="store", default=NA, type='character', help="Remove records within a radius around the geographic centroids of political countries and provinces"),
  make_option(c("-k", "--rmcountrycapitals"), action="store", default=NA, type='character', help="Remove records within a radius around country capitals"),
  make_option(c("-b", "--rminstitutions"), action="store", default=NA, type='character', help="Remove records in the vicinity of biodiversity institutions"),
  make_option(c("-u", "--rmurban"), action="store", default=NA, type='character', help="Remove records inside urban areas"),

  ## Spatial aggregation
  make_option(c("-r", "--resolution"), action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),
  make_option(c("--roundcoords"), action="store", default=2L, type='integer', help="Round spatial coordinates to the N decimal places, to reduce the dataset size (default, 2). To disable, set to a negative value."),

  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))



## Validation of the required argiments
if(is.na(opt$input)){
  cat("Input is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  stop("Output file is not specified.\n")
}


## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
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
SPECIESTREE  <- to_na( opt$speciestree )

COORDPREC      <- as.numeric( to_na(opt$coordprecision) )
COORDUNCRTMAX  <- as.numeric( to_na(opt$coorduncertainty) )
COORDUNCRTEXCL <- to_na(opt$coorduncertaintyexclude)

COUNTRY <- to_na( opt$country )
LATMIN  <- as.numeric( to_na(opt$latmin) )
LATMAX  <- as.numeric( to_na(opt$latmax) )
LONMIN  <- as.numeric( to_na(opt$lonmin) )
LONMAX  <- as.numeric( to_na(opt$lonmax) )

MINYEAR <- as.numeric(to_na( opt$minyear) )
EXTINCT <- to_na( opt$noextinct)
EXCLUDEHUMAN <- as.logical( opt$excludehuman )

BASISINCL <- to_na( opt$basisofrecordinclude )
BASISEXCL <- to_na( opt$basisofrecordexclude )

TERRESTRIAL <- to_na( opt$terrestrial )
CC_COUNTRY  <- to_na( opt$rmcountrycentroids )
CC_CAPITAL  <- to_na( opt$rmcountrycapitals )
CC_INSTIT   <- to_na( opt$rminstitutions )
CC_URBAN    <- to_na( opt$rmurban )

RESOLUTION  <- as.integer(opt$resolution)
ROUNDCOORDS <- as.numeric( opt$roundcoords )

CPUTHREADS <- as.numeric(opt$threads)
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input occurrences: ", INPUT, "\n", sep=""))

cat(paste("Selected phyla: ",    PHYLUM, "\n", sep = ""))
cat(paste("Selected classes: ",  CLASS,  "\n", sep = ""))
cat(paste("Selected orders: ",   ORDER,  "\n", sep = ""))
cat(paste("Selected families: ", FAMILY, "\n", sep = ""))
cat(paste("Selected genera: ",   GENUS,  "\n", sep = ""))
cat(paste("File with GBIF specieskeys: ", SPECIESKEYS,  "\n", sep = ""))
cat(paste("File with species in phylogenetic tree: ", SPECIESTREE, "\n", sep=""))

cat(paste("Coordinate precision threshold: ",                COORDPREC,      "\n", sep = ""))
cat(paste("Maximum allowed coordinate uncertainty: ",        COORDUNCRTMAX,  "\n", sep = ""))
cat(paste("Black-listed values of coordinate uncertainty: ", COORDUNCRTEXCL, "\n", sep = ""))

cat(paste("Country codes: ",     COUNTRY, "\n", sep = ""))
cat(paste("Minimum latitude: ",  LATMIN,  "\n", sep = ""))
cat(paste("Maximum latitude: ",  LATMAX,  "\n", sep = ""))
cat(paste("Minimum longitude: ", LONMIN,  "\n", sep = ""))
cat(paste("Maximum longitude: ", LONMAX,  "\n", sep = ""))

cat(paste("Basis of record to include: ", BASISINCL, "\n", sep=""))
cat(paste("Basis of record to exclude: ", BASISEXCL, "\n", sep=""))
cat(paste("Minimum year of occurrence: ", MINYEAR, "\n", sep=""))
cat(paste("List of extict species: ",     EXTINCT, "\n", sep=""))
cat(paste("Exclusion of human records: ", EXCLUDEHUMAN, "\n", sep=""))
cat(paste("Round coordinates: ",          ROUNDCOORDS, "\n", sep=""))

cat(paste("Terrestrial data: ",  TERRESTRIAL, "\n", sep=""))
cat(paste("Country and province centroids: ", CC_COUNTRY, "\n", sep=""))
cat(paste("Capitals: ",     CC_CAPITAL, "\n", sep=""))
cat(paste("Institutions: ", CC_INSTIT, "\n", sep=""))
cat(paste("Uraban areas: ", CC_URBAN, "\n", sep=""))

cat(paste("Spatial resolution: ",  RESOLUTION, "\n", sep=""))
cat(paste("Coordinate rounding: ", ROUNDCOORDS, "\n", sep=""))

cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output file: ", OUTPUT, "\n", sep=""))
cat("\n")

############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}


load_pckg("arrow")
load_pckg("data.table")
load_pckg("dplyr")
load_pckg("h3")
load_pckg("sf")
load_pckg("rgbif")
# load_pckg("plyr")


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

## Load specieskeys for taxa found in phylogenetic tree
if(!is.na(SPECIESTREE)){
  SPECIESTREE <- fread(file = SPECIESTREE, sep = ",")
  cat("Specieskeys in phylogenetic tree: ", length(unique(SPECIESTREE$specieskey)), "\n")
}


## Function to suppress function output
## (will be used for silencing of garbage collector)
quiet <- function(x) { 
  sink("/dev/null")
  on.exit(sink()) 
  invisible(force(x)) 
}


## Parameters for debugging
# INPUT <- "/mnt/Dat2/GBIF/Parquet/2022-08-01/"
# OUTPUT <- "Dataset_DOIs.txt"
# COORDPREC <- 0.1
# COORDUNCRTMAX <- 10000
# COORDUNCRTEXCL <- "301,3036,999,9999"
# PHYLUM <- NA
# CLASS <- NA
# ORDER <- NA
# FAMILY <- "Fabaceae"
# GENUS <- "Acacia"
# SPECIESKEYS <- NA
# COUNTRY <- "AU"
# LATMIN <- NA
# LATMAX <- NA
# LONMIN <- NA
# LONMAX <- NA
# MINYEAR <- 2005
# EXTINCT <- NA
# EXCLUDEHUMAN <- TRUE
# BASISINCL <- NA
# BASISEXCL <- "FOSSIL_SPECIMEN,LIVING_SPECIMEN"
# TERRESTRIAL <- " ~/.nextflow/assets/vmikk/phylonext/pipeline_data/Land_Buffered_025_dgr.RData"
# CC_COUNTRY  <- " ~/.nextflow/assets/vmikk/phylonext/pipeline_data/CC_CountryCentroids_buf_1000m.RData"
# CC_CAPITAL  <- " ~/.nextflow/assets/vmikk/phylonext/pipeline_data/CC_Capitals_buf_10000m.RData"
# CC_INSTIT   <- " ~/.nextflow/assets/vmikk/phylonext/pipeline_data/CC_Institutions_buf_100m.RData"
# CC_URBAN    <- " ~/.nextflow/assets/vmikk/phylonext/pipeline_data/CC_Urban.RData"
# RESOLUTION <- 4L
# ROUNDCOORDS <- 2L
# set_cpu_count(10)
# setDTthreads(threads = 10)


############################################## Intial filtering of Parquet data

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
  filter(!establishmentmeans %in% c("MANAGED", "INTRODUCED", "INVASIVE", "NATURALISED")) %>%
  filter(!is.na(decimallongitude)) %>% 
  filter(!is.na(decimallatitude)) %>% 
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>%
  filter(decimallatitude != decimallongitude)

## Coordinate precision filter
if(!is.na(COORDPREC)){
  cat("..Filtering by coordinate precision\n")
  dsf <- dsf %>% filter(coordinateprecision < COORDPREC | is.na(coordinateprecision))
}

## Coordinate uncertainty filter
if(!is.na(COORDUNCRTMAX)){
  cat("..Filtering by coordinate uncertainty\n")
  dsf <- dsf %>% filter(coordinateuncertaintyinmeters < COORDUNCRTMAX | is.na(coordinateuncertaintyinmeters))
}

## Coordinate uncertainty blacklist filter
if(!is.na(COORDUNCRTEXCL)){
  cat("..Filtering by coordinate uncertainty black-listed values\n")
  COORDUNCRTEXCL <- strsplit(x = COORDUNCRTEXCL, split = ",")[[1]]
  COORDUNCRTEXCL <- as.numeric(COORDUNCRTEXCL)
  dsf <- dsf %>% filter(!coordinateuncertaintyinmeters %in% COORDUNCRTEXCL)
}

## Basis of record filters
if(!is.na(BASISINCL) & !is.na(BASISEXCL)){
  cat("..Filtering by basis of record (inclusion and exclusion)\n")

  BASISINCL <- strsplit(x = BASISINCL, split = ",")[[1]]
  BASISEXCL <- strsplit(x = BASISEXCL, split = ",")[[1]]

  ## Check if selected values are not mutually exclusive
  if(length(intersect(BASISINCL, BASISEXCL)) > 0){
    stop("Mutually exclusive basis of record selected!\n")
  }

  dsf <- dsf %>% 
    filter( (!basisofrecord %in% BASISEXCL) & (basisofrecord %in% BASISINCL) )

} else if (!is.na(BASISINCL) & is.na(BASISEXCL)){
  cat("..Filtering by basis of record (inclusion only)\n")

  BASISINCL <- strsplit(x = BASISINCL, split = ",")[[1]]
  dsf <- dsf %>% filter( basisofrecord %in% BASISINCL )

} else if (is.na(BASISINCL) & ! is.na(BASISEXCL)){
  cat("..Filtering by basis of record (exclusion only)\n")

  BASISEXCL <- strsplit(x = BASISEXCL, split = ",")[[1]]
  dsf <- dsf %>% filter( ! basisofrecord %in% BASISEXCL )

} else {
  cat("..No filtering based on `Basis of record` field\n")
}


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
if(!is.na(SPECIESKEYS[[1]][1]) & is.na(SPECIESTREE[[1]][1])){
  cat("..Filtering by specieskeys\n")
  dsf <- dsf %>% filter(specieskey %in% SPECIESKEYS$specieskey)
}

## Specieskeys in phylogenetic tree
if(!is.na(SPECIESTREE[[1]][1])){
  cat("..Filtering by specieskeys in phylogenetic tree\n")
  dsf <- dsf %>% filter(specieskey %in% unique(SPECIESTREE$specieskey))
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
if(ROUNDCOORDS >= 0){
  cat("..Rounding coordinates\n")
  dsf <- dsf %>%
    mutate(
      decimallatitude  = round(decimallatitude,  ROUNDCOORDS),
      decimallongitude = round(decimallongitude, ROUNDCOORDS))
}


## Select columns and remove duplicated records
cat("..Column selection and unique record counts\n")
dsf <- dsf %>%
  select(decimallongitude, decimallatitude, datasetkey) %>%
  distinct()

## currently, arrow does not support mutating a column withour pulling data into R
## e.g.,
#  group_by(decimallongitude, decimallatitude) %>%
#  mutate(datasetkey = paste(sort(unique(datasetkey)), collapse = ";"))

## Collect the data
cat("..Collecting data\n")
datt <- dsf %>% collect()

## Convert to data.table
setDT(datt)

cat("Data table created\n")
cat("..Number of unique datasets (before spatial filtering): ",
  length(unique(datt$datasetkey)),
  "\n")


cat("Aggregating dataset keys by coordinates\n")
cat("..Number of records before aggregation: ", nrow(datt), "\n")

datt <- datt[ , .(
  datasetkey = paste(sort(unique(datasetkey)), collapse = ";")
  ), by = c("decimallongitude", "decimallatitude")]

cat("..Number of records after aggregation: ", nrow(datt), "\n")


################################################# Spatial filtering

cat("Spatial filtering\n")


## Remove non-terrestrial records 
if(!is.na(TERRESTRIAL)){
  cat("..Removing non-terrestrial records\n")
  
  ## Load land mask
  cat("...Loading land mask\n")
  TERRESTRIAL <- readRDS(TERRESTRIAL)

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points are on land
  cat("...Intersecting polygons and data points\n")
  land_intersect <- st_intersects(pts, TERRESTRIAL)
  land <- lengths(land_intersect) > 0

  non_terr <- sum(!land)
  cat("...Number of non-terrestrial points: ", non_terr, "\n")

  ## Remove outliers
  if(non_terr > 0){
    removed_nonterrestrial <- datt[ ! land ]   # outliers
    datt <- datt[ land ]
  } else {
    removed_nonterrestrial <- NULL   # no non-terrestrial samples found
  }

  rm(pts, TERRESTRIAL)
  quiet( gc() )
} else {
  removed_nonterrestrial <- NULL     # no terrestrial filtering was performed
}


## Remove country centroids and province centroids (similar to CoordinateCleaner::cc_cen)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_cen.R
## Default buffer = 1 km
if(!is.na(CC_COUNTRY)){
  cat("..Removing occurrences inside country and province centroids\n")

  ## Load polygons
  cat("...Loading polygons\n")
  CC_COUNTRY <- readRDS(CC_COUNTRY)
  cat("...Number of polygons provided: ", nrow(CC_COUNTRY), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("...Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_COUNTRY)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("...Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("...Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_COUNTRY <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_COUNTRY <- NULL   # no found
  }

  rm(pts, CC_COUNTRY)
  quiet( gc() )
} else {
  removed_CC_COUNTRY <- NULL     # no filtering was performed
}


## Filter country capitals (similar to CoordinateCleaner::cc_cap)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_cap.R
## Default buffer = 10 km
if(!is.na(CC_CAPITAL)){
  cat("..Removing occurrences within capitals\n")

  ## Load polygons
  cat("...Loading polygons\n")
  CC_CAPITAL <- readRDS(CC_CAPITAL)
  cat("...Number of polygons provided: ", nrow(CC_CAPITAL), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("...Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_CAPITAL)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("...Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("...Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_CAPITAL <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_CAPITAL <- NULL   # no found
  }

  rm(pts, CC_CAPITAL)
  quiet( gc() )
} else {
  removed_CC_CAPITAL <- NULL     # no filtering was performed
}


## Remove records in the vicinity of biodiversity institutions (similar to CoordinateCleaner::cc_inst)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_inst.R
## Default buffer = 100 m
if(!is.na(CC_INSTIT)){
  cat("..Removing occurrences in the vicinity of biodiversity institutions\n")

  ## Load polygons
  cat("...Loading polygons\n")
  CC_INSTIT <- readRDS(CC_INSTIT)
  cat("...Number of polygons provided: ", nrow(CC_INSTIT), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("...Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_INSTIT)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("...Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("...Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_INSTIT <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_INSTIT <- NULL   # no found
  }

  rm(pts, CC_INSTIT)
  quiet( gc() )
} else {
  removed_CC_INSTIT <- NULL     # no filtering was performed
}


## Remove records inside urban areas (similar to CoordinateCleaner::cc_urb)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_urb.R
if(!is.na(CC_URBAN)){
  cat("..Removing occurrences inside urban areas\n")

  ## Load polygons
  cat("...Loading polygons\n")
  CC_URBAN <- readRDS(CC_URBAN)
  cat("...Number of polygons provided: ", nrow(CC_URBAN), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("...Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_URBAN)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("...Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("...Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_URBAN <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_URBAN <- NULL   # no found
  }

  rm(pts, CC_URBAN)
  quiet( gc() )
} else {
  removed_CC_URBAN <- NULL     # no filtering was performed
}

cat("Data filtering finished.\n")
cat("Extracting unique dataset keys\n")

dataset_keys <- unique(datt$datasetkey)

## Split combined datasets
dataset_keys <- plyr::alply(.data = dataset_keys,
  .margins = 1, .fun = function(x){
    strsplit(x = x, split = ";")[[1]]
  })

dataset_keys <- unique( unlist(dataset_keys) )

cat(".. Number of unique dataset keys (after filtering): ",
  length(dataset_keys),
  "\n")


cat("Parsing dataset info with `rgbif` package\n")

dataset_dois <- plyr::alply(.data = dataset_keys,
  .margins = 1, .fun = function(x){
    cit <- try( rgbif::gbif_citation(x = x) )
    if(! "try-error" %in% class(cit)){
      res <- data.table(
        datasetkey = x,
        Title      = cit$citation$title,
        Text       = cit$citation$text,
        License    = cit$rights
        )
    } else {
      res <- data.table(
        datasetkey = x,
        Title      = NA,
        Text       = NA,
        License    = NA
        )
    }
    return(res)
  })

dataset_dois <- rbindlist(dataset_dois)
setorder(x = dataset_dois, datasetkey)

if(any(is.na(dataset_dois$Title))){
  cat("WARNING: information is missing for some of the datasets!\n")
  cat(".. Number of datasets with missing information: ",
    sum(is.na(dataset_dois$Title)),
    "\n")
}


cat("Exporting results\n")

fwrite(x = dataset_dois, file = OUTPUT, sep = "\t")


#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
