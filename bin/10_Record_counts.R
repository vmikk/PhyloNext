#!/usr/bin/env Rscript

## Count number of records per gridcell in the GBIF occurrence data

cat("Counting the number of records per grid cell\n")

## Usage:
# ./10_Record_counts.R \
#    --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet" \
#    --family "Fabaceae" \
#    --country "AU" \
#    --latmin -55.3228175 \
#    --latmax -9.0882278 \
#    --lonmin 72.2460938 \
#    --lonmax 168.2249543 \
#    --minyear 1945 \
#    --roundcoords 2 \
#    --terrestrial "pipeline_data/Land_Buffered_025_dgr.RData" \
#    --rmcountrycentroids "pipeline_data/CC_CountryCentroids_buf_1000m.RData" \
#    --rmcountrycapitals "pipeline_data/CC_Capitals_buf_10000m.RData" \
#    --rminstitutions "pipeline_data/CC_Institutions_buf_100m.RData" \
#    --rmurban "pipeline_data/CC_Urban.RData" \
#    --threads 10 \
#    --output "Record_counts"

## Note about rounding of coordinates:
## A value in decimal degrees to an accuracy of 
##   1 decimal place  is accurate to ~11.  km at the equator
##   2 decimal places is accurate to ~1.11 km at the equator
##   3 decimal places is accurate to ~111  m  at the equator

## Outputs:
# Record_counts_H3.RData
# Record_counts_H3.txt.gz
# Record_counts_Outliers.txt.gz


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

  ## Spatial filters (shapefile-based)
  make_option(c("-l", "--terrestrial"), action="store", default=NA, type='character', help="Remove non-terrestrial occurrences, provide land polygon in sf-format"),
  make_option(c("-c", "--rmcountrycentroids"), action="store", default=NA, type='character', help="Remove records within a radius around the geographic centroids of political countries and provinces"),
  make_option(c("-k", "--rmcountrycapitals"), action="store", default=NA, type='character', help="Remove records within a radius around country capitals"),
  make_option(c("-b", "--rminstitutions"), action="store", default=NA, type='character', help="Remove records in the vicinity of biodiversity institutions"),
  make_option(c("-u", "--rmurban"), action="store", default=NA, type='character', help="Remove records inside urban areas"),

  ## Spatial aggregation
  make_option(c("-r", "--resolution"), action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),
  make_option(c("--roundcoords"), action="store", default=2L, type='integer', help="Round spatial coordinates to the N decimal places, to reduce the dataset size (default, 2). To disable, set to a negative value."),

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
GENUS <- to_na( opt$genus )
SPECIESKEYS  <- to_na( opt$specieskeys )

COUNTRY <- to_na( opt$country )
LATMIN <- as.numeric( to_na(opt$latmin) )
LATMAX <- as.numeric( to_na(opt$latmax) )
LONMIN <- as.numeric( to_na(opt$lonmin) )
LONMAX <- as.numeric( to_na(opt$lonmax) )

MINYEAR <- as.numeric(to_na( opt$minyear) )
EXTINCT <- to_na( opt$noextinct)
EXCLUDEHUMAN <- as.logical( opt$excludehuman )

TERRESTRIAL <- to_na( opt$terrestrial )
CC_COUNTRY <- to_na( opt$rmcountrycentroids )
CC_CAPITAL <- to_na( opt$rmcountrycapitals )
CC_INSTIT  <- to_na( opt$rminstitutions )
CC_URBAN   <- to_na( opt$rmurban )

RESOLUTION <- as.integer(opt$resolution)
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

cat(paste("Country codes: ",     COUNTRY, "\n", sep = ""))
cat(paste("Minimum latitude: ",  LATMIN,  "\n", sep = ""))
cat(paste("Maximum latitude: ",  LATMAX,  "\n", sep = ""))
cat(paste("Minimum longitude: ", LONMIN,  "\n", sep = ""))
cat(paste("Maximum longitude: ", LONMAX,  "\n", sep = ""))

cat(paste("Minimum year of occurrence: ", MINYEAR, "\n", sep=""))
cat(paste("List of extict species: ", EXTINCT, "\n", sep=""))
cat(paste("Exclusion of human records: ", EXCLUDEHUMAN, "\n", sep=""))
cat(paste("Round coordinates: ", ROUNDCOORDS, "\n", sep=""))

cat(paste("Terrestrial data: ", TERRESTRIAL, "\n", sep=""))
cat(paste("Country and province centroids: ", CC_COUNTRY, "\n", sep=""))
cat(paste("Capitals: ", CC_CAPITAL, "\n", sep=""))
cat(paste("Institutions: ", CC_INSTIT, "\n", sep=""))
cat(paste("Uraban areas: ", CC_URBAN, "\n", sep=""))

cat(paste("Spatial resolution: ", RESOLUTION, "\n", sep=""))
cat(paste("Coordinate rounding: ", ROUNDCOORDS, "\n", sep=""))

cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
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
load_pckg("h3")
load_pckg("sf")

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


## Function to suppress function output
## (will be used for silencing of garbage collector)
quiet <- function(x) { 
  sink("/dev/null")
  on.exit(sink()) 
  invisible(force(x)) 
}


## Parameters for debugging
# INPUT <- "/mnt/Dat2/GBIF/Parquet/2022-07-01/"
# PHYLUM <- NA
# CLASS <- NA
# ORDER <- NA
# FAMILY <- "Fabaceae"
# GENUS <- "Acacia"
# COUNTRY <- "AU"
# LATMIN <- NA
# LATMAX <- NA
# LONMIN <- NA
# LONMAX <- NA
# MINYEAR <- 2005
# EXTINCT <- NA
# EXCLUDEHUMAN <- TRUE
# TERRESTRIAL <- "/home/mik/GitRepos_Forks/biodiverse-scripts/pipeline_data/Land_Buffered_025_dgr.RData"
# CC_COUNTRY <- "/home/mik/GitRepos_Forks/biodiverse-scripts/pipeline_data/CC_CountryCentroids_buf_1000m.RData"
# CC_CAPITAL <- "/home/mik/GitRepos_Forks/biodiverse-scripts/pipeline_data/CC_Capitals_buf_10000m.RData"
# CC_INSTIT <- "/home/mik/GitRepos_Forks/biodiverse-scripts/pipeline_data/CC_Institutions_buf_100m.RData"
# CC_URBAN <- "/home/mik/GitRepos_Forks/biodiverse-scripts/pipeline_data/CC_Urban.RData"
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
  select(kingdom, phylum, class, order, family, genus, species, decimallongitude, decimallatitude) %>%
  count(kingdom, phylum, class, order, family, genus, species, decimallongitude, decimallatitude)

    # gbifid = unique record ID
    # institutioncode, collectioncode = a lot of missing values

## Collect the data
cat("..Collecting data\n")
datt <- dsf %>% collect()

## Convert to data.table
setDT(datt)

cat("Data table created\n")

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





################################################# Spatial binning (H3) and Estimation of the number of records

cat("Spatial binning\n")

## H3 system
datt[ , H3 := h3::geo_to_h3(datt[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

cat("..Number of H3-based gridcells: ", length(unique(datt$H3)), "\n")

cat("Counting number of unique elements per grid cell\n")

## Estimate count summary
datt_h3 <- datt[, .(
  NumRecords = sum(n, na.rm = TRUE),
  UniqKingdoms = length(unique(kingdom)),
  UniqPhyla = length(unique(phylum)),
  UniqClasses = length(unique(class)),
  UniqFamilies = length(unique(family)),
  UniqGenera = length(unique(genus)),
  UniqSpecies = length(unique(species))
  ), by = "H3"]


## Replace actual coordinates with gridcell centroid coordinates
datt_h3 <- cbind(datt_h3, h3::h3_to_geo(datt_h3$H3))
setnames(datt_h3, c("lat","lng"), c("decimallatitude","decimallongitude"))


################################################# Prepare outliers data

cat("Preparing outliers data\n")


#### Add grid cell IDs of removed samples as attributes to the resulting data

OUTLIERS <- list()

## Non-terrestrial outliers
if(!is.null(removed_nonterrestrial)){
  removed_nonterrestrial[ , H3 := h3::geo_to_h3(removed_nonterrestrial[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]
  removed_nonterrestrial[ , OutlierType := "NonTerrestrial" ]
  OUTLIERS[[ "NonTerrestrial" ]] <- removed_nonterrestrial
  rm(removed_nonterrestrial)
}

## Country and province centroids
if(!is.null(removed_CC_COUNTRY)){
  removed_CC_COUNTRY[ , H3 := h3::geo_to_h3(removed_CC_COUNTRY[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]
  removed_CC_COUNTRY[ , OutlierType := "CC_COUNTRY" ]
  OUTLIERS[[ "CC_COUNTRY" ]] <- removed_CC_COUNTRY
  rm(removed_CC_COUNTRY)
}

## Capital centroids
if(!is.null(removed_CC_CAPITAL)){
  removed_CC_CAPITAL[ , H3 := h3::geo_to_h3(removed_CC_CAPITAL[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]
  removed_CC_CAPITAL[ , OutlierType := "CC_CAPITAL" ]
  OUTLIERS[[ "CC_CAPITAL" ]] <- removed_CC_CAPITAL
  rm(removed_CC_CAPITAL)
}

## Institution centroids
if(!is.null(removed_CC_INSTIT)){
  removed_CC_INSTIT[ , H3 := h3::geo_to_h3(removed_CC_INSTIT[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]
  removed_CC_INSTIT[ , OutlierType := "CC_INSTIT" ]
  OUTLIERS[[ "CC_INSTIT" ]] <- removed_CC_INSTIT
  rm(removed_CC_INSTIT)
}

## Urban areas
if(!is.null(removed_CC_URBAN)){
  removed_CC_URBAN[ , H3 := h3::geo_to_h3(removed_CC_URBAN[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]
  removed_CC_URBAN[ , OutlierType := "CC_URBAN" ]
  OUTLIERS[[ "CC_URBAN" ]] <- removed_CC_URBAN
  rm(removed_CC_URBAN)
}

if(length(OUTLIERS) > 0){
  OUTLIERS <- rbindlist(OUTLIERS)
  setcolorder(x = OUTLIERS, neworder = c("OutlierType", "H3", "n"))
} else {
  OUTLIERS <- data.table(OutlierType = NA, H3 = NA, n = NA)
  OUTLIERS <- OUTLIERS[-1,]
}


## Summary
n_outliers <- sum(OUTLIERS$n, na.rm = TRUE)
n_nonoutliers <- sum(datt_h3$NumRecords, na.rm = TRUE)
perc_outliers <- n_outliers / (n_outliers + n_nonoutliers) * 100

cat("..Total number of nonoutlier records: ", n_nonoutliers, "\n")
cat("..Total number of GBIF-records marked as outliers: ", n_outliers, "\n")
cat("..Percetage of GBIF-records marked as outliers: ", round(perc_outliers, 2), "%\n")
cat("..Total number of H3 cells with outliers: ", length(unique(OUTLIERS$H3)), "\n")
cat("..Total number of species in outliers: ", length(unique(OUTLIERS$species)), "\n")


################################################# Export data

## Export
cat("Exporting filtered occurrence data\n")

## Main data with the number of records
cat("..Main data in R-format\n")
saveRDS(
  object = datt_h3,
  file = paste0(OUTPUT, "_H3.RData"),
  compress = "xz")

cat("..Main data in tsv-format\n")
fwrite(x = datt_h3,
  file = paste0(OUTPUT, "_H3.txt.gz"),
  quote = F, sep = "\t",
  compress = "gzip")


## Export outlier removal summary
cat("..Outliers\n")
fwrite(x = OUTLIERS,
  file = paste0(OUTPUT, "_Outliers.txt.gz"),
  quote = F, sep = "\t",
  compress = "gzip")


##################### Session info

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
