#!/usr/bin/env Rscript

## Count number of records per gridcell in the GBIF occurrence data

cat("Counting the number of records per grid cell\n")
cat("Script name: 10_Record_counts.R\n")

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
  make_option("--maxyear", action="store", default=NA,   type='integer', help="Maximum year of occurrence"),
  make_option("--noextinct", action="store", default=NA, type='character', help="Remove extinct species (provide a file with extinct specieskeys)"),
  make_option("--excludehuman", action="store", default=TRUE, type='logical', help="Exclude human records (genus Homo)"),
  make_option("--basisofrecordinclude",  action="store", default=NA, type='character', help="Basis of record to include from the data"),
  make_option("--basisofrecordexclude", action="store", default="FOSSIL_SPECIMEN,LIVING_SPECIMEN", type='character', help="Basis of record to exclude from the data"),

  ## Coordinate precision and uncertainty filters
  make_option("--coordprecision", action="store", default=0.1, type='double', help="Coordinate precision threshold (max allowed value)"),
  make_option("--coorduncertainty", action="store", default=10000, type='double', help="Maximum allowed coordinate uncertainty, meters"),
  make_option("--coorduncertaintyexclude", action="store", default="301,3036,999,9999", type='character', help="Black-listed values of coordinate uncertainty"),

  ## Spatial filters (shapefile-based)
  make_option(c("-g", "--polygon"), action="store", default=NA, type='character', help="Custom area of interest (a file with polygons in GeoPackage format)"),
  make_option(c("-l", "--terrestrial"), action="store", default=NA, type='character', help="Remove non-terrestrial occurrences, provide land polygon in sf-format"),
  make_option(c("-c", "--rmcountrycentroids"), action="store", default=NA, type='character', help="Remove records within a radius around the geographic centroids of political countries and provinces"),
  make_option(c("-k", "--rmcountrycapitals"), action="store", default=NA, type='character', help="Remove records within a radius around country capitals"),
  make_option(c("-b", "--rminstitutions"), action="store", default=NA, type='character', help="Remove records in the vicinity of biodiversity institutions"),
  make_option(c("-u", "--rmurban"), action="store", default=NA, type='character', help="Remove records inside urban areas"),
  make_option(c("-w", "--wgsrpd"), action="store", default=NA, type='character', help="Path to the World Geographical Scheme for Recording Plant Distributions data (polygons in sf-format)"),
  make_option(c("-x", "--regions"), action="store", default=NA, type='character', help="Comma-separated list of WGSRPD regions"),

  ## Spatial aggregation
  make_option(c("-r", "--resolution"), action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),
  make_option(c("--roundcoords"), action="store", default=2L, type='integer', help="Round spatial coordinates to the N decimal places, to reduce the dataset size (default, 2). To disable, set to a negative value."),

  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option("--rcode", action="store", default="Shapefile_filters.R", type='character', help="File with R-code for spatial filtering"),

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

COORDPREC      <- as.numeric( to_na(opt$coordprecision) )
COORDUNCRTMAX  <- as.numeric( to_na(opt$coorduncertainty) )
COORDUNCRTEXCL <- to_na(opt$coorduncertaintyexclude)

COUNTRY <- to_na( opt$country )
LATMIN  <- as.numeric( to_na(opt$latmin) )
LATMAX  <- as.numeric( to_na(opt$latmax) )
LONMIN  <- as.numeric( to_na(opt$lonmin) )
LONMAX  <- as.numeric( to_na(opt$lonmax) )

MINYEAR <- as.numeric(to_na( opt$minyear) )
MAXYEAR <- as.numeric(to_na( opt$maxyear) )

EXTINCT      <- to_na( opt$noextinct)
EXCLUDEHUMAN <- as.logical( opt$excludehuman )

BASISINCL <- to_na( opt$basisofrecordinclude )
BASISEXCL <- to_na( opt$basisofrecordexclude )

POLYGON       <- to_na( opt$polygon )
WGSRPD        <- to_na( opt$wgsrpd )
WGSRPDREGIONS <- to_na( opt$regions )
TERRESTRIAL   <- to_na( opt$terrestrial )
CC_COUNTRY    <- to_na( opt$rmcountrycentroids )
CC_CAPITAL    <- to_na( opt$rmcountrycapitals )
CC_INSTIT     <- to_na( opt$rminstitutions )
CC_URBAN      <- to_na( opt$rmurban )

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
cat(paste("Maximum year of occurrence: ", MAXYEAR, "\n", sep=""))
cat(paste("List of extict species: ",     EXTINCT, "\n", sep=""))
cat(paste("Exclusion of human records: ", EXCLUDEHUMAN, "\n", sep=""))
cat(paste("Round coordinates: ",          ROUNDCOORDS, "\n", sep=""))

cat(paste("Custom polygons: ",  POLYGON,       "\n", sep=""))
cat(paste("WGSRPD data: ",      WGSRPD,        "\n", sep=""))
cat(paste("WGSRPD regions: ",   WGSRPDREGIONS, "\n", sep=""))
cat(paste("Terrestrial data: ", TERRESTRIAL,   "\n", sep=""))
cat(paste("Country and province centroids: ", CC_COUNTRY, "\n", sep=""))
cat(paste("Capitals: ",     CC_CAPITAL, "\n", sep=""))
cat(paste("Institutions: ", CC_INSTIT, "\n", sep=""))
cat(paste("Uraban areas: ", CC_URBAN, "\n", sep=""))

cat(paste("Spatial resolution: ",  RESOLUTION, "\n", sep=""))
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
# INPUT <- "/mnt/Dat2/GBIF/Parquet/2022-08-01/"
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
# MAXYEAR <- NA
# EXTINCT <- NA
# EXCLUDEHUMAN <- TRUE
# BASISINCL <- NA
# BASISEXCL <- "FOSSIL_SPECIMEN,LIVING_SPECIMEN"
# POLYGON <- NA
# WGSRPD <- NA
# WGSRPDREGIONS <- NA
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
  cat("..Filtering by collection date (min year)\n")
  dsf <- dsf %>% filter(year >= MINYEAR)
}
if(!is.na(MAXYEAR)){
  cat("..Filtering by collection date (max year)\n")
  dsf <- dsf %>% filter(year <= MAXYEAR)
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

cat("\nSpatial filtering\n")
cat("..Loading a file with external code: ", opt$rcode, "\n")
source(opt$rcode)
cat("Shapefile-based filtering finished.\n")


################################################# Spatial binning (H3) and Estimation of the number of records

cat("\nSpatial binning\n")

## H3 system
datt[ , H3 := h3::geo_to_h3(datt[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

cat("..Number of H3-based gridcells: ", length(unique(datt$H3)), "\n")


## Estimate number of records per species
cat("Counting number of records per grid cell per species\n")

datt_h3_sp <- datt[, .(
  NumRecords = sum(n, na.rm = TRUE)
  ), by = c("H3", "species")]


## Estimate global count summary
cat("Counting number of unique elements per grid cell\n")

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
cat("Exporting total record counts\n")

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



## Per-species counts
cat("Exporting per-species record counts\n")

cat("..Main data in R-format\n")
saveRDS(
  object = datt_h3_sp,
  file = paste0(OUTPUT, "_H3_PerSpecies.RData"),
  compress = "xz")

cat("..Main data in tsv-format\n")
fwrite(x = datt_h3_sp,
  file = paste0(OUTPUT, "_H3_PerSpecies.txt.gz"),
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
