#!/usr/bin/env Rscript

## Additional filtering of GBIF occurrences and spatial binning

cat("Species-level filtering of occurrence data and spatial binning\n")
cat("Script name: 11_Additional_filtering_and_aggregation.R\n")

## TO DO:
# - handle a case with single observation
# - handle a case with no observations
# - in case of multiple species and DBSCAN - do it separately for each species
# - environmental outliers?

## By default, occurrences are binned using H3 resolution 4
## (grid cell area ~ 1770 km2; hexagon edge length ~ 22.6 km)

## The DBSCAN algorithm requires 2 parameters:
# eps: specifies how close points should be to each other to be considered a part of a cluster.
#      If the distance between two points is lower or equal to this value, these points are considered neighbors.
# minPoints: the minimum number of points to form a dense region.


## Usage:
# ./11_Additional_filtering_and_aggregation.R \
#    --input "/tmp/Fabaceae_in_AU.parquet" \
#    --specieskey "5349398" \
#    --resolution "4" \
#    --terrestrial "pipeline_data/Land_Buffered_025_dgr.RData" \
#    --wgsrpd "pipeline_data/WGSRPD.RData" \
#    --regions "L2_Australia" \
#    --rmcountrycentroids "pipeline_data/CC_CountryCentroids_buf_1000m.RData" \
#    --rmcountrycapitals "pipeline_data/CC_Capitals_buf_10000m.RData" \
#    --rminstitutions "pipeline_data/CC_Institutions_buf_100m.RData" \
#    --rmurban "pipeline_data/CC_Urban.RData" \
#    --output "Fabaceae_in_AU"


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
  make_option(c("-g", "--polygon"), action="store", default=NA, type='character', help="Custom area of interest (a file with polygons in GeoPackage format)"),
  make_option(c("-l", "--terrestrial"), action="store", default=NA, type='character', help="Remove non-terrestrial occurrences, provide land polygon in sf-format"),
  make_option(c("-w", "--wgsrpd"), action="store", default=NA, type='character', help="Path to the World Geographical Scheme for Recording Plant Distributions data (polygons in sf-format)"),
  make_option(c("-x", "--regions"), action="store", default=NA, type='character', help="Comma-separated list of WGSRPD regions"),

  ## CoordinateCleaner-like filters
  make_option(c("-c", "--rmcountrycentroids"), action="store", default=NA, type='character', help="Remove records within a radius around the geographic centroids of political countries and provinces"),
  make_option(c("-k", "--rmcountrycapitals"), action="store", default=NA, type='character', help="Remove records within a radius around country capitals"),
  make_option(c("-b", "--rminstitutions"), action="store", default=NA, type='character', help="Remove records in the vicinity of biodiversity institutions"),
  make_option(c("-u", "--rmurban"), action="store", default=NA, type='character', help="Remove records inside urban areas"),

  ## DBSCAN options
  make_option(c("-d", "--dbscan"), action="store", default=FALSE, type='logical', help="Remove spatial outliers with density-based clustering"),
  make_option(c("-e", "--epsilon"), action="store", default=700, type='double', help="DBSCAN parameter epsilon, km"),
  make_option(c("-p", "--minpts"), action="store", default=3, type='double', help="DBSCAN min number of points"),

  ## Spatial aggregation
  make_option(c("-r", "--resolution"), action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),

  make_option(c("-t", "--threads"), action="store", default=2L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option("--rcode", action="store", default="Shapefile_filters.R", type='character', help="File with R-code for spatial filtering"),

  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output directory")
  )
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$input)){
  cat("Input is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output directory is not specified.\n", file=stderr())
  stop()
}

## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
INPUT      <- opt$input
SPECIESKEY <- to_na( opt$specieskey )

POLYGON     <- to_na( opt$polygon )
TERRESTRIAL <- to_na( opt$terrestrial )

WGSRPD        <- to_na( opt$wgsrpd )
WGSRPDREGIONS <- to_na( opt$regions )

CC_COUNTRY <- to_na( opt$rmcountrycentroids )
CC_CAPITAL <- to_na( opt$rmcountrycapitals )
CC_INSTIT  <- to_na( opt$rminstitutions )
CC_URBAN   <- to_na( opt$rmurban )

DBSCAN <- as.logical( opt$dbscan )
DBSCAN_EPS <- as.numeric(opt$epsilon)
DBSCAN_PTS <- as.numeric(opt$minpts)

RESOLUTION <- as.integer(opt$resolution)
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input occurrences: ", INPUT,      "\n", sep=""))
cat(paste("GBIF specieskey: ",   SPECIESKEY, "\n", sep=""))

cat(paste("Custom polygons: ",  POLYGON,       "\n", sep=""))
cat(paste("Terrestrial data: ", TERRESTRIAL,   "\n", sep=""))
cat(paste("WGSRPD data: ",      WGSRPD,        "\n", sep=""))
cat(paste("WGSRPD regions: ",   WGSRPDREGIONS, "\n", sep=""))

cat(paste("Country and province centroids: ", CC_COUNTRY, "\n", sep=""))
cat(paste("Capitals: ",     CC_CAPITAL, "\n", sep=""))
cat(paste("Institutions: ", CC_INSTIT,  "\n", sep=""))
cat(paste("Uraban areas: ", CC_URBAN,   "\n", sep=""))

cat(paste("Perform DBSCAN: ", DBSCAN, "\n", sep=""))
if(DBSCAN == TRUE){
  cat(paste("DBSCAN parameter epsilon: ",    DBSCAN_EPS, "\n", sep=""))
  cat(paste("DBSCAN min number of points: ", DBSCAN_PTS, "\n", sep=""))
}

cat(paste("Spatial resolution: ", RESOLUTION, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output directory: ", OUTPUT, "\n", sep=""))

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
# load_pckg("sp")
load_pckg("dbscan")

cat("\n")

## Set CPU thread pool
cat("Number of available CPU threads: ",  cpu_count(), "\n")
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
set_cpu_count(CPUTHREADS)           # for libarrow
setDTthreads(threads = CPUTHREADS)  # for data.table


## Function to suppress function output
## (will be used for silencing of garbage collector)
quiet <- function(x) { 
  sink("/dev/null")
  on.exit(sink()) 
  invisible(force(x)) 
}

############################################## Main pipeline

## Open dataset
cat("Loading Parquet data\n")
ds <- arrow::open_dataset(INPUT)

## Select species and collect the data
cat("..Collecting data\n")
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
NRECORDS <- nrow(datt)  # will be added as attribute to the results
cat("Total number of records: ", NRECORDS, "\n")

if(is.na(SPECIESKEY)){ cat("No specieskey is specified\n") }
cat("There are ", length(unique(datt$specieskey)), "unique specieskeys in the data.\n")



## Spatial filtering (shapefile-based)
cat("\nSpatial filtering\n")
cat("..Loading a file with external code: ", opt$rcode, "\n")
source(opt$rcode)
cat("Shapefile-based filtering finished.\n")


## Density-based outlier removal
if(DBSCAN == TRUE & nrow(datt) > 0){
  cat("Density-based outlier removal\n")

  ## Reproject lat and long such that the Euclidean distance would approximate of geographic distance
  cat("..Reprojecting coordinates\n")

  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Great circle distances, in km
  # dd1 <- as.dist(st_distance(pts) / 1000)

  ## Create an azimuthal equidistant (AEQD) projection for the center-point of the study area
  dfc <- st_transform(
    x = pts,
    crs = "+proj=aeqd +x_0=0 +y_0=0 +lon_0=0 +lat_0=0 +units=km")

  ## Extract transformed coordinates
  crd <- st_coordinates(dfc)

  ## Euclidean distance
  # dd2 <- dist(crd)

  ## Compare distances
  # dd1[1:10]
  # dd2[1:10]
  # dist(datt[, .(decimallongitude, decimallatitude)])[1:10]  # incorrect dist!


  ## Specify the minimum size of a cluster proportional to the total number of occurrences
  if(DBSCAN_PTS < 1){

    ## How many records do we have?
    nsamps <- nrow(pts)

    ## If MinPts parameter is specified as a fraction, translate it into integer count
    DBSCAN_PTS <- ceiling(nsamps * DBSCAN_PTS / 100)

    cat("...MinPts parameter for DBSCAN was specified as a floating point number\n")
    cat("...and corresponds to ", DBSCAN_PTS, " spatial points.\n")

    if(DBSCAN_PTS < 2){
      cat("....Specified MinPts is < 2 and is too low. Setting MinPts = 2\n")
      DBSCAN_PTS <- 2
    }
  }    # end of proportional DBSCAN_PTS


  cat("..Running DBSCAN\n")
  ## Hierarchical DBSCAN (HDBSCAN)
  # cl <- hdbscan(datt[, .(decimallongitude, decimallatitude)], minPts = DBSCAN_PTS)
  
  ## Non-hierarchical (with fixed epsilon)
  cl <- dbscan(crd, minPts = DBSCAN_PTS, eps = DBSCAN_EPS)

  ## Visualization 
  # hullplot(datt[, .(decimallongitude, decimallatitude)], cl)

  num_outliers <- sum(cl$cluster %in% 0)
  cat("..Number of potential spatial outliers: ", num_outliers, "\n")

  if(num_outliers > 0){
    
    ## Add cluster ID to the data
    # datt$ClusterID <- cl$cluster

    ## Extract outliers
    removed_dbscan <- datt[ which(cl$cluster == 0 ) ]

    ## Remove outliers
    datt <- datt[ which(cl$cluster != 0 ) ]
  
  } else {
    removed_dbscan <- NULL    # no outliers found
  }

} else { # end of DBSCAN
    removed_dbscan <- NULL    # DBSCAN-based filtering was not performed
}

## Check the data
check_nodata(datt)


######## Spatial aggregation
cat("Spatial aggregation using H3 system\n")

## H3 system
datt[ , H3 := h3::geo_to_h3(datt[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

## Aggregate species by species key and gridcell
datt_h3 <- unique(datt, by = c("specieskey", "H3"))

cat("..Number of H3-based gridcells: ", length(unique(datt_h3$H3)), "\n")

## Replace actual coordinates with gridcell centroid coordinates
uniq_h3 <- unique( datt_h3[, .(H3)] )
uniq_h3 <- cbind(uniq_h3, h3::h3_to_geo(uniq_h3$H3))

datt_h3[ , decimallongitude := NULL ]
datt_h3[ , decimallatitude := NULL ]
datt_h3 <- merge(x = datt_h3, y = uniq_h3, by = "H3", all.x = TRUE)
setnames(datt_h3, c("lat","lng"), c("decimallatitude","decimallongitude"))


#### Add grid cell IDs of removed samples as attributes to the resulting data

## Non-terrestrial outliers
if(!is.null(removed_nonterrestrial)){
  
  ## H3 binning of non-terrestrial outliers
  removed_nonterrestrial[ , H3 := h3::geo_to_h3(removed_nonterrestrial[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_nonterrestrial_h3") <- unique(removed_nonterrestrial$H3)
  attr(datt_h3, which = "removed_nonterrestrial_n") <- nrow(removed_nonterrestrial)
} else {
  attr(datt_h3, which = "removed_nonterrestrial_h3") <- NA 
  attr(datt_h3, which = "removed_nonterrestrial_n") <- 0
}


## WGSRPD-outliers
if(!is.null(removed_WGSRPD)){
  
  ## H3 binning of non-removed_WGSRPD outliers
  removed_WGSRPD[ , H3 := h3::geo_to_h3(removed_WGSRPD[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_WGSRPD_h3") <- unique(removed_WGSRPD$H3)
  attr(datt_h3, which = "removed_WGSRPD_n") <- nrow(removed_WGSRPD)
} else {
  attr(datt_h3, which = "removed_WGSRPD_h3") <- NA 
  attr(datt_h3, which = "removed_WGSRPD_n") <- 0
}


## Country and province centroids
if(!is.null(removed_CC_COUNTRY)){
  
  ## H3 binning of removed occurrences
  removed_CC_COUNTRY[ , H3 := h3::geo_to_h3(removed_CC_COUNTRY[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_CC_COUNTRY_h3") <- unique(removed_CC_COUNTRY$H3)
  attr(datt_h3, which = "removed_CC_COUNTRY_n") <- nrow(removed_CC_COUNTRY)
} else {
  attr(datt_h3, which = "removed_CC_COUNTRY_h3") <- NA 
  attr(datt_h3, which = "removed_CC_COUNTRY_n") <- 0
}


## Capital centroids
if(!is.null(removed_CC_CAPITAL)){
  
  ## H3 binning of removed occurrences
  removed_CC_CAPITAL[ , H3 := h3::geo_to_h3(removed_CC_CAPITAL[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_CC_CAPITAL_h3") <- unique(removed_CC_CAPITAL$H3)
  attr(datt_h3, which = "removed_CC_CAPITAL_n") <- nrow(removed_CC_CAPITAL)
} else {
  attr(datt_h3, which = "removed_CC_CAPITAL_h3") <- NA 
  attr(datt_h3, which = "removed_CC_CAPITAL_n") <- 0
}


## Institution centroids
if(!is.null(removed_CC_INSTIT)){
  
  ## H3 binning of removed occurrences
  removed_CC_INSTIT[ , H3 := h3::geo_to_h3(removed_CC_INSTIT[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_CC_INSTIT_h3") <- unique(removed_CC_INSTIT$H3)
  attr(datt_h3, which = "removed_CC_INSTIT_n") <- nrow(removed_CC_INSTIT)
} else {
  attr(datt_h3, which = "removed_CC_INSTIT_h3") <- NA 
  attr(datt_h3, which = "removed_CC_INSTIT_n") <- 0
}


## Urban areas
if(!is.null(removed_CC_URBAN)){
  
  ## H3 binning of removed occurrences
  removed_CC_URBAN[ , H3 := h3::geo_to_h3(removed_CC_URBAN[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_CC_URBAN_h3") <- unique(removed_CC_URBAN$H3)
  attr(datt_h3, which = "removed_CC_URBAN_n") <- nrow(removed_CC_URBAN)
} else {
  attr(datt_h3, which = "removed_CC_URBAN_h3") <- NA 
  attr(datt_h3, which = "removed_CC_URBAN_n") <- 0
}




## DBSCAN-based outliers
if(!is.null(removed_dbscan)){

  ## H3 binning of DBSCAN outliers
  removed_dbscan[ , H3 := h3::geo_to_h3(removed_dbscan[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

  attr(datt_h3, which = "removed_dbscan_h3") <- unique(removed_dbscan$H3)
  attr(datt_h3, which = "removed_dbscan_n") <- nrow(removed_dbscan)
  attr(datt_h3, which = "dbscan_minpts") <- DBSCAN_PTS
  attr(datt_h3, which = "dbscan_epsilon") <- DBSCAN_EPS

} else {
  attr(datt_h3, which = "removed_dbscan_h3") <- NA 
  attr(datt_h3, which = "removed_dbscan_n") <- 0
  attr(datt_h3, which = "dbscan_minpts") <- NA
  attr(datt_h3, which = "dbscan_epsilon") <- NA
}

## Add total number of records (before filtering and spatial aggregation)
attr(datt_h3, which = "number_of_records") <- NRECORDS


## Ignore NA values
length_no_na <- function(x){
  if(length(x) == 1){ 
    if(is.na(x)){ 
      res <- 0
    } else {
      res <- 1
    }
  } else {
    res <- length(x)
  }
  return(res)
}

## Summarize filtering results
OUTLIERS_REMOVED <- rbind(
  data.frame(OutlierType = "NonTerrestrial",   N = attr(datt_h3, which = "removed_nonterrestrial_n"), NH3 = length_no_na(attr(datt_h3, which = "removed_nonterrestrial_h3"))),
  data.frame(OutlierType = "WGSRPD",           N = attr(datt_h3, which = "removed_WGSRPD_n"),         NH3 = length_no_na(attr(datt_h3, which = "removed_WGSRPD_h3"))),
  data.frame(OutlierType = "CountryCentroids", N = attr(datt_h3, which = "removed_CC_COUNTRY_n"),     NH3 = length_no_na(attr(datt_h3, which = "removed_CC_COUNTRY_h3"))),
  data.frame(OutlierType = "Capitals",         N = attr(datt_h3, which = "removed_CC_CAPITAL_n"),     NH3 = length_no_na(attr(datt_h3, which = "removed_CC_CAPITAL_h3"))),
  data.frame(OutlierType = "Institutions",     N = attr(datt_h3, which = "removed_CC_INSTIT_n"),      NH3 = length_no_na(attr(datt_h3, which = "removed_CC_INSTIT_h3"))),
  data.frame(OutlierType = "UrbanAreas",       N = attr(datt_h3, which = "removed_CC_URBAN_n"),       NH3 = length_no_na(attr(datt_h3, which = "removed_CC_URBAN_h3"))),
  data.frame(OutlierType = "DBSCAN",           N = attr(datt_h3, which = "removed_dbscan_n"),         NH3 = length_no_na(attr(datt_h3, which = "removed_dbscan_h3")))
  )
OUTLIERS_REMOVED <- cbind(SpeciesKey = SPECIESKEY, OUTLIERS_REMOVED)
setDT(OUTLIERS_REMOVED)
OUTLIERS_REMOVED[ is.na(N), NH3 := 0 ]



## Export
cat("Exporting filtered occurrence data\n")

## Create output directory if it doesn't exist
if(!OUTPUT %in% "."){
  dir.create(path = OUTPUT, showWarnings = F, recursive = TRUE)
}

## Prepare output name
if(is.na(SPECIESKEY)){
  OUTFILE <- file.path(OUTPUT, "NoSpKey.RData")
} else {
  OUTFILE <- file.path(OUTPUT, paste0(SPECIESKEY, ".RData"))
}

saveRDS(
  object = datt_h3,
  file = OUTFILE,
  compress = "xz")


## Export as TSV
# fwrite(x = datt_h3,
#   file = "export.txt.gz",
#   quote = F, sep = "\t", compress = "gzip")


## Export outlier removal summary
fwrite(x = OUTLIERS_REMOVED,
  file = gsub(pattern = ".RData$", replacement = "_OutlierCounts.txt", x = OUTFILE),
  quote = F, sep = "\t")



#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
