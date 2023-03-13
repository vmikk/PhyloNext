#!/usr/bin/env Rscript

## Prepare a shapefile for spatially-constrained randomisations

cat("Preparing a shapefile with spatial polygons (to constrain Biodiverse randomizations)\n")
cat("Script name: 12_Prepare_SpatialConstraints.R\n")


## Usage example:
# ./12_Prepare_SpatialConstraints.R \
#   --input "Trimmed_occurrences.csv" \
#   --randconstrain "ZoogeographicRegions.gpkg" \
#   --threads 1 \
#   --output "shapefile"

## Output = shapefile

## NB. in GeoPackage file, the first column should contain polygon name.

## TO DO:
# - add support for disjunkt polygons? (same name, but multiple polygons) ?


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Species occurrence dataset"),
  make_option(c("-p", "--randconstrain"), action="store", default=NA, type='character', help="Polygons to perform spatially constrained randomization (GeoPackage format)"),
  make_option(c("-t", "--threads"), action="store", default=2L, type='integer', help="Number of CPU threads for arrow, default 2"),
  make_option(c("-o", "--output"),  action="store", default="shapefile", type='character', help="Output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$input)){
  stop("Input file is not specified.\n")
}
if(is.na(opt$randconstrain)){
  stop("GeoPackage file is not specified.\n")
}


## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
INPUT      <- opt$input
POLYG      <- opt$randconstrain
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT     <- opt$output

## Log assigned variables
cat(paste("Species occurrence dataset: ",   INPUT, "\n", sep=""))
cat(paste("Polygons: ",                     POLYG, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output file prefix: ",           OUTPUT, "\n", sep=""))

cat("\n")




############################################## Data for debugging

# INPUT      <- "Trimmed_occurrences.csv"
# POLYG      <- "ZoogeographicRegions.gpkg"
# CPUTHREADS <- 4
# OUTPUT     <- "shapefile"

############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("sf")
load_pckg("plyr")

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

## Load data
cat("Loading species occurrences\n")
datt <- fread(INPUT)
cat(".. Number of records detected: ", nrow(datt), "\n")

cat("Loading spatial polygons\n")
poly <- st_read(POLYG)
cat(".. Number of polygons detected: ", nrow(poly), "\n")

## Drop non-relevant information (keep only the first column)
poly <- poly[,1]
names(poly)[1] <- "Polygon"

## Check if polygon names are unique
polynames <- poly[,1][[1]]    # polygon names must by in the first column!
if(length(polynames) != length(unique(polynames))){
  cat("WARNING: Polygon names (based on the first columns in the GeoPackage file) must be unique!\n")
  cat("WARNING: Disjunkt polygons are not supported yet!\n")
}


## Inspect polygons
cat("Inspecting polygons\n")

cat(".. Checking CRS\n")
cat("... Geodetic CRS of the provided file is:\n")
st_crs(poly)

if(! identical(st_crs("EPSG:4326"), st_crs(poly))){
  cat("... Converting CRS to EPSG:4326\n")
  poly <- st_transform(poly, 4326)
}


cat(".. Preparing simple feature collection\n")
datt[, TmpID := 1:.N ]   # add record key
xy <- st_as_sf(
  x = datt[ , .(TmpID, decimallongitude, decimallatitude) ],
  coords = c("decimallongitude", "decimallatitude"),
  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Find which polygon a point belong to
cat(".. Spatial joining of points with polygons\n")
xy <- st_join(x = xy, y = poly)

## Check if some points are outside the provided polygons
if(any(is.na(xy$Polygon))){
  cat("WARNING: Some points in the occurrence dataset are outside of the provided polygons!\n")
  cat("WARNING: N = ", sum(is.na(xy$Polygon)), "\n")

  cat("WARNING: These points will be removed\n")
  xy <- xy[ ! is.na(xy$Polygon), ]
  if(nrow(xy) == 0){
    stop("No species occurrence belongs to the spatial polygons!\n")
  }
  cat("... Number of records inside polygons: ", nrow(xy), "\n")
}


cat(".. Subsetting polygons\n")

## Which polygons should we keep?
upoly <- unique(xy$Polygon)

cat("... Species occurrence records belongs to ", length(upoly), " polygons\n")
print(table(xy$Polygon))

poly <- poly[ poly$Polygon %in% upoly, ]


cat("Exporting shapefile\n")

st_write(poly, paste0(OUTPUT, ".shp"))


cat("\nAll done!\n")

#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
