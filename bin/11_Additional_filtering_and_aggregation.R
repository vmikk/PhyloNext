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
load_pckg("sf")
load_pckg("sp")
load_pckg("dbscan")

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


############################################## Main pipeline


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


## Remove non-terrestrial records 
if(!is.na(TERRESTRIAL)){
  cat("Removing non-terrestrial records\n")
  
  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points are on land
  land <- lengths(st_intersects(pts, TERRESTRIAL)) > 0

  non_terr <- sum(!land)
  cat("Number of non-terrestrial points = ", non_terr, "\n")

  ## Visualization
  # ggplot(data = pts) + 
  #   geom_sf(color = "red") + 
  #   geom_sf(data = TERRESTRIAL, fill = "grey")

  ## Remove outliers
  if(non_terr > 0){
    datt <- datt[ land ]
  }

  rm(pts)
}


## Density-based outlier removal
if(DBSCAN == TRUE){
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

    ## Remove outliers
    datt <- datt[ which(cl$cluster != 0 ) ]
  }

} # end of DBSCAN



######## Spatial aggregation
cat("Spatial aggregation using H3 system\n")

## H3 system
datt[ , H3 := h3::geo_to_h3(datt[, .(decimallatitude, decimallongitude)], res = RESOLUTION) ]

## Aggregate species by OTT IDs
datt_h3 <- unique(datt, by = c("specieskey", "H3"))

cat("Number of H3-based gridcells = ", length(unique(datt_h3$H3)), ")\n")

## Replace actual coordinates with gridcell centroid coordinates
uniq_h3 <- unique( datt_h3[, .(H3)] )
uniq_h3 <- cbind(uniq_h3, h3::h3_to_geo(uniq_h3$H3))

datt_h3[ , decimallongitude := NULL ]
datt_h3[ , decimallatitude := NULL ]
datt_h3 <- merge(x = datt_h3, y = uniq_h3, by = "H3", all.x = TRUE)
setnames(datt_h3, c("lat","lng"), c("decimallatitude","decimallongitude"))

