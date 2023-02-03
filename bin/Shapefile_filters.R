## This file contains code for spatial filtering (shapefile-based)
## which is shared among different workflow steps
## Sourced by the following scripts:
##  `10_Record_counts.R`
##  `11_Additional_filtering_and_aggregation.R`
##  `16_Derived_dataset.R`

## The following objects must be present in the environment:
# - datt
# - POLYGON
# - WGSRPD
# - WGSRPDREGIONS
# - TERRESTRIAL
# - CC_COUNTRY
# - CC_CAPITAL
# - CC_INSTIT
# - CC_URBAN


## Function to handle case with no observations
## (in the case if everithing was removed during the filtering)
check_nodata <- function(x){
  if(! nrow(x) > 0){ 
    cat("\nWARNING: no speciecies occurrences after the filtering\n")
  }
}


## Extract records belonging to a custom polygon
if(!is.na(POLYGON) & nrow(datt) > 0){
  cat("Extracting records belonging to a user-supplied polygons\n")

  ## Load polygons
  cat("..Loading GeoPackge file\n")
  POLYGON <- read_sf(POLYGON)

  cat("..Number of polygons detected: ", nrow(POLYGON), "\n")
  if(nrow(POLYGON) > 0){

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belongs to the polygon
  cat("..Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, POLYGON)
  poly <- lengths(poly_intersect) > 0
  non_poly <- sum(!poly)

  cat("..Number of points inside polygons: ",  sum(poly), "\n")
  cat("..Number of points outside polygons: ", non_poly,  "\n")

  ## Remove outliers
  if(sum(poly) > 0){ 
    datt <- datt[ poly ]
  } else {
  	cat("..ALL POINTS ARE OUTSIDE THE SUPPLIED POLYGONS!\n")
  	datt <- datt[ -c(1:nrow(datt)) ]
  }

  ## Check the data
  check_nodata(datt)
  
  ## Clean up
  rm(pts, poly_intersect, poly, non_poly, POLYGON)
  quiet( gc() )
  }
}


## Subset to WGSRPD regions
if(!is.na(WGSRPD) & !is.na(WGSRPDREGIONS) & nrow(datt) > 0){
  cat("Subsetting by World Geographical Regions\n")

  ## Load WGSRPD polygons
  cat("..Loading WGSRPD polygons\n")
  WGSRPD <- readRDS(WGSRPD)

  ## Split the selected regions (if there are more than one)
  WGSRPDREGIONS <- strsplit(x = WGSRPDREGIONS, split = ",")[[1]]

  ## Check if selected regions are in the shapefile
  in_polygons <- WGSRPDREGIONS %in% WGSRPD$LevelName
  if(any(!in_polygons)){
    cat("\n...! Check: ", WGSRPDREGIONS[!in_polygons], "!...\n")
    stop("ERROR: Some of the selected regions are not in the shapefile!\n")
  }

  ## Extract polygons of the selected regions
  POLY <- WGSRPD[ which(WGSRPD$LevelName %in% WGSRPDREGIONS),  ]

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("..Intersecting world polygons and data points\n")
  poly_intersect <- st_intersects(pts, POLY)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("..Number of points inside the selected polygons: ",  sum(poly), "\n")
  cat("..Number of points outside the selected polygons: ", non_poly,  "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_WGSRPD <- datt[ ! poly ]   # outliers
    datt <- datt[ poly ]
  } else {
    removed_WGSRPD <- NULL   # no non-WGSRPD samples found
  }

  check_nodata(datt)
  rm(pts, WGSRPD)
  quiet( gc() )
} else {
  removed_WGSRPD <- NULL     # no WGSRPD-filtering was performed
}



## Remove non-terrestrial records 
if(!is.na(TERRESTRIAL) & nrow(datt) > 0){
  cat("Removing non-terrestrial records\n")
  
  ## Load land mask
  cat("..Loading land mask\n")
  TERRESTRIAL <- readRDS(TERRESTRIAL)

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points are on land
  cat("..Intersecting polygons and data points\n")
  land_intersect <- st_intersects(pts, TERRESTRIAL)
  land <- lengths(land_intersect) > 0

  non_terr <- sum(!land)
  cat("..Number of non-terrestrial points: ", non_terr, "\n")

  ## Visualization
  # ggplot(data = pts) + 
  #   geom_sf(color = "red") + 
  #   geom_sf(data = TERRESTRIAL, fill = "grey")

  ## Remove outliers
  if(non_terr > 0){
    removed_nonterrestrial <- datt[ ! land ]   # outliers
    datt <- datt[ land ]
  } else {
    removed_nonterrestrial <- NULL   # no non-terrestrial samples found
  }

  check_nodata(datt)
  rm(pts, TERRESTRIAL)
  quiet( gc() )
} else {
  removed_nonterrestrial <- NULL     # no terrestrial filtering was performed
}



## Remove country centroids and province centroids (similar to CoordinateCleaner::cc_cen)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_cen.R
## Default buffer = 1 km
if(!is.na(CC_COUNTRY) & nrow(datt) > 0){
  cat("Removing occurrences inside country and province centroids\n")

  ## Load polygons
  cat("..Loading polygons\n")
  CC_COUNTRY <- readRDS(CC_COUNTRY)
  cat("..Number of polygons provided: ", nrow(CC_COUNTRY), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("..Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_COUNTRY)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("..Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("..Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_COUNTRY <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_COUNTRY <- NULL   # no found
  }

  check_nodata(datt)
  rm(pts, CC_COUNTRY)
  quiet( gc() )
} else {
  removed_CC_COUNTRY <- NULL     # no filtering was performed
}


## Filter country capitals (similar to CoordinateCleaner::cc_cap)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_cap.R
## Default buffer = 10 km
if(!is.na(CC_CAPITAL) & nrow(datt) > 0){
  cat("Removing occurrences within capitals\n")

  ## Load polygons
  cat("..Loading polygons\n")
  CC_CAPITAL <- readRDS(CC_CAPITAL)
  cat("..Number of polygons provided: ", nrow(CC_CAPITAL), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("..Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_CAPITAL)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("..Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("..Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_CAPITAL <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_CAPITAL <- NULL   # no found
  }

  check_nodata(datt)
  rm(pts, CC_CAPITAL)
  quiet( gc() )
} else {
  removed_CC_CAPITAL <- NULL     # no filtering was performed
}


## Remove records in the vicinity of biodiversity institutions (similar to CoordinateCleaner::cc_inst)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_inst.R
## Default buffer = 100 m
if(!is.na(CC_INSTIT) & nrow(datt) > 0){
  cat("Removing occurrences in the vicinity of biodiversity institutions\n")

  ## Load polygons
  cat("..Loading polygons\n")
  CC_INSTIT <- readRDS(CC_INSTIT)
  cat("..Number of polygons provided: ", nrow(CC_INSTIT), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("..Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_INSTIT)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("..Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("..Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_INSTIT <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_INSTIT <- NULL   # no found
  }

  check_nodata(datt)
  rm(pts, CC_INSTIT)
  quiet( gc() )
} else {
  removed_CC_INSTIT <- NULL     # no filtering was performed
}



## Remove records inside urban areas (similar to CoordinateCleaner::cc_urb)
## https://github.com/ropensci/CoordinateCleaner/blob/master/R/cc_urb.R
if(!is.na(CC_URBAN) & nrow(datt) > 0){
  cat("Removing occurrences inside urban areas\n")

  ## Load polygons
  cat("..Loading polygons\n")
  CC_URBAN <- readRDS(CC_URBAN)
  cat("..Number of polygons provided: ", nrow(CC_URBAN), "\n")

  ## Convert coordinates to sf class
  pts <- st_as_sf(
    x = datt[, .(decimallongitude, decimallatitude)],
    coords = c("decimallongitude", "decimallatitude"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## Check which points belong to the polygons
  cat("..Intersecting polygons and data points\n")
  poly_intersect <- st_intersects(pts, CC_URBAN)
  poly <- lengths(poly_intersect) > 0

  non_poly <- sum(!poly)
  cat("..Number of points inside the selected polygons: ", sum(poly), "\n")
  cat("..Number of points outside the selected polygons: ", non_poly, "\n")

  ## Remove outliers
  if(non_poly > 0){
    removed_CC_URBAN <- datt[ poly ]   # outliers
    datt <- datt[ ! poly ]
  } else {
    removed_CC_URBAN <- NULL   # no found
  }

  check_nodata(datt)
  rm(pts, CC_URBAN)
  quiet( gc() )
} else {
  removed_CC_URBAN <- NULL     # no filtering was performed
}


