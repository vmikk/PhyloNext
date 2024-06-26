#!/usr/bin/env Rscript

## Visualization of Biodiverse results (will be deprecated in the future release)

cat("Basic visualization of Biodiverse results\n")
cat("Script name: 14_Visualization.R\n")

## Usage example:
# ./14_Visualization.R \
#   --observed "02.Biodiverse_input/occ.bds.csv" \
#   --zscores "02.Biodiverse_results/RND_rand--z_scores--SPATIAL_RESULTS.csv" \
#   --variables "RICHNESS_ALL,PD,PD_P" \
#   --resolution 4 \
#   --threads 1 \
#   --plotz "raw" \
#   --world "pipeline_data/WorldMap_NaturalEarth_Medium.RData" \
#   --output "03.Plots" \
#   --format "pdf" --width 18 --height 18 --units "in"


## TO DO:
# - fix dingbats symbols in PDF


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-r", "--observed"), action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - observed indices"),
  make_option(c("-z", "--zscores"),  action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - Z-scores"),
  make_option(c("-v", "--variables"), action="store", default="RICHNESS_ALL,PD,PD_P", type='character', help="Diversity variables to plot (comma-separated entries)"),
  make_option(c("--resolution"),  action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),
  make_option(c("-p", "--plotz"), action="store", default="raw", type='character', help="Plot raw estimates or Z-scores"),
  make_option(c("-w", "--world"), action="store", default=NA, type='character', help="File with contour map of the world"),
  make_option(c("--antimeridianfix"), action="store", default=TRUE, type='logical', help="Fix H3 polygons that cross the antimeridian"),
  make_option(c("-t", "--threads"), action="store", default=1L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-f", "--format"), action="store", default="pdf", type='character', help="Image format (pdf, png, svg, jpg)"),
  make_option(c("-n", "--width"), action="store", default=18, type='double', help="Image size, width"),
  make_option(c("-m", "--height"), action="store", default=18, type='double', help="Image size, height"),
  make_option(c("-u", "--units"), action="store", default="in", type='character', help="Image size units (in, cm, mm, px)"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$observed)){
  stop("Input file with observed PD indices is not specified.\n")
}
if(is.na(opt$zscores)){
  stop("Input file with Z-scores is not specified.\n")
}
if(is.na(opt$output)){
  stop("Output directory is not specified.\n")
}

## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
INPUTO     <- opt$observed
INPUTZ     <- opt$zscores
VARIABLES  <- opt$variables
RESOLUTION <- as.integer(opt$resolution)
PLOTZ      <- opt$plotz
WORLD      <- to_na( opt$world )
ANTIFIX    <- as.logical( opt$antimeridianfix )

CPUTHREADS <- as.numeric(opt$threads)
FORMAT <- opt$format
WIDTH <- as.numeric(opt$width)
HEIGTH <- as.numeric(opt$height)
UNITS <- opt$units
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input file (observed indices): ", INPUTO,     "\n", sep=""))
cat(paste("Input file (Z-scores): ",         INPUTZ,     "\n", sep=""))
cat(paste("Indices to plot: ",               VARIABLES,  "\n", sep=""))
cat(paste("Spatial resolution: ",            RESOLUTION, "\n", sep=""))
cat(paste("Variable type to plot (raw or Z-scores): ", PLOTZ, "\n", sep=""))
cat(paste("Adding world map: ",             WORLD,      "\n", sep=""))
cat(paste("Antimeridian fix: ",             ANTIFIX,    "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output image format: ",          FORMAT,     "\n", sep=""))
cat(paste("Output image width: ",           WIDTH,      "\n", sep=""))
cat(paste("Output image height: ",          HEIGTH,     "\n", sep=""))
cat(paste("Output image size units: ",      UNITS,      "\n", sep=""))
cat(paste("Output directory: ",             OUTPUT,     "\n", sep=""))

cat("\n")


## Create output directory
if(!OUTPUT %in% "."){
  dir.create(path = OUTPUT, showWarnings = F, recursive = TRUE)
}

############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("sf")
load_pckg("h3")
load_pckg("data.table")
load_pckg("plyr")
load_pckg("ggplot2")

cat("\n")

## ggplot2 theme
theme_set(theme_classic(base_size = 14))

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

## Load world data
if(! is.na(WORLD) ){
cat("Loading world map\n")
  world <- readRDS(WORLD)
}

## Load input data
cat("Loading Biodiverse results\n")

## Raw indices
cat("..Observed indices\n")
res_r <- fread(INPUTO)

## Z-scores
cat("..Z-scores\n")
res_z <- fread(INPUTZ)


## Rename coordinates
setnames(x = res_r,
  old = c("Axis_1", "Axis_0"),
  new = c("Latitude", "Longitude"))

setnames(x = res_z,
  old = c("Axis_1", "Axis_0"),
  new = c("Latitude", "Longitude"))


## Test if the first column cotains valid H3 IDs
if( h3_is_valid(res_r[[1,1]]) ){
  colnames(res_r)[1] <- "H3"
  colnames(res_z)[1] <- "H3"
} else {
  ## Get H3 IDs for grid cells
  cat("H3 index was not found in the data\n")
  cat("..Indexing geo-coordinates\n")
  res_r[ , H3 := h3::geo_to_h3(res_r[, .(Latitude, Longitude)], res = RESOLUTION) ]
  res_z[ , H3 := h3::geo_to_h3(res_z[, .(Latitude, Longitude)], res = RESOLUTION) ]
}

## If there are multiple variables selected - split them
if(any(grepl(pattern = ",", x = VARIABLES))){
  VARIABLES <- strsplit(x = VARIABLES, split = ",")[[1]]
}


## Check if the selected index is in the tables
colz <- unique(c(colnames(res_r), colnames(res_z)))
if(any(!VARIABLES %in% colz)){
  cat("Some of the selected indices are not present in tables with results!\n")
  cat("Please check the spelling of index names or eneble their estimation in Biodiverse.\n")
  missing <- VARIABLES[ ! VARIABLES %in% colz ]
  cat("Indices missing: ", paste(missing, collapse = ", "), "\n")

  ## Exclude missing indices
  VARIABLES <- VARIABLES[ ! VARIABLES %in% missing ]
}

if(length(VARIABLES) == 0){
  stop("None of the selected indices were found in the results! Nothing to plot.\n")
}


## Get spatial points
# H3_pts <- h3_to_geo_sf(res_z$H3)
# plot(H3_pts)

## Get spatial polygons
cat("Preparing gridcell polygons\n")
if(PLOTZ %in% c("raw", "Raw", "RAW")){

  H3_poly <- h3_to_geo_boundary_sf(res_r$H3)
  # plot(H3_poly)

  cat("..Adding raw diversity estimates to polygons\n")
  H3_poly <- cbind(H3_poly, res_r[, ..VARIABLES])
}
if(PLOTZ %in% c("z", "Z", "z-scores", "Z-scores")){

  H3_poly <- h3_to_geo_boundary_sf(res_z$H3)
  # plot(H3_poly)
  
  cat("..Adding Z-scores of diversity estimates to polygons\n")
  H3_poly <- cbind(H3_poly, res_z[, ..VARIABLES])
}


## Fix H3 polygons that cross the antimeridian by cutting them in two
if(ANTIFIX == TRUE){
  cat("..Fixing antimeridian issue\n")
  H3_poly <- st_wrap_dateline(H3_poly, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
}


## Find plot limits
boxx <- st_bbox(H3_poly)
xx <- pretty(c(boxx["xmin"], boxx["xmax"]))
yy <- pretty(c(boxx["ymin"], boxx["ymax"]))


## Plotting function
plot_function <- function(varname){

  if(is.na(WORLD)){
    PP <- ggplot(H3_poly) +
      geom_sf(aes_string(fill = varname), color = NA) +
      scale_fill_distiller(palette = "Spectral") + 
      ggtitle(varname) +
      xlim(xx[1], xx[length(xx)]) +
      ylim(yy[1], yy[length(yy)])
  }

  if(!is.na(WORLD)){
    PP <- ggplot(H3_poly) +
      geom_sf(data = world, fill = "grey95", color = "grey80") + 
      geom_sf(aes_string(fill = varname), color = NA) +
      scale_fill_distiller(palette = "Spectral") + 
      ggtitle(varname) +
      xlim(xx[1], xx[length(xx)]) +
      ylim(yy[1], yy[length(yy)])
  }

  return(PP)
}


## Function to export plot
plot_export <- function(varname){

  cat("..Plotting ", varname, "\n")

  ## Create a plot
  PP <- plot_function(varname)

  ## Make output file name
  OUTNAME <- file.path(OUTPUT, paste0(varname, ".", FORMAT))

  ## Export plot
  ggsave(filename = OUTNAME, plot = PP,
    width = WIDTH, height = HEIGTH, units = UNITS)

}


## Plot all variables
cat("Start plotting\n")

a_ply(.data = VARIABLES, .margins = 1, .fun = plot_export)

cat("Plotting finished\n")


## Export simple feature collection
cat("Exporting sf object\n")

saveRDS(object = H3_poly,
  file = file.path(OUTPUT, "H3_polygons.RData"),
  compress = "xz")


#####################

cat("All done!\n")

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
