#!/usr/bin/env Rscript

## Visualization of Biodiverse results

## Usage example:
# ./14_Visualization.R \
#   --observed "02.Biodiverse_input/occ.bds.csv" \
#   --zscores "02.Biodiverse_results/RND_rand--z_scores--SPATIAL_RESULTS.csv" \
#   --threads 1 \
#   --variables "PHYLO_RPD1" \
#   --world "pipeline_data/WorldMap_NaturalEarth_Medium.RData" \
#   --output "03.Plots"

############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-r", "--observed"), action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - observed indices"),
  make_option(c("-z", "--zscores"),  action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - Z-scores"),
  make_option(c("-v", "--variables"), action="store", default="PHYLO_RPD1", type='character', help="Diversity variables to plot"),
  make_option(c("-w", "--world"), action="store", default=NA, type='character', help="File with contour map of the world"),
  make_option(c("-t", "--threads"), action="store", default=1L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-f", "--format"), action="store", default="pdf", type='character', help="Image format (pdf, png, svg, jpg)"),
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
  if(x %in% "NA"){ x <- NA }
  return(x)
}

## Assign variables
INPUTO <- opt$observed
INPUTZ <- opt$zscores
VARIABLES <- opt$variables
WORLD <- to_na( opt$world )

CPUTHREADS <- as.numeric(opt$threads)
FORMAT <- opt$format
OUTPUT <- opt$output


## Log assigned variables
cat(paste("Input file (observed indices): ", INPUTO, "\n", sep=""))
cat(paste("Input file (Z-scores): ", INPUTZ, "\n", sep=""))
cat(paste("Variables to plot: ", VARIABLES, "\n", sep=""))
cat(paste("Adding world map: ", WORLD, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output image format: ", FORMAT, "\n", sep=""))
cat(paste("Output directory: ", OUTPUT, "\n", sep=""))

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

## The first columns should be a gridcell
colnames(res_r)[1] <- "H3"
colnames(res_z)[1] <- "H3"


## Check if the selected index is in the tables
#                                             ///// update for multiple inds
if(!VARIABLES %in% colnames(res_r)){
  stop("Selected index is not present in the table with observed indices!\n")
}
if(!VARIABLES %in% colnames(res_z)){
  stop("Selected index is not present in the table with Z-scores!\n")
}


## Get spatial points
# H3_pts <- h3_to_geo_sf(res_z$H3)
# plot(H3_pts)

## Get spatial polygons
cat("Preparing gridcell polygons\n")
H3_poly <- h3_to_geo_boundary_sf(res_z$H3)
# plot(H3_poly)


## Add diversity indices to the data
cat("Adding diversity estimates to polygons\n")

H3_poly <- cbind(H3_poly, res_z[, ..VARIABLES])


## Find plot limits
boxx <- st_bbox(H3_poly)
xx <- pretty(c(boxx["xmin"], boxx["xmax"]))
yy <- pretty(c(boxx["ymin"], boxx["ymax"]))


if(is.na(WORLD)){
  PP <- ggplot(H3_poly) +
    geom_sf(aes_string(fill = VARIABLES), color = NA) +
    scale_fill_distiller(palette = "Spectral") + 
    ggtitle(VARIABLES) +
    xlim(xx[1], xx[length(xx)]) +
    ylim(yy[1], yy[length(yy)])
}

if(!is.na(WORLD)){
  PP <- ggplot(H3_poly) +
    geom_sf(data = world, fill = "grey95", color = "grey80") + 
    geom_sf(aes_string(fill = VARIABLES), color = NA) +
    scale_fill_distiller(palette = "Spectral") + 
    ggtitle(VARIABLES) +
    xlim(xx[1], xx[length(xx)]) +
    ylim(yy[1], yy[length(yy)])
}

OUTNAME <- file.path(OUTPUT, paste0(VARIABLES, ".", FORMAT))

ggsave(filename = OUTNAME,
  plot = PP, width = 18, height = 18)



#####################

cat("Plotting finished\n")

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
