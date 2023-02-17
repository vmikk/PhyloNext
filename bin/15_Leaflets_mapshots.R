#!/usr/bin/env Rscript

## Create screenshots of a certain diversity/endemism metric based on Leaflet map

cat("Screenshots of Leaflet maps\n")
cat("Script name: 15_Leaflets_mapshots.R\n")

## Usage example:
# ./15_Leaflets_mapshots.R \
#   --input "Leaflet_object.RData" \
#   --metric "PD" \
#   --allmetrics "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE,PE_WE,SES_PE_WE,CANAPE,Redundancy"


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"),  action="store", default=NA, type='character', help="Input data with Leaflet object"),
  make_option(c("-m", "--metric"), action="store", default=NA, type='character', help="Diversity or endemicity metric name (e.g., `PD`)"),
  make_option(c("-a", "--allmetrics"), action="store", default=NA, type='character', help="Comma-separated list of all metric present in the Leaflet object"),

  make_option(c("-z", "--zoom"),  action="store", default=2,    type='integer', help="Zoom factor of the image"),
  make_option(c("--width"),       action="store", default=3840, type='integer', help="Image width"),
  make_option(c("--height"),      action="store", default=2160, type='integer', help="Image height")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$input)){
  stop("Input file with with Leaflet object is not specified.\n")
}
if(is.na(opt$metric)){
  stop("Diversity or endemicity metric name is not specified.\n")
}

## Assign variables
INPUT  <- opt$input
METRIC <- opt$metric
ALLMET <- opt$allmetrics

ZOOM   <- opt$zoom
WIDTH  <- opt$width
HEIGHT <- opt$height

## Log assigned variables
cat(paste("Input file: ",   INPUT,  "\n", sep=""))
cat(paste("Metric name: ",  METRIC, "\n", sep=""))
cat(paste("Metrics present in Leaflet: ",  ALLMET, "\n", sep=""))

cat(paste("Zoom factor: ",  ZOOM,   "\n", sep=""))
cat(paste("Image width: ",  WIDTH,  "\n", sep=""))
cat(paste("Image height: ", HEIGHT, "\n", sep=""))

cat("\n")


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}



load_pckg("leaflet")
load_pckg("mapview")
load_pckg("htmlwidgets")
load_pckg("chromote")
load_pckg("webshot2")
# load_pckg("webshot")

cat("\n")

## Which browser is used?
# chromote:::find_chrome()


############################################## Prepare data

## Parameters for debugging
# INPUT  <- "Leaflet_object.RData"
# METRIC <- "PD"
# ALLMET <- "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE"

## Load input data
cat(".. Loading leaflet object\n")
lf <- readRDS(INPUT)

## Prepare a vector of diversity index names for plotting
if(is.null(ALLMET)){
  ## Extract variable names from the object attributes
  VARIABLES <- attr(lf, "VARIABLES")
} else {
  ## Use user-supplied metirc list
  ## Split all metrics into vector
  VARIABLES <- strsplit(x = ALLMET, split = ",")[[1]]

}


## Check the number of grid cells in leaflet object,
## and specify the time to wait before taking screenshot.
## Sometimes a longer delay is needed for all gridcells to display properly.
cat(".. Finding the number of grid cells in the leaflet object\n")
NGRIDS <- length(lf$x$calls[[2]]$args[[1]])
cat("... Number of grids: ", NGRIDS, "\n")

if(NGRIDS < 500){
  DELAY <- 10
} else if (NGRIDS < 1000) {
  DELAY <- 5
} else if (NGRIDS < 2000) {
  DELAY <- 10
} else if (NGRIDS < 3000) {
  DELAY <- 15
} else if (NGRIDS < 4000) {
  DELAY <- 20
} else {
  DELAY <- 30
}
cat("... Setting approximate Chromote delay to: ", DELAY, " seconds \n")


## Function to export chloropleth as png
render_variable <- function(m, varname = "PD",
  zoom = 2, vwidth = 3840, vheight = 2160,
  delay = 1){
  # m <- lf

  ## Set coordinate bounds
  # m$x$limits
  # m %>% fitBounds(-44, -8, 113, 155)

  ## Variables are listed in the second last call
  # m$x$calls[[ length(m$x$calls) - 1 ]]

  ## Hide all variables, and select the one to show
  cat(".. Selecting variable to plot\n")
  m <- m %>% hideGroup(VARIABLES) %>% showGroup(group = varname)

  ## Remove legend
  # m <- m %>% clearControls()

  ## Remove controls
  # removeControl(m, layerId)

  ## Export plot
  # mapshot(
  #   m,
  #   url = NULL,
  #   file = "file.png",
  #   remove_controls = c("zoomControl", "layersControl", "homeButton", "scaleBar", "drawToolbar", "easyButton")
  # )
  # webshot v1 returns error

  cat(".. Saving temporary HTML widget\n")
  tmpfile <- paste0("temp_", varname, ".html")

  saveWidget(widget = m,
    file = tmpfile,
    selfcontained = FALSE)
  
  cat(".. Taking a screenshot\n")
  
  ## With webshot2 & Chrome
  webshot2::webshot(url = tmpfile,
    file = paste0(varname, ".png"),
    cliprect = "viewport",
    delay = delay,
    zoom = zoom,
    vwidth = vwidth,
    vheight = vheight)

  ## With webshot1 & PhantomJS
  # webshot::webshot(url = "temp.html",
  #   file = paste0(varname, ".png"),
  #   cliprect = "viewport",
  #   delay = 1,
  #   zoom = zoom,
  #   vwidth = vwidth,
  #   vheight = vheight)

  ## Remove tmp files
  cat(".. Removing temporary files\n")
  system(paste0("rm ", tmpfile))
  system(paste0("rm -r ", paste0("temp_", varname), "_files"))

}


if(NGRIDS > 5000){

  cat("\nWARNING: The number of cells is too high! Mapshot/Chromote will likely to fail\n")
  cat(".. Skipping mapshot creation\n")

} else {

  ## Run the function
  render_variable(
    m = lf,
    varname = METRIC,
    zoom = ZOOM, vwidth = WIDTH, vheight = HEIGHT,
    delay = DELAY)

  cat("Plotting finished\n")

}



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





######## manual chromote

# vwidth = 3840, vheight = 2160, zoom = 2


# url <- webshot2:::file_url(url)
# chromote <- chromote:::default_chromote_object()

# cliprect <- c(0, 0, vwidth, vheight)
# expand   <- list(expand)




# b <- chromote$new_session()
# b$Browser$getVersion()
# # b$view()



# s <- NULL
# p <- chromote$new_session(wait_ = FALSE,
#       width = vwidth,
#       height = vheight
#     )$
#     then(function(session) {
#       s <<- session

#       if (!is.null(useragent)) {
#         s$Network$setUserAgentOverride(userAgent = useragent)
#       }
#       res <- s$Page$loadEventFired(wait_ = FALSE)
#       s$Page$navigate(url, wait_ = FALSE)
#       res
#     })$
#     then(function(value) {
#       if (delay > 0) {
#         promise(function(resolve, reject) {
#           later(
#             function() {
#               resolve(value)
#             },
#             delay
#           )
#         })
#       } else {
#         value
#       }
#     })$
#     then(function(value) {
#         s$screenshot(
#           filename = file, selector = "html", cliprect = cliprect,
#           expand = expand, scale = zoom,
#           show = FALSE, wait_ = FALSE
#         )
#     })$
#     then(function(value) {
#       message(url, " screenshot completed")
#       normalizePath(value)
#     })$
#     finally(function() {
#       s$close()
#     })
