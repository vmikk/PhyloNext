#!/usr/bin/env Rscript

## Visualization of Biodiverse results using intaractive maps (Leaflet-based)

## Usage example:
# ./15_Leaflets.R \
#   --observed "02.Biodiverse_input/occ.bds.csv" \
#   --zscores "02.Biodiverse_results/RND_rand--z_scores--SPATIAL_RESULTS.csv" \
#   --variables "RICHNESS_ALL,PD,zPD,PD_P,zPD_P" \
#   --palette "quantile" \
#   --color "RdYlBu" \
#   --bins 5 \
#   --output "03.Plots/Choropleth.html"

## for Z-score-based variables, add `z` prefix (e.g., zPD)


## TO DO:
# - handle NA values (e.g., if Richness == 1)


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-r", "--observed"), action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - observed indices"),
  make_option(c("-z", "--zscores"),  action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - Z-scores"),
  make_option(c("-v", "--variables"), action="store", default="RICHNESS_ALL,PD,zPD,PD_P,zPD_P", type='character', help="Diversity variables to plot (comma-separated entries)"),
  make_option(c("-p", "--palette"), action="store", default="quantile", type='character', help="Color palette type"),
  make_option(c("-c", "--color"), action="store", default="RdYlBu", type='character', help="Color gradient scheme"),
  make_option(c("-b", "--bins"), action="store", default=5L, type='integer', help="Number of color bins"),
  # make_option(c("-t", "--threads"), action="store", default=1L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default="Choropleth.html", type='character', help="Output file name")
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
  stop("Output file is not specified.\n")
}

## Assign variables
INPUTO <- opt$observed
INPUTZ <- opt$zscores
VARIABLES <- opt$variables
PALETTE <- opt$palette
COLOR <- opt$color
BINS <- as.numeric( opt$bins )
OUTPUT <- opt$output

## Log assigned variables
cat(paste("Input file (observed indices): ", INPUTO, "\n", sep=""))
cat(paste("Input file (Z-scores): ", INPUTZ, "\n", sep=""))
cat(paste("Indices to plot: ", VARIABLES, "\n", sep=""))
cat(paste("Color palette type: ", PALETTE, "\n", sep=""))
cat(paste("Color gradient scheme: ", COLOR, "\n", sep=""))
cat(paste("Number of color bins: ", BINS, "\n", sep=""))
cat(paste("Output file: ", OUTPUT, "\n", sep=""))

# CPUTHREADS <- as.numeric(opt$threads)
# cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))


cat("\n")


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
load_pckg("leaflet")
load_pckg("mapview")

cat("\n")


## Set CPU thread pool
# cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
# setDTthreads(threads = CPUTHREADS)  # for data.table

# ## Start local cluster
# if(CPUTHREADS > 1){
#   load_pckg("doFuture")
#   registerDoFuture()
#   plan(multicore, workers = CPUTHREADS)
#   options(future.globals.maxSize = 1e10)
# 
#   parall <- TRUE
# } else {
#   parall <- FALSE
# }


############################################## Prepare data

# INPUTO <- "/mnt/Dat2/GBIF/Fabaceae_Nextflow_220324/02.Biodiverse_results/RND_SPATIAL_RESULTS.csv"
# INPUTZ <- "/mnt/Dat2/GBIF/Fabaceae_Nextflow_220324/02.Biodiverse_results/RND_rand--z_scores--SPATIAL_RESULTS.csv"
# VARIABLES <- "RICHNESS_ALL,PD,zPD,PD_P,zPD_P"


## Load input data
cat("Loading Biodiverse results\n")

## Raw indices
cat("..Observed indices\n")
res_r <- fread(INPUTO)

## Z-scores
cat("..Z-scores\n")
res_z <- fread(INPUTZ)

## The first columns should be a gridcell ID
colnames(res_r)[1] <- "H3"
colnames(res_z)[1] <- "H3"

## Remove redundant column
res_r[, Axis_0 := NULL ]
res_z[, Axis_0 := NULL ]

## Rename Z-scores (add `z` prefix)
colnames(res_z)[-1] <- paste0("z", colnames(res_z)[-1])

## Merge the data into a single table
res <- merge(x = res_r, y = res_z, by = "H3", all.x = TRUE)


## If there are multiple variables selected - split them
if(any(grepl(pattern = ",", x = VARIABLES))){
  VARIABLES <- strsplit(x = VARIABLES, split = ",")[[1]]
}


## Check if the selected index is in the tables
colz <- colnames(res)
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


## Prepare spatial polygons
cat("Preparing gridcell polygons\n")
H3_poly <- h3_to_geo_boundary_sf(res$H3)

cat("..Adding diversity estimates to polygons\n")
H3_poly <- cbind(H3_poly, res[, ..VARIABLES])

plot(H3_poly)


## Find plot limits
# boxx <- st_bbox(H3_poly)
# xx <- pretty(c(boxx["xmin"], boxx["xmax"]))
# yy <- pretty(c(boxx["ymin"], boxx["ymax"]))



############################################## Leaflet choropleths

cat("Creating leaflet map\n")

## Set value for the minZoom and maxZoom settings
leaflet(options = leafletOptions(minZoom = 0, maxZoom = 18))


## Function to create a label with single variable
single_label <- function(num, name){

  sprintf(
    "<strong>%s:</strong> %.3g<br/>",
    rep(name, times = length(num)),
    num
    )
}
# single_label(H3_poly$PD, name = "PD")

## Create labels for all variables
cat("..Generating polygon labels\n")
labels <- alply(.data = VARIABLES, .margins = 1, .fun = function(v){
  single_label(num = H3_poly[[ v ]], name = v)
  })

## Concatenate labels from all variables
labels <- do.call(paste0, labels)

## Mark labels as HTML
labels <- labels %>% lapply(htmltools::HTML)


## Create a map, add OpenStreetMap map tiles as background
cat("..Building basemap\n")
m <- leaflet() %>% addTiles()


## Color palette
# pal <- colorQuantile(palette = "RdYlBu", domain = H3_poly$PD, n = 5, reverse = TRUE)
#
# bins <- c(0, 10, 20, 50, 100, 200, 500, 1000, Inf)
# pal <- colorBin("YlOrRd", domain = H3_poly$PD, bins = 7, reverse = TRUE) # , bins = bins)
#
# pal <- colorNumeric(palette = "YlGnBu", domain = H3_poly$PD, reverse = TRUE)



cat("..Preparing color palettes\n")

## Function to create color palette
gen_color_palette <- function(x, type = "quantile", col = "RdYlBu", nbins = 5, rev = TRUE){

  ## Bin numeric data via the quantile function
  if(type %in% "quantile"){
    pal <- colorQuantile(palette = col, domain = x, n = nbins, reverse = rev, na.color = "#808080")
  }

  ## Equal-interval binning (based on `cut` function)
  if(type %in% "equal"){
    pal <- colorBin(palette = col, domain = x, bins = nbins, reverse = rev, na.color = "#808080")
  }

  ## Simple linear mapping from continuous numeric data to an interpolated palette
  if(type %in% "continuous"){
    pal <- colorNumeric(palette = col, domain = x, reverse = rev, na.color = "#808080")
  }

  return(pal)
}
# e.g., gen_color_palette(1:10, type = "quantile", nbins = 5)
#       gen_color_palette(1:10, type = "equal",    nbins = 5)
#       gen_color_palette(1:10, type = "continuous")

## Make color palette for all variables
pals <- alply(.data = VARIABLES, .margins = 1,
  .fun = function(v, ...){ gen_color_palette(x = H3_poly[[ v ]], ...) }, 
  type = PALETTE, col = COLOR, nbins = BINS, rev = TRUE)

names(pals) <- VARIABLES






cat("..Adding polygons\n")

# ## Add PD polygons to the map
# m <- m %>% 
#   addPolygons(data = H3_poly,
#     fillColor = ~pal(PD),
#     group = "PD",
#     opacity = 0.8,
#     fillOpacity = 0.8,
#     weight = 0.3, color = "white", dashArray = "1",
#     highlightOptions = highlightOptions(
#       weight = 2, color = "#777777", dashArray = "1",
#       opacity = 0.8, bringToFront = TRUE),
#     label = labels,
#     labelOptions = labelOptions(
#       style = list("font-weight" = "normal", padding = "3px 8px"),
#       textsize = "10px", direction = "auto")
#     ) %>%
#   addLegend("bottomright", pal = pal, values = H3_poly$PD,
#     title = "PD", group = "PD",  opacity = 1
#     )


for(v in VARIABLES){

  cat("... ", v, "\n")

  m <- m %>% 
    addPolygons(data = H3_poly,
      fillColor = ~ pals[[v]]( H3_poly[[v]] ),
      group = v,
      opacity = 0.8,
      fillOpacity = 0.8,
      weight = 0.3, color = "white", dashArray = "1",
      highlightOptions = highlightOptions(
        weight = 2, color = "#777777", dashArray = "1",
        opacity = 0.8, bringToFront = TRUE),
      label = labels,
      labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "10px", direction = "auto")
      ) %>%
    addLegend("bottomright", pal = pals[[v]], values = H3_poly[[v]],
      title = v, group = v,  opacity = 1
      )

}
rm(v)




## Add variable selector
cat("..Adding variable selector\n")
m <- m %>%
  addLayersControl(
    overlayGroups = c(VARIABLES),
    options = layersControlOptions(collapsed = FALSE)
    )

## Hide all vars except the first one
cat("..Hiding variables\n")
m <- m %>% 
  hideGroup(VARIABLES[-1])




cat("..Exporting the results\n")

## Use mapdeck as the rendering platform instead of leaflet
# mapviewOptions(platform = "mapdeck")

## Save a leaflet map as .html
mapshot(
  m,
  url = OUTPUT,
  file = NULL,
  selfcontained = TRUE
)

# remove_controls = c("zoomControl", "layersControl",
#    "homeButton", "scaleBar", "drawToolbar", "easyButton")


# library(htmlwidgets)
# saveWidget(m, file="m.html")


cat("Plotting finished\n")


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
