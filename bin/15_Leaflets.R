#!/usr/bin/env Rscript

## Visualization of Biodiverse results using intaractive maps (Leaflet-based)

cat("Leaflet-based visualization of Biodiverse results\n")
cat("Script name: 15_Leaflets.R\n")

## Usage example:
# ./15_Leaflets.R \
#   --observed "RND_SPATIAL_RESULTS.csv" \
#   --sesscores "RND_rand--z_scores--SPATIAL_RESULTS.csv" \
#   --sigscores "RND_rand--SPATIAL_RESULTS.csv" \
#   --reccounts "Record_counts_H3.RData" \
#   --variables "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE,PE_WE,SES_PE_WE,CANAPE,Redundancy" \
#   --palette "quantile" \
#   --color "RdYlBu" \
#   --bins 5 \
#   --output "03.Plots/Choropleth.html"

## For a list of indices available in Biodiverse, see:
# https://github.com/shawnlaffan/biodiverse/wiki/IndicesDevVersion

## For Standardized-Effect-Size-based variables, add `SES_` prefix (e.g., `SES_PD`)

## For CANAPE (categorical analysis of neo- and paleoendemism; Mishler et al., 2014), add `CANAPE` to the variables list
# CANAPE is able to distinguish different types of centres of endemism, and can thus give insights 
# into different evolutionary and ecological processes that may be responsible for these patterns. 
# - The centres of paleo-endemism indicate places where there are over-representation of long branches that are rare across the landscape.
# - The centres of neo-endemism indicate an area where there is an over-representation of short branches that are rare on the landscape.
# - Mixture of both paleo-endemism and neo-endemism
# - Super-endemic sites

## H3 grid cell index is represented as a 15-character hexadecimal string (e.g., `830021fffffffff`)

## Output:
# - Choropleth.html = Leaflet-based interactive map
# - Biodiverse_results_merged.txt = Diversity metrics (per grid cell) in a tabular format
# - Diversity_estimates.gpkg = Polygons with diversity metrics in GeoPackage format

## TO DO:
# - handle NA values (e.g., if Richness == 1)


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-r", "--observed"),  action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - observed indices"),
  make_option(c("-s", "--sesscores"), action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - SES-scores"),
  make_option(c("-q", "--sigscores"), action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - Randomization p-values"),
  make_option(c("-n", "--reccounts"), action="store", default=NA, type='character', help="File with the total number of (filtered) records per H3 cell"),
  make_option(c("--resolution"),      action="store", default=4L, type='integer', help="Spatial resolution of the H3 Geospatial Indexing System"),
  make_option(c("-v", "--variables"), action="store", default="RICHNESS_ALL,PD,SES_PD,PD_P,SES_PD_P", type='character', help="Diversity variables to plot (comma-separated entries)"),
  make_option(c("--canapesuper"),     action="store", default=TRUE, type='logical', help="Include the `superendemism` class in CANAPE results (default, FALSE)"),
  make_option(c("-p", "--palette"), action="store", default="quantile", type='character', help="Color palette type"),
  make_option(c("-c", "--color"), action="store", default="RdYlBu", type='character', help="Color gradient scheme for the diversity indices (except for SES, CANAPE, and redundancy metrics)"),
  make_option(c("-b", "--bins"), action="store", default=5L, type='integer', help="Number of color bins for quantile palette"),
  make_option(c("--colorses"), action="store", default="threat", type='character', help="Color scheme for standardized effect sizes, SES (default, `threat`; alternative - `hotspots`"),
  make_option(c("-j", "--redundancy"), action="store", default=0, type='double', help="Redundancy threshold for hiding the grid cells with low number of records (disabled by default)"),
  make_option(c("--shortid"), action="store", default=TRUE, type='logical', help="Shorten H3 index name of grid cell labels on the map"),
  make_option(c("--antimeridianfix"), action="store", default=TRUE, type='logical', help="Fix H3 polygons that cross the antimeridian"),
  # make_option(c("-t", "--threads"), action="store", default=1L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default="Choropleth.html", type='character', help="Output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$observed)){
  stop("Input file with observed PD indices is not specified.\n")
}
if(is.na(opt$sesscores)){
  stop("Input file with SES-scores is not specified.\n")
}
if(is.na(opt$sigscores)){
  stop("Input file with radnomization-based p-values is not specified.\n")
}
if(is.na(opt$reccounts)){
  stop("Input file with the total number of GBIF-records per H3 cell is not specified.\n")
}
if(is.na(opt$output)){
  stop("Output file is not specified.\n")
}


## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
INPUTR      <- opt$observed       # observed results (raw index values)
INPUTS      <- opt$sesscores      # standardized index values (SES)
INPUTP      <- opt$sigscores      # randomisations for each index in SPATIAL_RESULTS
NRECORDS    <- opt$reccounts      # total number of GBIF records per H3-cell
VARIABLES   <- opt$variables
RESOLUTION  <- as.integer(opt$resolution)
CANAPESUPER <- as.logical( opt$canapesuper)
PALETTE     <- opt$palette
COLOR       <- opt$color
BINS        <- as.numeric( opt$bins )
COLORSES    <- opt$colorses
REDUNDANCYTRSH <- as.numeric(to_na( opt$redundancy ))
SHORTID <- as.logical( opt$shortid )
ANTIFIX <- as.logical( opt$antimeridianfix )
OUTPUT  <- opt$output

## Check the redundancy range
if(!is.na(REDUNDANCYTRSH)){
  if(! (0 <= REDUNDANCYTRSH & REDUNDANCYTRSH <= 1) ){
    stop("Redundacy threshold should be in the [0,1] range.\n")
  }
}

## Log assigned variables
cat(paste("Input file (observed indices): ",  INPUTR,     "\n", sep=""))
cat(paste("Input file (SES-scores): ",        INPUTS,     "\n", sep=""))
cat(paste("Input file (p-values): ",          INPUTP,     "\n", sep=""))
cat(paste("Input file (number of records): ", NRECORDS,   "\n", sep=""))
cat(paste("Spatial resolution: ",             RESOLUTION, "\n", sep=""))
cat(paste("Indices to plot: ",                VARIABLES,  "\n", sep=""))
cat(paste("CANAPE superendemism enabled: ",   CANAPESUPER,"\n", sep=""))
cat(paste("Color palette type: ",             PALETTE,    "\n", sep=""))
cat(paste("Color gradient scheme: ",          COLOR,      "\n", sep=""))
cat(paste("Number of color bins: ",           BINS,       "\n", sep=""))
cat(paste("SES color palette: ",              COLORSES,   "\n", sep=""))
cat(paste("Redundancy threshold: ",           REDUNDANCYTRSH, "\n", sep=""))
cat(paste("Display short H3 index names: ",   SHORTID,    "\n", sep=""))
cat(paste("Antimeridian fix: ",               ANTIFIX,    "\n", sep=""))
cat(paste("Output file: ",                    OUTPUT,     "\n", sep=""))

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

## Parameters for debugging
# INPUTR <- "RND_SPATIAL_RESULTS.csv"
# INPUTS <- "RND_rand--z_scores--SPATIAL_RESULTS.csv"
# INPUTP <- "RND_rand--SPATIAL_RESULTS.csv"    # "RND_rand--SPATIAL_RESULTS.csv"
# NRECORDS <- "Record_counts_H3.RData"
# VARIABLES <- "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE,PE_WE,SES_PE_WE,CANAPE,Redundancy"
# REDUNDANCYTRSH <- 0.7
# PALETTE <- "quantile"
# COLOR <- "RdYlBu"
# BINS <- 5
# SHORTID <- TRUE
# ANTIFIX <- TRUE

## Load input data
cat("Loading Biodiverse results\n")

## Raw indices
cat("..Observed indices\n")
res_r <- fread(INPUTR)

## SES-scores
cat("..SES-scores\n")
res_s <- fread(INPUTS)

## P-values
cat("..P-values\n")
res_p <- fread(INPUTP)


## Test if the first column cotains valid H3 IDs
if( h3_is_valid(res_r[[1,1]]) ){
  colnames(res_r)[1] <- "H3"
  colnames(res_z)[1] <- "H3"
  colnames(res_p)[1] <- "H3"
} else {
  ## Get H3 IDs for grid cells
  cat("H3 index was not found in the data\n")
  cat("..Indexing geo-coordinates\n")
  res_r[ , H3 := h3::geo_to_h3(res_r[, .(Axis_0, Axis_1)], res = RESOLUTION) ]
  res_s[ , H3 := h3::geo_to_h3(res_s[, .(Axis_0, Axis_1)], res = RESOLUTION) ]
  res_p[ , H3 := h3::geo_to_h3(res_p[, .(Axis_0, Axis_1)], res = RESOLUTION) ]
}

## Remove redundant columns
cat("Removing redundant columns\n")
res_r[, c("ELEMENT", "Axis_0", "Axis_1") := NULL ]
res_s[, c("ELEMENT", "Axis_0", "Axis_1") := NULL ]
res_p[, c("ELEMENT", "Axis_0", "Axis_1") := NULL ]


## Rename SES-scores (add `SES_` prefix)
setnames(res_s,
  old = colnames(res_s)[ ! colnames(res_s) %in% "H3" ],
  new = paste0("SES_", colnames(res_s)[ ! colnames(res_s) %in% "H3" ]))

## Merge the data into a single table
cat("Merging data into a single table\n")
res <- merge(x = res_r, y = res_s, by = "H3", all.x = TRUE)

## Clean up
rm(res_r, res_s)

## If there are multiple variables selected - split them
if(any(grepl(pattern = ",", x = VARIABLES))){
  VARIABLES <- strsplit(x = VARIABLES, split = ",")[[1]]
}

## ?? Remove "monomorphic" variables
## They could cause an error with color gradients and binning

## CANAPE (categorical analysis of neo- and paleoendemism; Mishler et al., 2014)
if("CANAPE" %in% VARIABLES){
  cat("Preparing data for CANAPE analysis\n")
  do_CANAPE <- TRUE

  required_vars <- c(
    "P_PHYLO_RPD2", 
    "P_PD_P", "P_PE_WE_P", "P_PD_P_per_taxon",
    "P_PHYLO_RPE2",
    "P_PHYLO_RPE_NULL2")

  ## Check if we have all the required variables in the data
  if(any(! required_vars %in% colnames(res_p) )){
    cat("..WARNING: some variables required for the CANAPE analysis are missing from the data!\n")
    cat(paste(required_vars[ ! required_vars %in% colnames(res_p) ], collapse = ", "), "\n")
    cat("..Please add `calc_phylo_rpd2,calc_phylo_rpe2` to the `indices` parameter of the pipeline\n")
    cat("..Skipping the CANAPE analysis\n")
    do_CANAPE <- FALSE
  } else {
    
    cat("..Running the CANAPE analysis\n")

    ## Significance-testing and thresholds are based on the code by Nunzio Knerr and Shawn Laffan
    # https://github.com/NunzioKnerr/biodiverse_pipeline/blob/0e3eb237870a3a72a0be9939062123eee488e960/R_release/load_biodiverse_results_and_report_on_CANAPE_by_taxa.R#L38
    ## Two-tailed test for RPD
    sig1 <- function(x){
      if (x >= 0.99) {
        return("VeryHighlySig")
      } else if (x >= 0.975){
        return ("HighlySig")
      } else if (x <= 0.01){
        return ("VerySigLow")
      } else if (x <= 0.025){
        return ("SigLow")
      } else {
        return("NotSig")
      }
    }

    ## Two-pass test for RPE
    sig2 <- function(x, y, z, with_super = TRUE){
      # x = P_PE_WE_P; y = P_PHYLO_RPE_NULL2; z = P_PHYLO_RPE2

      if (is.na(x)) { x = 0   }
      if (is.na(y)) { y = 0   }
      if (is.na(z)) { z = 0.5 }

      if(with_super == TRUE){
        if (x < 0.95 & y < 0.95) {
          return("NotSignificant")
        } else if (z <= 0.025){
          return ("Neo_endemism")
        } else if (z >= 0.975){
          return ("Paleo_endemism")
        } else if (x >= 0.99 & y >= 0.99){
          return ("Super_endemism")
        } else {
          return("Mixed_endemism")
        }
      }

      if(with_super == FALSE){
        if (x <= 0.95 & y <= 0.95) {
           return ("NotSignificant")
         } else if (z < 0.025) {
           return ("Neo_endemism")         # RPE < 0.025
         } else if (z > 0.975) {
           return ("Paleo_endemism")       # RPE > 0.975
         } else {
           return ("Mixed_endemism")
         }
      }
    } # end of sig2

    ## Vectorize the functions
    sig1 <- Vectorize(sig1)
    sig2 <- Vectorize(sig2)

    ## Step 1
    # cat("...Step 1 - Parsing significance tests\n")
    # canape_data <- data.table(
    #   H3 = res_p[[ "H3" ]],
    #   P_PHYLO_RPD1_SIG = sig1( res_p[[ "P_PHYLO_RPD1" ]] ),
    #   P_PHYLO_RPD2_SIG = sig1( res_p[[ "P_PHYLO_RPD2" ]] ),
    #   P_PD_P_SIG = sig1( res_p[[ "P_PD_P" ]] ),
    #   P_PE_WE_P_SIG = sig1( res_p[[ "P_PE_WE_P" ]] ),
    #   P_PD_P_per_taxon_SIG = sig1( res_p[[ "P_PD_P_per_taxon" ]] ),
    #   P_PHYLO_RPE2_ONE_STEP_SIG = sig1( res_p[[ "P_PHYLO_RPE2" ]] )
    #   )

    ## Step 2
    # cat("...Step 2 - Inferring endemism type\n")

    ## Null model of PD evenly distributed across terminals, but with the same range per terminal and where ancestral nodes are of zero length
    ## P_PHYLO_RPE1_SIG
    # canape_data$CANAPE <- sig2(
    #   res_p[[ "P_PE_WE_P" ]],
    #   res_p[[ "P_PHYLO_RPE_NULL1" ]],
    #   res_p[[ "P_PHYLO_RPE1" ]])

    ## Null model where PE is calculated using a tree where all branches are of equal length
    ## P_PHYLO_RPE2_SIG
    # canape_data$CANAPE <- sig2(
    #   res_p[[ "P_PE_WE_P" ]],
    #   res_p[[ "P_PHYLO_RPE_NULL2" ]],
    #   res_p[[ "P_PHYLO_RPE2" ]])

    cat("...Inferring endemism type\n")
  
    canape_data <- data.table(
      H3 = res_p[[ "H3" ]],
      CANAPE = sig2(
          res_p[[ "P_PE_WE_P" ]],
          res_p[[ "P_PHYLO_RPE_NULL2" ]],
          res_p[[ "P_PHYLO_RPE2" ]],
          with_super = CANAPESUPER)
      )

    canape_data$CANAPE <- factor(canape_data$CANAPE,
      levels = c("Neo_endemism", "Paleo_endemism", "NotSignificant", "Mixed_endemism", "Super_endemism"))

    ## Add CANAPE to the main table
    res <- merge(x = res, y = canape_data, by = "H3", all.x = TRUE)

  } # end of CANAPE

} else {
  cat("Skipping the CANAPE analysis\n")
  do_CANAPE <- FALSE
}


### Estimate sampling redundancy for each cell (how well the are is sampled) 
###  = 1 - (Richness / Number of specimens)
### (see Mishler et al., 2020; DOI: 10.1111/jse.12590) 
cat("Estimating sampling redundancy\n")
if(!"RICHNESS_ALL" %in% colnames(res)){
  cat("..WARNING: species richness is missing from the estimated indices!\n")
  cat("..Skipping redundancy estimation\n")

  ## Remove Redundancy variable from the list (if it's there)
  VARIABLES <- VARIABLES[ ! VARIABLES %in% "Redundancy" ]

} else {

  ## Loading number of records
  cat("..Loading a file with the total number of GBIF-records per H3-cell\n")
  NRECORDS <- readRDS(NRECORDS)

  ## Add N to the main table
  cat("..Adding N records to the main table\n")
  res <- merge(x = res, y = NRECORDS[, .(H3, NumRecords)],
               by = "H3", all.x = TRUE)
  
  ## Calc the redundancy
  cat("..Estimating the index\n")
  res[, Redundancy := ( 1 - (RICHNESS_ALL / NumRecords) )]

  ## Add the index to the VARIABLES
  # VARIABLES <- c(VARIABLES, "Redundancy")
}


## Check if the selected index is in the tables
colz <- colnames(res)
if(any(!VARIABLES %in% colz)){
  cat("Some of the selected indices are not present in tables with results!\n")
  cat("Please check the spelling of index names or enable their estimation in Biodiverse.\n")
  missing <- VARIABLES[ ! VARIABLES %in% colz ]
  cat("Indices missing: ", paste(missing, collapse = ", "), "\n")

  ## Exclude missing indices
  VARIABLES <- VARIABLES[ ! VARIABLES %in% missing ]
}

## Check the data
if(length(VARIABLES) == 0){
  stop("None of the selected indices were found in the results! Nothing to plot.\n")
}

## Remove gridcell with the redundancy index below the specified threshold
if(REDUNDANCYTRSH > 0 & "Redundancy" %in% colnames(res)){
  cat("Removing grid cells with low redundancy index\n")
  cat("..The specified threshold is ", REDUNDANCYTRSH, "\n")
  cat("..There are ", sum(res$Redundancy <= REDUNDANCYTRSH), " grid cells to be removed\n")
  cat("..And ", sum(res$Redundancy > REDUNDANCYTRSH), " grid cells will be preserved\n")

  if(any(res$Redundancy <= REDUNDANCYTRSH)){
    res <- res[ Redundancy > REDUNDANCYTRSH ]
  }
}

## Check the data
if(!nrow(res) > 0){
  stop("There are no grid cells to display!\n")
}


cat("Preparing data table for export\n")

## Add grid cell coordinates
cat(".. Adding geo-coordinates for grid cell centers\n")
res[, c("Latitude", "Longitude") := as.data.table(h3::h3_to_geo(H3)) ]

## Reorder columns
cat(".. Reordering columns\n")
setcolorder(res, c("H3", "Latitude", "Longitude"))

## Export the data
cat("Exporting the data table\n")
fwrite(x = res, file = "Biodiverse_results_merged.txt", sep = "\t")


## Prepare spatial polygons
cat("Preparing gridcell polygons\n")
H3_poly <- h3_to_geo_boundary_sf(res$H3)

cat("..Adding diversity estimates to polygons\n")

if("NumRecords" %in% colnames(res)){
  vars <- c(VARIABLES, "NumRecords")
  H3_poly <- cbind(H3_poly, res[, ..vars])
} else {
  H3_poly <- cbind(H3_poly, res[, ..VARIABLES])
}

## Assign rownames as H3 grid cell IDs
cat("..Adding H3 cell names\n")
if(SHORTID == TRUE){
  ## Truncate fff-tail of index name (e.g., `830021fffffffff` -> `830021f`)
  rownames(H3_poly) <- gsub(pattern = "f+$", replacement = "f", x = res$H3)
} else {
  rownames(H3_poly) <- res$H3
}

cat("..Exporting polygons with divsity estimates in GeoPackage format\n")
## (using sf + GPKG driver of GDAL library)
st_write(
  obj = H3_poly,
  dsn = "Diversity_estimates.gpkg",
  layer = "diversity_estimates")


## Fix H3 polygons that cross the antimeridian by cutting them in two
if(ANTIFIX == TRUE){
  cat("..Fixing antimeridian issue\n")
  H3_poly <- st_wrap_dateline(H3_poly, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
}


# plot(H3_poly)


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

  if( is.numeric(num) ){
    res <- sprintf(
      "<strong>%s:</strong> %.3g<br/>",
      rep(name, times = length(num)),
      num
      )
  } else {
    res <- sprintf(
      "<strong>%s:</strong> %s<br/>",
      rep(name, times = length(num)),
      num
      )
  }
  return(res)
}
# single_label(H3_poly$PD, name = "PD")           # numeric labels
# single_label(H3_poly$CANAPE, name = "CANAPE")   # factor labels

## Create labels for all variables
cat("..Generating polygon labels\n")
labels <- alply(.data = VARIABLES, .margins = 1, .fun = function(v){
  single_label(num = H3_poly[[ v ]], name = v)
  })

## Add number of GBIF-records to the labels
labels <- c(labels, 
  list( single_label(num = H3_poly[[ "NumRecords" ]], name = "Number of records"))
  )

## Add H3 index to the labels
labels <- c(labels, 
  list( single_label(num = rownames(H3_poly), name = "H3 index") )
  )


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
gen_color_palette <- function(x, type = "quantile", col = "RdYlBu", colses = "threat", nbins = 5, rev = TRUE){

  ## Bin numeric data via the quantile function
  if(type %in% "quantile"){

    ## To fix `'breaks' are not unique` error, check the number of potential bins in advance
    testbins <- quantile(x,
      probs = seq(0, 1, length.out = nbins + 1),
      na.rm = TRUE, names = FALSE)
    
    nbins <- length(unique(testbins)) - 1

    pal <- colorQuantile(palette = col, domain = x, n = nbins, reverse = rev, na.color = "#808080")
    attr(pal, which = "newbins") <- nbins
  }

  ## Equal-interval binning (based on `cut` function)
  if(type %in% "equal"){
    pal <- colorBin(palette = col, domain = x, bins = nbins, reverse = rev, na.color = "#808080")
    attr(pal, which = "newbins") <- nbins
  }

  ## Simple linear mapping from continuous numeric data to an interpolated palette
  if(type %in% "continuous"){
    pal <- colorNumeric(palette = col, domain = x, reverse = rev, na.color = "#808080")
    attr(pal, which = "newbins") <- 999
  }

  ## SES-score mapping (symmetric around zero)
  if(type %in% "ses"){

    ## Hotspot-type colors (high SES values = red, low = blue)
    if(colses %in% "hotspots"){
    ses_colors <- c(
      HighlyNegative = "#27408B",
      Negative       = "#4876FF",
      NotSignificant = "#FAFAD2",
      Positive       = "#FF0000",
      HighlyPositive = "#8B0000"
      )
    }

    ## Threat-type colors (high SES values = blue, low = red), as in Mishler et al., 2014
    if(colses %in% "threat"){
    ses_colors <- c(
      HighlyNegative = "#8B0000",
      Negative       = "#FF0000",
      NotSignificant = "#FAFAD2",
      Positive       = "#4876FF",
      HighlyPositive = "#27408B"
      )
    }

    ses_bins <-c(-1000, -2.58, -1.96, 1.96, 2.58, 1000)

    pal <- colorBin(palette = ses_colors, bins = ses_bins,
      na.color = "#808080", domain = c(-1000, 1000))
    
    attr(pal, which = "newbins") <- 5
    # plot(-4:4, col = pal(-4:4), cex = 2, pch = 16) # test
  }

  ## CANAPE-style mapping (neo / paleo / mixed / super)
  if(type %in% "canape"){

    canape_colors <- c(
      Neo_endemism   = "#FF0000",
      Paleo_endemism = "#4876FF",
      NotSignificant = "#FAFAD2",
      Mixed_endemism = "#CB7FFF",
      Super_endemism = "#9D00FF"
      )

    pal <- colorFactor(
      palette = canape_colors,
      levels = names(canape_colors),
      na.color = "#808080")

    attr(pal, which = "newbins") <- 5
    # plot(1:length(canape_colors), col = pal(names(canape_colors)), cex = 2, pch = 16) # test
  }

  ## Redundancy-style mapping (yellow-brown, [0,1])
  if(type %in% "redundancy"){
    pal <- colorNumeric(palette = "YlOrBr", domain = seq(0, 1, 0.01), reverse = FALSE, na.color = "#808080")
    attr(pal, which = "newbins") <- 999
    # plot(seq(0,1,0.05), col = pal(seq(0,1,0.05)), cex = 2, pch = 16) # test
  }

  return(pal)
}
# e.g., gen_color_palette(1:10, type = "quantile", nbins = 5)
#       gen_color_palette(1:10, type = "equal",    nbins = 5)
#       gen_color_palette(1:10, type = "continuous")

## Make color palette for all variables
pals <- list()

VARIABLES_ses <- grep(pattern = "^SES_", x = VARIABLES, value = TRUE)
VARIABLES_raw <- grep(pattern = "^SES_", x = VARIABLES, value = TRUE, invert = TRUE)
VARIABLES_raw <- VARIABLES_raw[ ! VARIABLES_raw %in% c("CANAPE", "Redundancy") ]

## Colors for "raw" variables
if(length(VARIABLES_raw) > 0){

  pals_raw <- alply(.data = VARIABLES_raw, .margins = 1,
    .fun = function(v, ...){ gen_color_palette(x = H3_poly[[ v ]], ...) }, 
    type = PALETTE, col = COLOR, nbins = BINS, rev = TRUE)

  names(pals_raw) <- VARIABLES_raw
  pals <- c(pals, pals_raw)
  rm(pals_raw)
}


## Colors of SES-scores
if(length(VARIABLES_ses) > 0){

  pals_ses <- alply(.data = VARIABLES_ses, .margins = 1,
    .fun = function(v, ...){ gen_color_palette(x = H3_poly[[ v ]], ...) }, 
    type = "ses", colses = COLORSES)

  names(pals_ses) <- VARIABLES_ses
  pals <- c(pals, pals_ses)
  rm(pals_ses)
}


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



## Shortcut to add polygons with legend to the map
add_polygons_with_legend <- function(m, v, pal, ses_labels = FALSE){
  res <- m %>% 
      addPolygons(data = H3_poly,
        fillColor = ~ pal( H3_poly[[v]] ),
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
      addLegend("bottomright", pal = pal, values = H3_poly[[v]],
        title = v, group = v,  opacity = 1)

  ## Relable upper bounds of effect sizes
  ## (unfortunately, addLegend(labels =...) does not work with `pal` argument)
  if(ses_labels == TRUE){
    
    newlabs <- c("< -2.58", "< -1.96", "-1.96 - 1.96", "> 1.96", "> 2.58")
    newlabs <- gsub(pattern = " - ", replacement = " &ndash; ", x = newlabs)

    slotid <- length(res$x$calls)
    if(
      "labels" %in% names(res$x$calls[[ slotid ]]$args[[1]]) &                # `labels` slot is present
      length(res$x$calls[[ slotid ]]$args[[1]]$labels) == length(newlabs)     # and has the same number of categories
      ){
      res$x$calls[[ slotid ]]$args[[1]]$labels <- newlabs
    }
  }

  return(res)
}
## Example:
# add_polygons_with_legend(m, "PD", pal = pals[[ "PD" ]])


## Loop throug all variables
## Except "CANAPE" (it has a categorical color scheme) and "Redundancy"
for(v in VARIABLES[ ! VARIABLES %in% c("CANAPE", "Redundancy") ]){

  cat("... ", v, "\n")

  if(v %in% VARIABLES_raw & PALETTE %in% "quantile"){
    if(attr(pals[[v]], "newbins") != BINS){
    cat(".... number of bins was adjusted to ", attr(pals[[v]], "newbins"), "\n")
    }
  }

  ## SES lables?
  if(v %in% VARIABLES_ses){
    ses_lab <- TRUE
  } else {
    ses_lab <- FALSE
  }

  tmp <- try( add_polygons_with_legend(m = m, v = v, pal = pals[[v]], ses_labels = ses_lab) )

  ## Quantile palette may fail
  if(("try-error" %in% class(tmp) | attr(pals[[v]], "newbins") == 1) & PALETTE %in% "quantile"){
   cat(".... Warning: quantile palette failed, trying continuous palette.\n") 

   ## Generate new color palette (continuous)
   tmppal <- gen_color_palette(x = H3_poly[[ v ]],
    type = "continuous", col = COLOR, nbins = BINS, rev = TRUE)

   tmp <- add_polygons_with_legend(m = m, v = v, pal = tmppal)

   rm(tmppal)
  } # end of `try-error`

  m <- tmp
  rm(tmp, ses_lab)

}   # end of loop
rm(v)


## Add CANAPE to the plot
if("CANAPE" %in% VARIABLES){

  cat("... CANAPE\n")

  canape_pal <- gen_color_palette(
    x = H3_poly[[ "CANAPE" ]],
    type = "canape")

  m <- m %>% 
    addPolygons(data = H3_poly,
      fillColor = ~ canape_pal( H3_poly[[ "CANAPE" ]] ),
      group = "CANAPE",
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
    addLegend("bottomright", pal = canape_pal, values = H3_poly[[ "CANAPE" ]],
      title = "CANAPE", group = "CANAPE",  opacity = 1)

} # end of CANAPE


## Add Redundancy index to the plot
if("Redundancy" %in% VARIABLES){

  cat("... Redundancy\n")

  redundancy_pal <- gen_color_palette(
    x = H3_poly[[ "Redundancy" ]],
    type = "redundancy")

  m <- m %>% 
    addPolygons(data = H3_poly,
      fillColor = ~ redundancy_pal( H3_poly[[ "Redundancy" ]] ),
      group = "Redundancy",
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
    addLegend("bottomright", pal = redundancy_pal, values = H3_poly[[ "Redundancy" ]],
      title = "Sampling redundancy", group = "Redundancy",  opacity = 1)

}


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

cat("...Exporting Choropleth map in HTML\n")

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

cat("...Exporting Leaflet object\n")

attr(m, which = "VARIABLES") <- VARIABLES

saveRDS(object = m, file = "Leaflet_object.RData", compress = "xz")


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
