#!/usr/bin/env Rscript

## Filter GBIF occurrence data

cat("Filtering GBIF species occurrence data\n")
cat("Script name: 10_Filter_occurrences.R\n")

## Input data could be located in the cloud.
## E.g., to query GBIF AWS snapshot use:
##  --input "s3://gbif-open-data-eu-central-1/occurrence/2022-01-01/occurrence.parquet"

## Usage:
# ./10_Filter_occurrences.R \
#    --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet" \
#    --family "Fabaceae" \
#    --country "AU" \
#    --latmin -55.3228175 \
#    --latmax -9.0882278 \
#    --lonmin 72.2460938 \
#    --lonmax 168.2249543 \
#    --minyear 1945 \
#    --roundcoords 2 \
#    --threads 10 \
#    --noccurrences 30 \
#    --output "Fabaceae_in_AU"

## Note about rounding of coordinates (`--roundcoords 2`):
## A value in decimal degrees to an accuracy of 2 decimal places is accurate to ~1.11 km at the equator

# TO DO:
# - add country name validation
# - check if the input folder exists (except of 's3')

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

  ## Polygons (to find a bounding box if no coordinates are provided)
  make_option(c("-g", "--polygon"), action="store", default=NA, type='character', help="Custom area of interest (a file with polygons in GeoPackage format)"),
  make_option(c("-w", "--wgsrpd"),  action="store", default=NA, type='character', help="Path to the World Geographical Scheme for Recording Plant Distributions data (polygons in sf-format)"),
  make_option(c("-x", "--regions"), action="store", default=NA, type='character', help="Comma-separated list of WGSRPD regions"),

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

  make_option(c("--roundcoords"), action="store", default=2L, type='integer', help="Round spatial coordinates to the N decimal places, to reduce the dataset size (default, 2). To disable, set to a negative value."),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-n", "--noccurrences"), action="store", default=30, type='double', help="Occurrence threshold (used for parquet partitioning)"),
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

POLYGON       <- to_na( opt$polygon )
WGSRPD        <- to_na( opt$wgsrpd )
WGSRPDREGIONS <- to_na( opt$regions )

MINYEAR <- as.numeric(to_na( opt$minyear) )
MAXYEAR <- as.numeric(to_na( opt$maxyear) )

EXTINCT      <- to_na( opt$noextinct )
EXCLUDEHUMAN <- as.logical( opt$excludehuman )

BASISINCL <- to_na( opt$basisofrecordinclude )
BASISEXCL <- to_na( opt$basisofrecordexclude )

ROUNDCOORDS  <- as.numeric( opt$roundcoords )
CPUTHREADS   <- as.numeric( opt$threads )
OCCURRENCES  <- as.numeric( opt$noccurrences )
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

cat(paste("Custom polygons: ",  POLYGON,       "\n", sep=""))
cat(paste("WGSRPD data: ",      WGSRPD,        "\n", sep=""))
cat(paste("WGSRPD regions: ",   WGSRPDREGIONS, "\n", sep=""))

cat(paste("Basis of record to include: ", BASISINCL, "\n", sep=""))
cat(paste("Basis of record to exclude: ", BASISEXCL, "\n", sep=""))
cat(paste("Minimum year of occurrence: ", MINYEAR, "\n", sep=""))
cat(paste("Maximum year of occurrence: ", MAXYEAR, "\n", sep=""))
cat(paste("List of extict species: ",     EXTINCT, "\n", sep=""))
cat(paste("Exclusion of human records: ", EXCLUDEHUMAN, "\n", sep=""))

cat(paste("Round coordinates: ", ROUNDCOORDS, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Occurrence threshold for parquet partitioning: ", OCCURRENCES, "\n", sep=""))
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
  colnames(SPECIESKEYS)[1] <- "specieskey"
  SPECIESKEYS <- unique(na.omit(SPECIESKEYS))
  cat("Specieskey list loaded. Number of records: ", nrow(SPECIESKEYS), "\n")
}


############################################## Main pipeline


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


## Country-based filtering
if(!is.na(COUNTRY)){
  cat("..Filtering by Country\n")
  COUNTRY <- strsplit(x = COUNTRY, split = ",")[[1]]
  dsf <- dsf %>% filter(countrycode %in% COUNTRY)
}



### Filtering by coordinates

## If no coordinates specified,
## but filtering by custom polygons or WGSRPD is required,
## find a coordinate box
any_coordinates <- sum(!is.na(LATMIN), !is.na(LATMAX), !is.na(LONMIN), !is.na(LONMAX))
any_polygons <- sum(!is.na(POLYGON), !is.na(WGSRPDREGIONS))
if(any_coordinates == 0 & any_polygons > 0){

  cat("Coordinate box is not specified, but spatial polygons are provided for filtering\n")
  cat("Infering coordinate box from polygons\n")

  ## User-supplied polygons
  if(!is.na(POLYGON)){
    cat("..Loading GeoPackge file\n")
    POLYGON <- read_sf(POLYGON)
    poly_bbox <- st_bbox(POLYGON)  # x = longitude, y = latitude
    poly_latmin <- poly_bbox$ymin
    poly_latmax <- poly_bbox$ymax
    poly_lonmin <- poly_bbox$xmin
    poly_lonmax <- poly_bbox$xmax
  } else {
    poly_latmin <- poly_latmax <- poly_lonmin <- poly_lonmax <- NA
  }

  ## WGSRPD polygons
  if(!is.na(WGSRPD) & !is.na(WGSRPDREGIONS)){
    cat("..Loading WGSRPD polygons\n")
    WGSRPD <- readRDS(WGSRPD)
    WGSRPDREGIONS <- strsplit(x = WGSRPDREGIONS, split = ",")[[1]]
    WGSRPD <- WGSRPD[ which(WGSRPD$LevelName %in% WGSRPDREGIONS),  ]
    wgsrpd_bbox <- st_bbox(WGSRPD)
    wgsrpd_latmin <- wgsrpd_bbox$ymin
    wgsrpd_latmax <- wgsrpd_bbox$ymax
    wgsrpd_lonmin <- wgsrpd_bbox$xmin
    wgsrpd_lonmax <- wgsrpd_bbox$xmax
  } else {
    wgsrpd_latmin <- wgsrpd_latmax <- wgsrpd_lonmin <- wgsrpd_lonmax <- NA
  }

  cat("..Finding coordinate bounding box:\n")
  LATMIN <- min(c(poly_latmin, wgsrpd_latmin), na.rm = TRUE)
  LATMAX <- max(c(poly_latmax, wgsrpd_latmax), na.rm = TRUE)
  LONMIN <- min(c(poly_lonmin, wgsrpd_lonmin), na.rm = TRUE)
  LONMAX <- max(c(poly_lonmax, wgsrpd_lonmax), na.rm = TRUE)

  cat(paste("...Minimum latitude from polygons: ",  LATMIN,  "\n", sep = ""))
  cat(paste("...Maximum latitude from polygons: ",  LATMAX,  "\n", sep = ""))
  cat(paste("...Minimum longitude from polygons: ", LONMIN,  "\n", sep = ""))
  cat(paste("...Maximum longitude from polygons: ", LONMAX,  "\n", sep = ""))

}

## Filter by coordinates
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
cat("..Rounding coordinates\n")
if(ROUNDCOORDS >= 0){
  cat("..Rounding coordinates\n")
  dsf <- dsf %>%
    mutate(
      decimallatitude  = round(decimallatitude,  ROUNDCOORDS),
      decimallongitude = round(decimallongitude, ROUNDCOORDS))
}

## Select columns and remove duplicated records
cat("..Column selection\n")
dsf <- dsf %>%
  select(specieskey, species, decimallongitude, decimallatitude) %>%
  distinct()



## Count number of records by species
cat("Counting number of records per species\n")
sp_counts <- dsf %>%
  count(specieskey) %>%
  collect() %>%
  mutate(Partition = case_when(n <= OCCURRENCES ~ "low",
                              n > OCCURRENCES ~ "high"))

smr <- table(sp_counts$Partition)
if("low" %in% names(smr)){
  cat("..Number of species with low number of occurrences: ", smr[["low"]], "\n")
}
if("high" %in% names(smr)){
  cat("..Number of species with high number of occurrences: ", smr[["high"]], "\n")
}

## Export species counts
fwrite(x = sp_counts,
  file = "SpeciesCounts.txt",
  sep = "\t")


## Group by species  -- there would be too many partitions (default arrow max = 1024)
# cat("Groupping by species\n")
# dsf <- dsf %>%
#   group_by(specieskey)


## Add partition ID to the dataset
dsf <- dsf %>%
  left_join(sp_counts[, c("specieskey", "Partition")]) %>%
  group_by(Partition)


## Export species occurrence, partition by species
cat("Exporting filtered occurrence data\n")

# outdir <- paste0(OUTPUT, ".parquet")
outdir <- OUTPUT
if(! outdir %in% "."){ dir.create(outdir, showWarnings = FALSE) }

dir_high <- paste0(outdir, "/Partition=high")
dir_low <- paste0(outdir, "/Partition=low")
dir.create(dir_high, showWarnings = FALSE)
dir.create(dir_low, showWarnings = FALSE)

write_dataset(
  dataset = dsf,
  path = outdir,
  format = "parquet",
  basename_template = "occ_{i}")




##################### Session info

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
