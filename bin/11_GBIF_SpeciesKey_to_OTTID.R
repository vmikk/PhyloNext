#!/usr/bin/env Rscript

## Prepare file with Open Tree species IDs (ott_ids) based on GBIF species keys
## In the output, first column must have the header 'ott_id' and contain OpenTree taxon ids
## Output will be used as input for `induced_synth_subtree_from_csv.py` script

cat("Species-level filtering of occurrence data and spatial binning\n")
cat("Script name: 11_GBIF_SpeciesKey_to_OTTID.R\n")

## Usage:
# ./11_GBIF_SpeciesKet_to_OTTID.R \
#    --input "SpeciesCounts.txt" \
#    --taxgroup "Mammals" \
#    --threads 1

## Input  = SpeciesCounts.txt (based on `10_Record_counts.R` results)
## Output = Species_OTT.csv


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the directory with pre-filtered GBIF occurrences in Parquet format"),
  make_option(c("-g", "--taxgroup"), action="store", default="All life", type='character', help="Taxonomy group in OpenTree"),
  make_option(c("-t", "--threads"), action="store", default=1L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default="Species_OTT.csv", type='character', help="Output file")
  )
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$input)){
  cat("Input is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output file is not specified.\n", file=stderr())
  stop()
}

## Assign variables
INPUT      <- opt$input
TAXGROUP   <- opt$taxgroup
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT     <- opt$output

## Modify tax group argument (no spaces are allowed in arguments)
# if(TAXGROUP %in% "All_life"){ TAXGROUP <- "All life" }
TAXGROUP <- gsub(pattern = "_", replacement = " ", x = TAXGROUP)


## Log assigned variables
cat(paste("Input occurrences: ", INPUT, "\n", sep=""))
cat(paste("Taxonomy group in OpenTree: ", TAXGROUP, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output file: ", OUTPUT, "\n", sep=""))

cat("\n")

############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("rotl")     # Interface to the 'Open Tree of Life' API
load_pckg("rgbif")

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

## Open dataset
cat("Loading input data\n")
datt <- fread(file = INPUT, sep = "\t")

cat("..Number of specieskeys in the data: ", nrow(datt), "\n")

cat("\nParsing GBIF taxonomy with `rgbif` package\n")

## Obtain species names based on GBIF specieskeys
spnames <- alply(.data = datt$specieskey, .margins = 1, .fun = function(x){
  # e.g., x <- 5219243

  ## Lookup details for a specific name in all taxonomies in GBIF
  tx <- try( rgbif::name_usage(key = x, limit = 5) )

  # TO DO - check for fuzzy matches?
  setDT(tx$data)
  tx$data <- tx$data[ key == x ]

  if("try-error" %in% class(tx) | nrow(tx$data) == 0){
    tax <- data.table(kingdom = NA, phylum = NA, class = NA, order = NA, family = NA, genus = NA, species = NA)
  } else {
    tax <- data.table(
      kingdom = tx$data$kingdom[1], 
      phylum  = tx$data$phylum[1], 
      class   = tx$data$class[1], 
      order   = tx$data$order[1], 
      family  = tx$data$family[1], 
      genus   = tx$data$genus[1],
      species = tx$data$species[1]
      )
  }

  res <- data.table(
    specieskey = x,
    tax)

 return(res)
}, .parallel = parall)

spnames <- rbindlist(spnames, fill = TRUE)


cat("..Number of species found: ", length(na.omit(unique(spnames$species))), "\n")
cat("..Number of specieskeys without matches: ", sum(is.na(spnames$species)), "\n")

if(any(is.na(spnames$species))){
  missing_spkeys <- spnames[ is.na(species) ]$specieskey
  cat("..Specieskeys without matches: \n")
  cat(paste(missing_spkeys, collapse = " "), "\n")

  spnames <- spnames[ ! is.na(species) ]
}


cat("\nParsing OTT taxonomy\n")


## Function to split vector into n chunks
chunk <- function(x, n){
  if(n > 1) { res <- split(x, cut(seq_along(x), n, labels = FALSE)) }
  if(n == 1){ res <- list(); res[[1]] <- x }
  return(res)
}

chunk_table <- function(x, n){
  if(n == 1){ res <- list(); res[[1]] <- x }
  if(n > 1) {
    cc <- chunk(1:nrow(x), n = n)
    res <- llply(.data = cc, .fun = function(r){ x[ r, ] })
  }
  return(res)
}


## Currently, `rotl::tnrs_match_names` is limited to 10,000 names for exact matches
## -> We need to split species list into chunks
NCHUNKS <- (nrow(spnames) %/% 10000) + 1
spnames <- chunk_table(x = spnames, n = NCHUNKS)


## Function to get OTT ids
get_ott_ids <- function(x, fuzzy = FALSE, group = "All life", progr = "text"){

  res <- adply(
    .data = x,
    .margins = 1,
    .fun = function(sp){
      
      rz <- rotl::tnrs_match_names(
        names = sp,
        context_name = group,
        do_approximate_matching = fuzzy)

      rz <- data.table(
          verbatim_name = sp,
          ott_id = rz$ott_id
          # ,
          # species = rz$unique_name,
          # is_synonym = rz$is_synonym,
          # flags = rz$flags,
          # number_matches = rz$number_matches
          )
      return(rz)
    },
    .progress = progr, .parallel = parall
    )

  setDT(res)
  res[, X1 := NULL ]
  return(res)
}
# e.g., get_ott_ids(spnames[[1]]$species)



cat("\nMatching species names to the Open Tree Taxonomy with `rotl` package\n")

species_ott <- llply(.data = spnames, .fun = function(x){
  get_ott_ids(
    x = x$species,
    fuzzy = FALSE, group = TAXGROUP, progr = "none")
  })

species_ott <- rbindlist(species_ott)

## Add OTT IDs to the species table
spnames <- rbindlist(spnames)
spnames <- merge(x = spnames, y = species_ott,
  by.x = "species", by.y = "verbatim_name", all.x = TRUE)

setcolorder(x = spnames, neworder = c(
  "ott_id", "specieskey",
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
  ))


cat("..Number of species found in OTT: ", length(na.omit(unique(spnames$ott_id))), "\n")
cat("..Number of species without matches in OTT: ", sum(is.na(spnames$ott_id)), "\n")

if(any(is.na(spnames$ott_id))){
  cat("WARNING: there are species without OTT IDs in the data! They will be excluded.\n")
  cat("..GBIF specieskeys without matches in OTT: ", paste(spnames[ is.na(ott_id) ]$specieskey, collapse = ", "), "\n")
  spnames <- spnames[ ! is.na(ott_id) ]
}

if(any(duplicated(spnames$ott_id))){
  cat("WARNING: there are non-unique OTT IDs in the data!\n")
  cat("..Number of duplicated OTT IDs: ", sum(duplicated(spnames$ott_id)), "\n")
  cat("..Duplicates: ", paste(spnames$ott_id[ duplicated(spnames$ott_id) ], collapse = ", "), "\n")
}

## Export results
cat("\nExporting results\n")

fwrite(x = spnames, file = OUTPUT, quote = "auto", sep = ",")

cat("..All done!\n")

#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
