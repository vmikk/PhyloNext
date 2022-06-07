#!/usr/bin/env Rscript

## Prepare input data for Biodiverse

## TO DO:
# - add option to retrive tree via Open Tree API
# - export GBIF-OTT ID matches
# - export species IDs that are not in tree
# - check tree nodes (in case if tree is not strictly species-level)
# - handle the case without phylogentic tree (do we need it?)
# - use pre-computed OTT-vs-GIBF IDs table (instead of OTT API)

## Name matching with Open Tree of life:
# - currently fuzzy-matching is disabled (only exact matches are used)

## Usage example:
# ./12_Prepare_Biodiverse_input.R \
#   --input "filtered_data" \
#   --phytree "phy_trees/Fabaceae_Australian.nwk" \
#   --taxgroup "Rosids" \
#   --threads 12 \
#   --output "Biodiverse_input"


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the directory with filtered data"),
  make_option(c("-f", "--inputfile"), action="store", default=NA, type='character', help="Text file with paths to the filtered data"),
  make_option(c("-p", "--phytree"), action="store", default=NA, type='character', help="Phylogenetic tree (pre-computed in Newick format or get it from API)"),
  make_option(c("-g", "--taxgroup"), action="store", default="All life", type='character', help="Taxonomy group in OpenTree"),

  make_option(c("-t", "--threads"), action="store", default=2L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(sum(is.na(opt$input), is.na(opt$inputfile)) < 1){
  stop("Input is not specified.\n")
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
INPUT <- opt$input
INPUTFILE <- opt$inputfile
PHYTREE <- to_na( opt$phytree )
TAXGROUP <- opt$taxgroup
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT <- opt$output

## Modify tax group argument (no spaces are allowed in arguments)
# if(TAXGROUP %in% "All_life"){ TAXGROUP <- "All life" }
TAXGROUP <- gsub(pattern = "_", replacement = " ", x = TAXGROUP)


## Log assigned variables
cat(paste("Input directory with filtered occurrences: ", INPUT, "\n", sep=""))
cat(paste("Input file with paths to the filtered occurrences: ", INPUTFILE, "\n", sep=""))
cat(paste("Pre-computed phylogenetic tree: ", PHYTREE, "\n", sep=""))
cat(paste("Taxonomy group in OpenTree: ", TAXGROUP, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
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

load_pckg("data.table")
load_pckg("plyr")
load_pckg("rotl")     # Interface to the 'Open Tree of Life' API
# load_pckg("rgbif")
load_pckg("ape")

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

## Load phylogenetic tree
if(!is.na(PHYTREE)){
  
  ## Function to extract file extension
  file_extension <- function(x){
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }

  if(file_extension(PHYTREE) %in% c("nwk", "newick", "Newick", "tre", "tree")){

    ## Load Newick tree
    cat("Loading pre-computed phylogenetic tree\n")
    TREE <- ape::read.tree(PHYTREE)
  
  } else if(PHYTREE %in% "API") {
    
    cat("Fetching phylogenetic tree from Open Tree API\n")
    # ////////////////////////////////////////  -- TO DO

  } else {

    stop("Unknown tree option is specified. Please check `--phytree` argument.\n")

  }
}



## Load and combine filtered data
cat("Loading filtered occurrences\n")

## Directory with *.RData files
if(!is.na(INPUT)){
  fls <- list.files(
    path = INPUT,
    pattern = ".RData$",
    full.names = T, recursive = T)
}

## Text file with paths to *.RData files
if(!is.na(INPUTFILE)){
  fls <- read.table(
    file = INPUTFILE,
    header = F, sep = "\t", col.names = "FilePath")$FilePath
}

cat("..Files found: ", length(fls), "\n")
if(length(fls) == 0){
  stop("ERROR: no filtered data found in the speciefied path.\n")
}

datt <- alply(.data = fls, .margins = 1,
  .fun = function(z){ readRDS(z) },
  .progress = "none", .parallel = parall)

## Count number of gridcells-outliers
cat("..Counting the number of outliers\n")

outliers_terrestrial <- laply(.data = datt, .fun = function(z){
  length( na.omit(attr(z, "removed_nonterrestrial_h3")) )
  })

outliers_dbscan <- laply(.data = datt, .fun = function(z){
  length( na.omit(attr(z, "removed_dbscan_h3")) )
  })

cat("...The number of excluded non-terrestrial gridcells: ", sum(outliers_terrestrial), "\n")
cat("...The number of excluded DBSCAN-based gridcells-outliers: ", sum(outliers_dbscan), "\n")


cat("..Merging filtered occurrences from different specieskeys\n")
datt <- rbindlist(datt)

## Check data
if(!nrow(datt) > 1){
  stop("No species occurrences found. Please check the `--input` argument.\n")
}

## Data summary
cat("..Total number of records: ", nrow(datt), "\n")
cat("..Number of unique specieskeys: ", length(unique(datt$specieskey)), "\n")
cat("..Number of unique gridcells: ", length(unique(datt$H3)), "\n")



######################################
###################################### Map GBIF IDs to Open Tree of Life (OTT) IDs
######################################

## Unique species IDs
species_uniq <- unique(datt[, .(specieskey, species)])

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
          id = rz$ott_id,
          species = rz$unique_name,
          is_synonym = rz$is_synonym,
          flags = rz$flags,
          number_matches = rz$number_matches
          )
      return(rz)
    },
    .progress = progr, .parallel = parall
    )

  setDT(res)
  res[, X1 := NULL ]
  colnames(res)[-1] <- paste0("ott_", colnames(res)[-1])
  return(res)
}


## To list available OTT tax groups run
# rotl::tnrs_contexts()
#
## as of rotl v.3.0.12 these contexts are available
# Possible contexts:
#    Animals 
#       Birds, Tetrapods, Mammals, Amphibians, Vertebrates 
#       Arthropods, Molluscs, Nematodes, Platyhelminthes, Annelids 
#       Cnidarians, Arachnids, Insects 
#    Fungi 
#       Basidiomycetes, Ascomycetes 
#    All life 
#    Bacteria 
#       SAR group, Archaea, Excavata, Amoebozoa, Centrohelida 
#       Haptophyta, Apusozoa, Diatoms, Ciliates, Forams 
#    Land plants 
#       Hornworts, Mosses, Liverworts, Vascular plants, Club mosses 
#       Ferns, Seed plants, Flowering plants, Monocots, Eudicots 
#       Rosids, Asterids, Asterales, Asteraceae, Aster 
#       Symphyotrichum, Campanulaceae, Lobelia 

cat("Matching species names to the Open Tree Taxonomy\n")

species_ott <- get_ott_ids(x = species_uniq$species,
  fuzzy = FALSE, group = TAXGROUP, progr = "none")

cat("..Species names with OTT IDs found", sum(!is.na(species_ott$ott_id)), "\n")
cat("..Species names without matches to Open Tree", sum(is.na(species_ott$ott_id)), "\n")


# species_exact <- species_ott[ approximate_match == FALSE ]
# cat("..Number of exact matches", nrow(species_exact), "\n")

# species_fuzzy <- species_ott[ approximate_match != FALSE ]
# cat("..Number of fuzzy matches", nrow(species_fuzzy), "\n")


## Add OTT IDs to the species table
species_uniq <- merge(
  x = species_uniq,
  y = species_ott[ , .(verbatim_name, ott_id)],
  by.x = "species", by.y = "verbatim_name", all.x = TRUE)

## Remove species without OTT matches
species_uniq <- species_uniq[ !is.na(ott_id) ]

## Check for duplicates in OTT IDs
if(any(duplicated(species_uniq$ott_id))){
  cat("Warning: There are OTT IDs matching multiple GBIF specieskeys. These records will be merged\n")

  ## Inspect the records
  # species_uniq[ OTT %in% species_uniq$OTT[ duplicated(species_uniq$OTT) ] ]
}


## Check for duplicates in species keys
# species_uniq[ specieskey %in% species_uniq$specieskey[ duplicated(species_uniq$specieskey) ] ]

## Remove duplicates
species_uniq <- unique(species_uniq, by = c("ott_id", "specieskey") )

# cat("Exporting OTT IDs\n")
# saveRDS(object = species, file = "species_OTT.RData", compress = "xz")



######################################
###################################### Trim the phylogenetic tree
######################################

cat("Preparing phylogenetic tree\n")

## In the tree, tip names have "ott" prefix (e.g., "ott114")
## Add this prefix to OTT IDs
species_uniq[, OTT := paste0("ott", ott_id)]

## Summary
uniq_otts <- unique(species_uniq$OTT)                 # there could be duplicates
in_tree_total <- sum(uniq_otts %in% TREE$tip.label)
in_tree_percent <- round(in_tree_total / length(uniq_otts) * 100, 1)

cat("..Phylogenetic tree contains ", length(TREE$tip.label), "tips\n")
cat("..In total, there are ", length(unique(species_uniq$specieskey)), " unique specieskeys (with OTT matches) in GBIF data.\n")
cat("..Of them, ", in_tree_total, "(", in_tree_percent, "%) are present at the tips of the tree.\n")


## These OTT IDs should be preserved (in occurrence data and on a tree)
otts_to_keep <- uniq_otts[ uniq_otts %in% TREE$tip.label ]

## Subset tree
## NB. only tree tips are considered now, but there could be named nodes!
cat("..Timming phylogentic tree:\n")
cat("..Number of tips to remove from the tree: ", length(TREE$tip.label) - length(otts_to_keep), "\n")
tree_trimmed <- ape::keep.tip(phy = TREE, tip = otts_to_keep)


## Export tree in Nexus format
cat("..Exporting trimmed phylogentic tree for Biodiverse\n")

ape::write.nexus(tree_trimmed,
  file = file.path(OUTPUT, "Trimmed_tree.nex"))


######################################
###################################### Subset species
######################################

cat("Preparing occurrence data\n")

records_before_filtering <- nrow(datt)

## Add OTT IDs to occurrences
cat("..Adding OTT IDs to occurrences\n")

if(any(duplicated(species_uniq$OTT))){
  cat("WARNING: duplicated")
}

datt <- merge(
  x = datt,
  y = species_uniq[, .(specieskey, OTT)],
  by = "specieskey", all.x = TRUE)


## Remove taxa
cat("..Removing taxa without OTT IDs\n")

# datt <- datt[ !is.na(OTT) ]
datt <- datt[ OTT %in% otts_to_keep  ]

records_after_filtering <- nrow(datt)
records_precent <- round(records_after_filtering / records_before_filtering * 100, 1)

cat("..After removal of not-in-tree taxa, dataset is comprised of ",
  records_after_filtering, "(", records_precent, "%) entries.\n")


## Remove duplicated IDs
if(any(duplicated(species_uniq$ott_id))){
  cat("..Merging duplicated OTT IDs\n")
  datt <- unique(datt, by = c("OTT", "H3"))
  cat("..Number of records without duplicates: ", nrow(datt), "\n")
}

# any(is.na(datt$OTT))
# any(datt$OTT %in% "")


## Export data (in long format) for Biodiverse
cat("..Exporting trimmed occurrence data for Biodiverse\n")
fwrite(
  x = datt[, .(OTT, H3, specieskey, species, decimallatitude, decimallongitude)],
  file = file.path(OUTPUT, "Trimmed_occurrences.csv"),
  sep = ",", quote = T, row.names = FALSE, col.names = TRUE)



## Coordinates of the observed gridcells
cat("Exporting gridcell coordinates\n")
uniq_h3 <- unique(datt[, .(H3, decimallatitude, decimallongitude)])

cat("..Number of unique gridcells after trimming: ", nrow(uniq_h3), "\n")

fwrite(
  x = uniq_h3,
  file = file.path(OUTPUT, "H3_GridCell_Centres.csv"),
  sep = "\t", quote = F, row.names = FALSE, col.names = TRUE)



#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
