#!/usr/bin/env Rscript

## Aggregate Biodiverse randomization results (from different threads)
## and estimate Z-scores


## Usage example:
# ./13_Aggregate_Biodiverse_results.R \
#   --observed "obs.csv" \
#   --randomized "rnd.txt" \
#   --output "Biodiverse_input" \
#   --threads 1 \
############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-r", "--observed"), action="store", default=NA, type='character', help="Input file (CSV) with Biodiverse results - observed indices"),
  make_option(c("-z", "--randomized"),  action="store", default=NA, type='character', help="Input file (TXT) with paths to Biodiverse randomized results `rand--SPATIAL_RESULTS.csv`"),

  make_option(c("-t", "--threads"), action="store", default=2L, type='integer', help="Number of CPU threads for arrow, default 4"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$observed)){
  stop("Input is not specified.\n")
}
if(is.na(opt$randomized)){
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
OBSERVED <- opt$observed
RANDOMIZED <- opt$randomized
CPUTHREADS <- as.numeric(opt$threads)
OUTPUT <- opt$output

## Log assigned variables
cat(paste("Input, observed values: ", OBSERVED, "\n", sep=""))
cat(paste("Input, randomized values: ", RANDOMIZED, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))
cat(paste("Output file: ", OUTPUT, "\n", sep=""))

cat("\n")


############################################## Load packages

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("h3")
load_pckg("plyr")

cat("\n")


## Set CPU thread pool
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table

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


############################################## Main pipeline

cat("Loading input data...\n")

## Load observed values
obs <- fread(OBSERVED)
obs[, Axis_0 := NULL ]

## Load randomization results  - "*_rand--SPATIAL_RESULTS.csv"
rndfiles <- read.delim(file = RANDOMIZED, header = FALSE, sep = "\t", col.names = "FileID")
datt <- alply(.data = rndfiles$FileID, .margins = 1, .fun = fread)

## Bind randomization chunks together
datt <- rbindlist(datt)

## Remove useless columns
datt[, Axis_0 := NULL ]


cat("Processing data...\n")

## Reshape data to long format
datl <- melt(data = datt, id.vars = c("ELEMENT"),
  variable.name = "IndexName", value.name = "Value")

## Add index types
datl[, IndexType := tstrsplit(x = IndexName, split = "_", keep = 1)]

## Estimate the number of randomization iterations
## beacause of undefined values, it could be different for different samples and indices
niter <- datl[ 
  IndexType %in% "Q" , 
  .(NIter = sum(Value, na.rm = TRUE)), 
  by = c("ELEMENT", "IndexName")]

niter[, IndexName := gsub(pattern = "^Q_", replacement = "", x = IndexName) ]


## Prepare data for Z-score estimation
ZSCOREL <- datl[ 
  IndexType %in% c("SUMX", "SUMXX"),
  .(Sum = sum(Value, na.rm = TRUE)),
  , by = c("ELEMENT", "IndexName", "IndexType")]

## Remove index name prefix
ZSCOREL[, IndexName := gsub(pattern = "^SUMX_|^SUMXX_", replacement = "", x = IndexName) ]

## Reshape to wide format
ZSCORE <- dcast(data = ZSCOREL,
  formula = ELEMENT + IndexName ~ IndexType,
  fun.aggregate = sum, fill = NA, value.var = "Sum")


## Reshape observed values into long format
obsl <- melt(data = obs, id.vars = "ELEMENT",
  variable.name = "IndexName", value.name = "Observed")

## Add observed values to the table with randomizations
if(any(!unique(ZSCORE$IndexName) %in% unique(obsl$IndexName))){
  cat("WARNING: index names mismatches with observed data.\n")
}

ZSCORE <- merge(x = ZSCORE, y = obsl, by = c("ELEMENT", "IndexName"), all.x = TRUE)


## There are some indices without SUMX and SUMXX
# ii <- unique(obsl$IndexName)[ !unique(obsl$IndexName) %in% unique(ZSCORE$IndexName) ]
# obsl[ IndexName %in% ii ]
##     e.g., PMPD1_MAX  PMPD1_MEAN PMPD1_MIN  PMPD1_RMSD PNTD1_MAX  PNTD1_MEAN PNTD1_MIN  PNTD1_RMSD


## Add number of iterations
if(any(!unique(ZSCORE$IndexName) %in% unique(niter$IndexName))){
  cat("WARNING: index names mismatches with randomization-based data.\n")
}

ZSCORE <- merge(x = ZSCORE, y = niter, by = c("ELEMENT", "IndexName"), all.x = TRUE)


cat("Estimating Z-scores...\n")

## Estimate variance of null-distribuition (from randomized values)
ZSCORE[ , VAR := ((SUMXX/NIter) - (SUMX/NIter)^2)  ]

## Define tolerance threshold (a value less than tol is considered equal to zero)
## Biodiverse uses tol = 1e-10  (https://github.com/shawnlaffan/biodiverse/commit/5a1ca6aaa11c95ab870e8553f4d7750c37af9efd)
tol = .Machine$double.eps              # 2.220446e-16
ZSCORE[ abs(VAR) < tol  , VAR := 0 ]

## SD of null-distribuition
ZSCORE[ , SD := sqrt(VAR) ]

## Estimate Z-scores
ZSCORE[ , ZScore := (Observed - (SUMX/NIter)) / SD ]

## Replace -Inf and Inf by zero
ZSCORE[ is.infinite(ZScore), ZScore := 0 ]

## Replace NaN by zero
ZSCORE[ is.nan(ZScore), ZScore := 0 ]


## How to deal with NAs??


## Add coordinates to the table
coords <- data.table(ELEMENT = unique(obs$ELEMENT))
crd <- h3::h3_to_geo(h3_index = coords$ELEMENT)
coords <- cbind(coords, crd)
rm(crd)

ZSCORE <- merge(x = ZSCORE, y = coords, by = "ELEMENT", all.x = TRUE)

## Reorder columns
setcolorder(ZSCORE, c("ELEMENT", "lat", "lng", "IndexName"))


## Export results
cat("Exporting data...\n")

fwrite(x = ZSCORE, file = OUTPUT, quote = F, sep = "\t")



## Compare with Biodiverse Z-scores
# zzl <- melt(data = zz, id.vars = "ELEMENT", variable.name = "IndexName", value.name = "BZ")
# ZSCORE <- merge(x = ZSCORE, y = zzl, by = c("ELEMENT", "IndexName"), all.x = TRUE)
# ggplot(data = ZSCORE, aes(x = BZ, y = ZScore)) + 
#   geom_point() + annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf) + facet_wrap(~ IndexName, scales = "free") +
#   labs(x = "Biodiverse Z-score", y = "Manual Z-score")



cat("All done.\n")

#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
