#!/usr/bin/env nextflow
/*

========================================================================================
    GBIF phylogenetic diversity pipeline
========================================================================================
    Version: v0.1
    License: MIT
    Github : https://github.com/vmikk/biodiverse-scripts
    Website: TBA
    Slack  : TBA
----------------------------------------------------------------------------------------
*/

// TO DO:
// - specify default path to pipeline data and scripts
// - specify threads for each process in a config: https://www.nextflow.io/docs/latest/process.html#cpus
// - specify Dockerimages: https://www.nextflow.io/docs/latest/process.html#container
// - split the pipeline into workflows?
// - fix publishDir: https://www.nextflow.io/docs/latest/process.html#publishdir


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.1'

//// Initialize parameters, set default values

// Filtering, stage I - "10_Filter_occurrences.R"
params.input = false
params.outdir = "$baseDir/results"
params.phylum = "NA"
params.class = "NA"
params.order = "NA"
params.family = "NA"
params.country = "NA"
params.latmin = "NA"
params.latmax = "NA"
params.lonmin = "NA"
params.lonmax = "NA"
params.minyear = 1945
params.noextinct = "NA"
params.roundcoords = true

params.terrestrial = "pipeline_data/Land_Buffered_025_dgr.RData"


// Filtering, stage II - "11_Additional_filtering_and_aggregation.R"
params.h3resolution = 4
params.dbscan = false
params.dbscannoccurrences = 30
params.dbscanepsilon = 700
params.dbscanminpts = 3
params.terrestrial = "NA"

// Filtered data aggregation - "12_Prepare_Biodiverse_input.R"
params.phytree = "NA"
params.taxgroup = "All life"

// Help message flag
params.helpMsg = false



// Pipeline help message
def helpMsg() {
    log.info"""
    ====================================================================
    GBIF phylogenetic diversity pipeline :  Version ${version}
    ====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run main.nf --input .... --outdir ....
    
    Options:
    REQUIRED:
        --input               Path to the directory with parquet files (GBIF occurrcence dump)
        --outdir              The output directory where the results will be saved
    OPTIONAL:
        --phylum              phylum ...
        --class               class ...
        --order               order ...
        --family              family ...
        --country             country ...
        --latmin              latmin ...
        --latmax              latmax ...
        --lonmin              lonmin ...
        --lonmax              lonmax ...
        --minyear             minyear ...
        --noextinct           noextinct ...
        --roundcoords         roundcoords ...
        --h3resolution        h3resolution ...
        --dbscan              dbscan ...
        --dbscannoccurrences  dbscannoccurrences ...
        --dbscanepsilon       dbscanepsilon ...
        --dbscanminpts        dbscanminpts ...
        --terrestrial         terrestrial ...
    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit(0)
}

// Check if input path was provided
if (params.input == false) {
    println( "Please provide the directory with input data wuth `--input`")
    exit(1)
}


// Print the parameters to the console and to the log
log.info """
        GBIF phylogenetic diversity pipeline
        ===========================================
        GBIF occurrence dump:     ${params.input}
        Output path:              ${params.outdir}
        H3 spatial resolution:    ${params.h3resolution}
        Spatial outliers removal: ${params.dbscan}
        """
        .stripIndent()

if(params.dbscan == true){
    log.info "DBSCAN epsilon:           ${params.dbscanepsilon}".stripIndent()
    log.info "DBSCAN minptsl:           ${params.dbscanminpts}".stripIndent()
} 

log.info "\n"


