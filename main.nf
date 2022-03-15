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
//      change copy to move??


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.1'

//// Initialize parameters, set default values

params.scripts_path = "${baseDir}/bin"
params.data_path = "${baseDir}/pipeline_data"

// Filtering, stage I - "10_Filter_occurrences.R"
params.input = false
params.outdir = "${baseDir}/results"
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
params.dbscannoccurrences = 30

// Filtering, stage II - "11_Additional_filtering_and_aggregation.R"
params.h3resolution = 4
params.dbscan = false
params.dbscanepsilon = 700
params.dbscanminpts = 3
params.terrestrial = params.data_path + "/Land_Buffered_025_dgr.RData"


// Filtered data aggregation - "12_Prepare_Biodiverse_input.R"
params.phytree = "NA"
params.taxgroup = "All_life"

// Help message flag
params.helpMsg = false


// Number of CPUs to use at different steps --> configure as ${task.cpus}
// params.cpusfilt1 = 10
// params.cpusfilt2l = 5
// params.cpusfilt2h = 1
// params.cpusbioprep = 10


// Define output paths for different steps
out_flt1 = params.outdir + "/00.filtered1.parquet"
out_flt2 = params.outdir + "/01.filtered2"
out_biod = params.outdir + "/02.Biodiverse_input"
out_logs = params.outdir + "/logs"


// Pipeline help message
def helpMsg() {
    log.info"""
    ====================================================================
    GBIF phylogenetic diversity pipeline :  Version ${version}
    ====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run main.nf --input ... --outdir ...
    
    Options:
    REQUIRED:
        --input               Path to the directory with parquet files (GBIF occurrcence dump)
        --outdir              The output directory where the results will be saved
    OPTIONAL:
        --phylum              Phylum to analyze (multiple comma-separated values allowed); e.g., "Chordata"
        --class               Class to analyze (multiple comma-separated values allowed); e.g., "Mammalia"
        --order               Order to analyze (multiple comma-separated values allowed); e.g., "Carnivora"
        --family              Family to analyze (multiple comma-separated values allowed); e.g., "Felidae,Canidae"
        --country             Country code, ISO 3166 (multiple comma-separated values allowed); e.g., "DE,PL,CZ"
        --latmin              Minimum latitude of species occurrences (decimal degrees); e.g., 5.1
        --latmax              Maximum latitude of species occurrences (decimal degrees); e.g., 15.5
        --lonmin              Minimum longitude of species occurrences (decimal degrees); e.g., 47.0
        --lonmax              Maximum longitude of species occurrences (decimal degrees); e.g., 55.5
        --minyear             Minimum year of record's occurrences; default, 1945
        --noextinct           File with extinct species specieskeys for their removal
        --roundcoords         Logical, round spatial coordinates to two decimal places, to reduce the dataset size (default, TRUE)
        --h3resolution        Spatial resolution of the H3 geospatial indexing system; e.g., 4
        --dbscan              Logical, remove spatial outliers with density-based clustering; e.g., "false"
        --dbscannoccurrences  Minimum species occurrence to perform DBSCAN; e.g., 30
        --dbscanepsilon       DBSCAN parameter epsilon, km; e.g., "700"
        --dbscanminpts        DBSCAN min number of points; e.g., "3"
        --terrestrial         Land polygon for removal of non-terrestrial occurrences; e.g., "pipeline_data/Land_Buffered_025_dgr.RData"
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
        Land-mask:                ${params.terrestrial}
        Spatial outliers removal: ${params.dbscan}
        """
        .stripIndent()

if(params.dbscan == true){
    log.info "DBSCAN epsilon:           ${params.dbscanepsilon}".stripIndent()
    log.info "DBSCAN minptsl:           ${params.dbscanminpts}".stripIndent()
} 

log.info "==========================================="
log.info "\n"


// Input directory with parquet files (GBIF dump dir)
input_ch = Channel.value(params.input)


// Occurrence filtering, stage I
process occ_filter {

    publishDir "$params.outdir", mode: 'copy'
    cpus 10

    input:
      val input

    output:
      val "${out_flt1}/Partition=low", emit: part_low
      val "${out_flt1}/Partition=high", emit: part_high
      path "spp.txt", emit: spp

    script:
    """
    Rscript ${params.scripts_path}/10_Filter_occurrences.R \
      --input ${input} \
      --phylum ${params.phylum} \
      --class ${params.class} \
      --order ${params.order} \
      --family ${params.family} \
      --country ${params.country} \
      --latmin ${params.latmin} \
      --latmax ${params.latmax} \
      --lonmin ${params.lonmin} \
      --lonmax ${params.lonmax} \
      --minyear ${params.minyear} \
      --noextinct ${params.noextinct} \
      --roundcoords ${params.roundcoords} \
      --threads ${task.cpus} \
      --noccurrences ${params.dbscannoccurrences} \
      --output ${out_flt1}

    ## Prepare species list for DBSCAN
    awk '\$3 ~ /high/ {print \$1 }' SpeciesCounts.txt > spp.txt

    """
}


// Outlier filtering, stage II - without DBSCAN, all low abundant species
process outl_low {

    publishDir "$params.outdir/01.filtered2", mode: 'copy'
    cpus 5

    input:
      path(part_low)

    output:
      path "NoSpKey.RData", emit: lowabsp

    script:
    """
    Rscript ${params.scripts_path}/11_Additional_filtering_and_aggregation.R \
      --input "${part_low}" \
      --dbscan false \
      --resolution ${params.h3resolution} \
      --terrestrial ${params.terrestrial} \
      --threads ${task.cpus} \
      --output ${out_flt2}

    cp ${out_flt2}/NoSpKey.RData NoSpKey.RData
    """
}

//  The default workflow
workflow {

    // Run stage-I filtering
    occ_filter(input_ch)
}


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}