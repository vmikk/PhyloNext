#!/usr/bin/env nextflow
/*
Pipeline for 
========================================================================================
    vmikk/biodiverse
========================================================================================
    Github : https://github.com/vmikk/biodiverse-scripts
    Website: TBA
    Slack  : TBA
----------------------------------------------------------------------------------------
*/

// Enable DSL2 syntax
nextflow.enable.dsl = 2


//// Initialize parameters, set default values

// Filtering, stage I - "10_Filter_occurrences.R"
params.input = false
params.output = "$baseDir/results"
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

// Filtering, stage II - "11_Additional_filtering_and_aggregation.R"
params.specieskey = "NA"
params.h3resolution = 4
params.dbscan = false
params.dbscannoccurrences = 30
params.dbscanepsilon = 700
params.dbscanminpts = 3
params.terrestrial = "NA"

// Filtered data aggregation - "12_Prepare_Biodiverse_input.R"
params.phytree = "NA"
params.taxgroup = "All life"


// Print the parameters to the console
println "-----------------------------------------------------"
println "GBIF occurrence dump:     " + params.input
println "Output path:              " + params.output
println "H3 spatial resolution:    " + params.h3resolution
println "Spatial outliers removal: " + params.dbscan
println "DBSCAN epsilon:           " + params.dbscanepsilon
println "DBSCAN minpts:            " + params.dbscanminpts

println "-----------------------------------------------------"


