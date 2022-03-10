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

nextflow.enable.dsl = 2


// Initialize parameters, set a default values
params.input = false
params.output = "results"
params.h3resolution = 4
params.dbscan = false
params.dbscanepsilon = 700
params.dbscanminpts = 3


// Print the parameters to the console
println "-----------------------------------------------------"
println "GBIF occurrence dump:     " + params.input
println "Output path:              " + params.output
println "H3 spatial resolution:    " + params.h3resolution
println "Spatial outliers removal: " + params.dbscan
println "DBSCAN epsilon:           " + params.dbscanepsilon
println "DBSCAN minpts:            " + params.dbscanminpts

println "-----------------------------------------------------"


