#!/usr/bin/env nextflow
/*

========================================================================================
    PhyloNext: GBIF phylogenetic diversity pipeline
========================================================================================
    License: MIT
    Github : https://github.com/vmikk/phylonext
    Website: https://phylonext.github.io/
----------------------------------------------------------------------------------------
*/


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Parameters are defined in the `conf/params.config` file

// // Biodiverse args
// biodiverse_args = "function=" + params.randname + " max_iters=" + params.iterations

// How many randomization iterations should be per Biodiverse thread?
if( params.biodiversethreads > 1 ) {
    randomization_chunks = (0..<params.biodiversethreads)
    iterations_per_thread = params.iterations / params.biodiversethreads
    iterations_per_thread = Math.ceil(iterations_per_thread).toInteger()
}
else {
    randomization_chunks = [ 1 ]
    iterations_per_thread = params.iterations
}
biodiverse_args = "function=" + params.randname + " max_iters=" + iterations_per_thread


// OTT-derived tree should have Latin tip names
if(params.phytree == null & params.phylabels == "OTT"){
  println("No user-supplied phylogenetic tree was provided.")
  println("Tree will be automatically fetched from the Open Tree of Life.")
  println("!! Please set `--phylabels` parameter to `Latin`")
  exit(1)
}

// Optional input files
specieskeys        = params.specieskeys        ? file(params.specieskeys)        : file("${params.outdir}/no_file0")
noextinct          = params.noextinct          ? file(params.noextinct)          : file("${params.outdir}/no_file1")
terrestrial        = params.terrestrial        ? file(params.terrestrial)        : file("${params.outdir}/no_file2")
rmcountrycentroids = params.rmcountrycentroids ? file(params.rmcountrycentroids) : file("${params.outdir}/no_file3")
rmcountrycapitals  = params.rmcountrycapitals  ? file(params.rmcountrycapitals)  : file("${params.outdir}/no_file4")
rminstitutions     = params.rminstitutions     ? file(params.rminstitutions)     : file("${params.outdir}/no_file5")
rmurban            = params.rmurban            ? file(params.rmurban)            : file("${params.outdir}/no_file6")
wgsrpd             = params.wgsrpd             ? file(params.wgsrpd)             : file("${params.outdir}/no_file7")
phytree            = params.phytree            ? file(params.phytree)            : file("${params.outdir}/no_file8")
world              = params.world              ? file(params.world)              : file("${params.outdir}/no_file9")
polygon            = params.polygon            ? file(params.polygon)            : file("${params.outdir}/no_file10")


// Pipeline help message
def helpMsg() {
    log.info"""
    =====================================================================
    PhyloNext: GBIF phylogenetic diversity pipeline :  Version ${workflow.manifest.version}
    =====================================================================

    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run vmikk/phylonext -r main --input ... --outdir ...

    Options:
    REQUIRED:
        --input               Path to the directory with parquet files (GBIF occurrcence dump)
        --outdir              The output directory where the results will be saved
    OPTIONAL:
        --phylum              Phylum to analyze (multiple comma-separated values allowed); e.g., "Chordata"
        --classis             Class to analyze (multiple comma-separated values allowed); e.g., "Mammalia"
        --order               Order to analyze (multiple comma-separated values allowed); e.g., "Carnivora"
        --family              Family to analyze (multiple comma-separated values allowed); e.g., "Felidae,Canidae"
        --genus               Genus to analyze (multiple comma-separated values allowed); e.g., "Felis,Canis,Lynx"
        --specieskeys         Custom list of GBIF specieskeys (file with a single column, with header)

        --phytree             Custom phylogenetic tree
        --taxgroup            Specific taxonomy group in Open Tree of Life (default, "All_life")
        --phylabels           Type of tip labels on a phylogenetic tree ("OTT" or "Latin")
        --maxage              Manually assign root age for a tree obtained from Open Tree of Life; e.g., 127
        --phyloonly           Prune Open Tree tips for which there are no phylogenetic inputs; logical, default, false

        --country             Country code, ISO 3166 (multiple comma-separated values allowed); e.g., "DE,PL,CZ"
        --latmin              Minimum latitude of species occurrences (decimal degrees); e.g., 5.1
        --latmax              Maximum latitude of species occurrences (decimal degrees); e.g., 15.5
        --lonmin              Minimum longitude of species occurrences (decimal degrees); e.g., 47.0
        --lonmax              Maximum longitude of species occurrences (decimal degrees); e.g., 55.5
        --minyear             Minimum year of record's occurrences; default, 1945
        --maxyear             Maximum year of record's occurrences; default, none
        --coordprecision      Coordinate precision threshold (less than maximum allowed value; default, 0.1)
        --coorduncertainty    Maximum allowed coordinate uncertainty, meters (default, 10000)
        --coorduncertaintyexclude Black list of coordinate uncertainty values (default, "301,3036,999,9999")
        --basisofrecordinclude Basis of record to include from the data; e.g., "PRESERVED_SPECIMEN"
        --basisofrecordexclude Basis of record to exclude from the data; e.g., "FOSSIL_SPECIMEN,LIVING_SPECIMEN"
        --polygon             Custom area of interest (a file with polygons in GeoPackage format)
        --wgsrpd              Polygons of World Geographical Regions; e.g., "pipeline_data/WGSRPD.RData"
        --regions             Names of World Geographical Regions; e.g., "L1_EUROPE,L1_ASIA_TEMPERATE"
        --noextinct           File with extinct species specieskeys for their removal (file with a single column, with header)
        --excludehuman        Logical, exclude genus "Homo" from occurrence data (default, true)
        --roundcoords         Numeric, round spatial coordinates to N decimal places, to reduce the dataset size (default, 2; set to negative to disable rounding)
        --h3resolution        Spatial resolution of the H3 geospatial indexing system; e.g., 4

        --dbscan              Logical, remove spatial outliers with density-based clustering; e.g., "false"
        --dbscannoccurrences  Minimum species occurrence to perform DBSCAN; e.g., 30
        --dbscanepsilon       DBSCAN parameter epsilon, km; e.g., "700"
        --dbscanminpts        DBSCAN min number of points; e.g., "3"

        --terrestrial         Land polygon for removal of non-terrestrial occurrences; e.g., "pipeline_data/Land_Buffered_025_dgr.RData"
        --rmcountrycentroids  Polygons with country and province centroids; e.g., "pipeline_data/CC_CountryCentroids_buf_1000m.RData"
        --rmcountrycapitals   Polygons with country capitals; e.g., "pipeline_data/CC_Capitals_buf_10000m.RData"
        --rminstitutions      Polygons with biological institutuions and museums; e.g., "pipeline_data/CC_Institutions_buf_100m.RData"
        --rmurban             Polygons with urban areas; e.g., "pipeline_data/CC_Urban.RData"

        --deriveddataset      Prepare a list of DOIs for the datasets used (default, true)

        --indices             Comma-seprated list of diversity and endemism indices; e.g., "calc_richness,calc_pd,calc_pe"
        --randname            Randomisation scheme type; e.g., "rand_structured"
        --iterations          Number of randomisation iterations; e.g., 1000
        --biodiversethreads   Number of Biodiverse threads; e.g., 10
        --randconstrain       Polygons to perform spatially constrained randomization (GeoPackage format)

    Leaflet interactive visualization:
        --leaflet_var         Variables to plot; e.g., "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE,PE_WE,SES_PE_WE,CANAPE,Redundancy"
        --leaflet_canapesuper Include the `superendemism` class in CANAPE results (default, false)
        --leaflet_color       Color scheme for continuous variables (default, "RdYlBu")
        --leaflet_palette     Color palette for continuous variables (default, "quantile")
        --leaflet_bins        Number of color bins for continuous variables (default, 5)
        --leaflet_sescolor    Color scheme for standardized effect sizes, SES (default, "threat"; alternative - "hotspots)
        --leaflet_redundancy  Redundancy threshold for hiding the grid cells with low number of records (default, 0 = display all grid cells)
        --leaflet_zoom        Zoom factor of the image for Leaflet-based screenshots (default, 2)
        --leaflet_width       Image width of Leaflet-based screenshots (default, 3840)
        --leaflet_height      Image height of Leaflet-based screenshots (default, 2160)

    Static visualization:
        --plotvar             Variables to plot (multiple comma-separated values allowed); e.g., "RICHNESS_ALL,PD,PD_P"
        --plottype            Plot type
        --plotformat          Plot format (jpg,pdf,png)
        --plotwidth           Plot width (default, 18 inches)
        --plotheight          Plot height (default, 18 inches)
        --plotunits           Plot size units (in,cm)
        --world               World basemap

    NEXTFLOW-SPECIFIC:
        -qs                   Queue size (max number of processes that can be executed in parallel); e.g., 8
        -w                    Path to the working directory to store intermediate results (default, "./work")
        -resume               Execute the pipeline using the cached results.<br>Useful to continue executions that was stopped by an error
        -profile              Configuration profile; e.g., "docker"
        -params-file          Parameter file in YAML or JSON format (e.g., "Mammals.yaml")
        -c / -C               Configuration file (`-C` ignores all default values) (default, "nextflow.config")
    """.stripIndent()
}
// Show help msg
if (params.helpMsg | params.help){
    helpMsg()
    exit(0)
}

// Check if input path was provided
if (params.input == null) {
    println( "Please provide the directory with input data with `--input`")
    exit(1)
}


// Print the parameters to the console and to the log
log.info """
        ====================================================================
        PhyloNext: GBIF phylogenetic diversity pipeline, Version ${workflow.manifest.version}
        ====================================================================
        GBIF occurrence dump:     ${params.input}
        Output path:              ${params.outdir}
        Launch directory:         ${workflow.launchDir}
        Working directory:        ${workflow.workDir}
        """
        .stripIndent()

log.info "===================================================================="
log.info "\n"


// Occurrence filtering, stage I
process occ_filter {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/00.Filtered1.parquet", mode: 'copy'
    // cpus 10

    input:
      path(input)
      path(noextinct)
      path(polygon)
      path(wgsrpd)

    output:
      path "Partition_low",     emit: part_low,  type: "dir", optional: true
      path "Partition_high",    emit: part_high, type: "dir", optional: true
      path "spp.txt",           emit: spp
      path "SpeciesCounts.txt", emit: spp_all

    script:

    filter_specieskeys = params.specieskeys ? "--specieskeys $specieskeys" : ""
    filter_extinct     = params.noextinct   ? "--noextinct $noextinct"     : ""
    filter_polygon     = params.polygon     ? "--polygon $polygon"         : ""
    filter_wgsrpd      = params.wgsrpd      ? "--wgsrpd $wgsrpd --regions ${params.regions}" : ""

    """
    10_Filter_occurrences.R \
      --input   ${input} \
      --phylum  ${params.phylum} \
      --class   ${params.classis} \
      --order   ${params.order} \
      --family  ${params.family} \
      --genus   ${params.genus} \
      --country ${params.country} \
      --latmin  ${params.latmin} \
      --latmax  ${params.latmax} \
      --lonmin  ${params.lonmin} \
      --lonmax  ${params.lonmax} \
      --minyear ${params.minyear} \
      --maxyear ${params.maxyear} \
      --coordprecision          ${params.coordprecision} \
      --coorduncertainty        ${params.coorduncertainty} \
      --coorduncertaintyexclude ${params.coorduncertaintyexclude} \
      --basisofrecordinclude ${params.basisofrecordinclude} \
      --basisofrecordexclude ${params.basisofrecordexclude} \
      ${filter_specieskeys} \
      ${filter_extinct} \
      ${filter_polygon} \
      ${filter_wgsrpd} \
      --excludehuman ${params.excludehuman} \
      --roundcoords  ${params.roundcoords} \
      --threads      ${task.cpus} \
      --noccurrences ${params.dbscannoccurrences} \
      --output "."

    ## Prepare species list for DBSCAN
    awk '\$3 ~ /high/ {print \$1 }' SpeciesCounts.txt > spp.txt

    ## Remove `equal` sign (causes a problem with azcopy)
    mv 'Partition=low'  Partition_low
    mv 'Partition=high' Partition_high

    ## Check the size of output directories, if empty - remove them
    if [ -d 'Partition_low'  ]; then rmdir --ignore-fail-on-non-empty 'Partition_low'  ; fi
    if [ -d 'Partition_high' ]; then rmdir --ignore-fail-on-non-empty 'Partition_high' ; fi

    """
}


// Counting the total number of records per H3 cell ("sampling effort" for the redundancy estimator)
process record_count {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/01.NumRecords", mode: 'copy'
    // cpus 10

    input:
      path(input)
      path(noextinct)
      path(polygon)
      path(wgsrpd)
      path(terrestrial)
      path(rmcountrycentroids)
      path(rmcountrycapitals)
      path(rminstitutions)
      path(rmurban)

    output:
      path "Record_counts_H3.RData",  emit: n_recr
      path "Record_counts_H3.txt.gz", emit: n_rect
      path "Record_counts_H3_PerSpecies.RData",  emit: n_spr
      path "Record_counts_H3_PerSpecies.txt.gz", emit: n_spt
      path "Record_counts_Outliers.txt.gz", emit: n_outl, optional: true

    script:

    filter_specieskeys  = params.specieskeys        ? "--specieskeys $specieskeys" : ""
    filter_extinct      = params.noextinct          ? "--noextinct $noextinct"     : ""
    filter_polygon      = params.polygon            ? "--polygon $polygon"         : ""
    filter_terrestrial  = params.terrestrial        ? "--terrestrial $terrestrial" : ""
    filter_country      = params.rmcountrycentroids ? "--rmcountrycentroids $rmcountrycentroids" : ""
    filter_capitals     = params.rmcountrycapitals  ? "--rmcountrycapitals $rmcountrycapitals"   : ""
    filter_institutions = params.rminstitutions     ? "--rminstitutions $rminstitutions" : ""
    filter_urban        = params.rmurban            ? "--rmurban $rmurban" : ""
    filter_wgsrpd       = params.wgsrpd             ? "--wgsrpd $wgsrpd --regions ${params.regions}" : ""

    """
    10_Record_counts.R \
      --input   ${input} \
      --phylum  ${params.phylum} \
      --class   ${params.classis} \
      --order   ${params.order} \
      --family  ${params.family} \
      --genus   ${params.genus} \
      --country ${params.country} \
      --latmin  ${params.latmin} \
      --latmax  ${params.latmax} \
      --lonmin  ${params.lonmin} \
      --lonmax  ${params.lonmax} \
      --minyear ${params.minyear} \
      --maxyear ${params.maxyear} \
      --coordprecision          ${params.coordprecision} \
      --coorduncertainty        ${params.coorduncertainty} \
      --coorduncertaintyexclude ${params.coorduncertaintyexclude} \
      --basisofrecordinclude ${params.basisofrecordinclude} \
      --basisofrecordexclude ${params.basisofrecordexclude} \
      ${filter_specieskeys} \
      ${filter_extinct} \
      --excludehuman ${params.excludehuman} \
      ${filter_polygon} \
      ${filter_wgsrpd} \
      ${filter_terrestrial} \
      ${filter_country} \
      ${filter_capitals} \
      ${filter_institutions} \
      ${filter_urban} \
      --roundcoords ${params.roundcoords} \
      --resolution  ${params.h3resolution} \
      --threads     ${task.cpus} \
      --rcode       \$(which "Shapefile_filters.R") \
      --output      "Record_counts"

    """
}


// Outlier filtering, stage II - without DBSCAN, all low abundant species
process outl_low {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/01.Filtered2_by_species", mode: 'copy'
    // cpus 5

    input:
      path(part_low)
      path(polygon)
      path(wgsrpd)
      path(terrestrial)
      path(rmcountrycentroids)
      path(rmcountrycapitals)
      path(rminstitutions)
      path(rmurban)

    output:
      path "NoSpKey.RData", emit: lowabsp, optional: true
      path "NoSpKey_OutlierCounts.txt", emit: outlierslow, optional: true

    script:

    filter_polygon      = params.polygon            ? "--polygon $polygon"         : ""
    filter_terrestrial  = params.terrestrial        ? "--terrestrial $terrestrial" : ""
    filter_country      = params.rmcountrycentroids ? "--rmcountrycentroids $rmcountrycentroids" : ""
    filter_capitals     = params.rmcountrycapitals  ? "--rmcountrycapitals $rmcountrycapitals" : ""
    filter_institutions = params.rminstitutions     ? "--rminstitutions $rminstitutions" : ""
    filter_urban        = params.rmurban            ? "--rmurban $rmurban" : ""
    filter_wgsrpd       = params.wgsrpd             ? "--wgsrpd $wgsrpd --regions ${params.regions}" : ""

    """
    11_Additional_filtering_and_aggregation.R \
      --input "${part_low}" \
      --dbscan false \
      --resolution ${params.h3resolution} \
      ${filter_polygon} \
      ${filter_wgsrpd} \
      ${filter_terrestrial} \
      ${filter_country} \
      ${filter_capitals} \
      ${filter_institutions} \
      ${filter_urban} \
      --threads ${task.cpus} \
      --rcode   \$(which "Shapefile_filters.R") \
      --output "."

    """
}


// Outlier filtering, stage II - with DBSCAN, independently by species
process outl_high {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/01.Filtered2_by_species", mode: 'copy'
    // cpus 1

    // Add species ID to the log file
    tag "$sp"

    // Iterate for each species from `species_ch` channel
    input:
      val sp
      path(part_high)
      path(polygon)
      path(wgsrpd)
      path(terrestrial)
      path(rmcountrycentroids)
      path(rmcountrycapitals)
      path(rminstitutions)
      path(rmurban)

    output:
      path "${sp}.RData", emit: sp, optional: true
      path "${sp}_OutlierCounts.txt", emit: outliers, optional: true

    script:

    filter_polygon      = params.polygon            ? "--polygon $polygon"         : ""
    filter_terrestrial  = params.terrestrial        ? "--terrestrial $terrestrial" : ""
    filter_country      = params.rmcountrycentroids ? "--rmcountrycentroids $rmcountrycentroids" : ""
    filter_capitals     = params.rmcountrycapitals  ? "--rmcountrycapitals $rmcountrycapitals" : ""
    filter_institutions = params.rminstitutions     ? "--rminstitutions $rminstitutions" : ""
    filter_urban        = params.rmurban            ? "--rmurban $rmurban" : ""
    filter_wgsrpd       = params.wgsrpd             ? "--wgsrpd $wgsrpd --regions ${params.regions}" : ""

    """
    11_Additional_filtering_and_aggregation.R \
      --input     "${part_high}" \
      --specieskey ${sp} \
      --dbscan     ${params.dbscan} \
      --epsilon    ${params.dbscanepsilon} \
      --minpts     ${params.dbscanminpts} \
      --resolution ${params.h3resolution} \
      ${filter_polygon} \
      ${filter_wgsrpd} \
      ${filter_terrestrial} \
      ${filter_country} \
      ${filter_capitals} \
      ${filter_institutions} \
      ${filter_urban} \
      --threads ${task.cpus} \
      --rcode   \$(which "Shapefile_filters.R") \
      --output "."

    """
}


// // Create a file with paths to all chunks with filtered results
// // If there are a lot of files - it triggers and error:
// // String too long. The given string is 265799 Unicode code units long, but only a maximum of 65535 is allowed.
// process filtered_filelist {
//
//     // container image is required for Cloud only
//     label "container_r"
//
//     input:
//     val spp
//
//     output:
//     path "filtered_results.txt", emit: FLT
//
//     shell:
//     $/
//     echo "${spp}" \
//       | sed -z 's/, /\n/g; s/^\[//; s/\]//' \
//       > filtered_results.txt
//     /$
// }


// Prepare a list of OpenTree IDs (to obtain a phylogenetic tree)
process prep_ott_ids {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/02.OTT_tree", mode: 'copy'
    // cpus 1

    input:
      path spp_all

    output:
      path "Species_OTT.csv", emit: spp_ott

    script:
    """
    11_GBIF_SpeciesKey_to_OTTID.R \
      --input    ${spp_all} \
      --taxgroup ${params.taxgroup} \
      --threads  ${task.cpus} \
      --output   "Species_OTT.csv"
    """
}


// Get an induced subtree from synthetic phylogenetic tree for a set of OpenTree IDs
process get_ott_tree {

    label "container_ott"
    queue "custom_pool"

    publishDir "$params.outdir/02.OTT_tree", mode: 'copy'
    // publishDir "$params.tracedir", pattern: "citations.txt", mode: 'copy', saveAs: { filename -> "OTT_tree_citations.txt" }
    // cpus 1

    input:
      path spp_ott

    output:
      path "ott_label_dated_tree.tre", emit: tree
      path "ottid_dated_tree.tre",     emit: treedated
      path "labelled_tree.tre",        emit: treelabeled
      path "ottlabel.txt",             emit: labels
      path "conflict_annot.tre",       emit: conflicts, optional: true
      path "support_annot.tre",        emit: support,   optional: true
      path "citations.txt"
      path "date_citations.txt"
      path "synth.log"

    script:
    max_age    = params.maxage    ? "--max-age ${params.maxage}" : ""
    phylo_only = params.phyloonly ? "--phylo-only" : ""

    """
    induced_synth_subtree_from_csv.py \
      --query ${spp_ott} \
      --output_dir "\$(pwd)" \
      --label_format "name" \
      ${max_age} \
      ${phylo_only}

    """
}


// Merge filtered species occurrences and prep data for Biodiverse
// Add record counts (taking into account DBSCAN-based outliers)
process merge_occ {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/02.Biodiverse_input", mode: 'copy'
    // cpus 10

    input:
      path(spp, stageAs: "FilteredData/*")
      path phytree
      path counts

    output:
      path "H3_GridCell_Centres.csv", emit: h3coords
      path "Trimmed_occurrences.csv", emit: occurrences
      path "Trimmed_tree.nex", emit: tree
      path "Trimmed_tree.nwk", emit: treenwk

    script:
    """
    12_Prepare_Biodiverse_input.R \
      --input "FilteredData" \
      --phytree ${phytree} \
      --phylabels ${params.phylabels} \
      --taxgroup  ${params.taxgroup} \
      --counts    ${counts} \
      --threads   ${task.cpus} \
      --output    "."

    # --inputfile ${spp} \

    """

    // Nextflow gets the names of the input files from `spp` (`flt_ch` channel)
    // and copies/symlinks them to the current dir --> `--input` arg = current dir
}



// Prepare a list of DOIs for the GBIF datasets
// Acrivated by `--deriveddataset` param
process derived_datasets {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.tracedir", mode: 'copy'
    // cpus 1

    input:
      path(input)
      path(noextinct)
      path(polygon)
      path(wgsrpd)
      path(terrestrial)
      path(rmcountrycentroids)
      path(rmcountrycapitals)
      path(rminstitutions)
      path(rmurban)
      path(speciesinphylotree)

    output:
      path "Dataset_DOIs.txt", emit: doi

    script:

    filter_specieskeys  = params.specieskeys        ? "--specieskeys $specieskeys" : ""
    filter_extinct      = params.noextinct          ? "--noextinct $noextinct"     : ""
    filter_polygon      = params.polygon            ? "--polygon $polygon"         : ""
    filter_terrestrial  = params.terrestrial        ? "--terrestrial $terrestrial" : ""
    filter_country      = params.rmcountrycentroids ? "--rmcountrycentroids $rmcountrycentroids" : ""
    filter_capitals     = params.rmcountrycapitals  ? "--rmcountrycapitals $rmcountrycapitals"   : ""
    filter_institutions = params.rminstitutions     ? "--rminstitutions $rminstitutions" : ""
    filter_urban        = params.rmurban            ? "--rmurban $rmurban" : ""
    filter_wgsrpd       = params.wgsrpd             ? "--wgsrpd $wgsrpd --regions ${params.regions}" : ""

    """
    16_Derived_dataset.R \
      --input   ${input} \
      --phylum  ${params.phylum} \
      --class   ${params.classis} \
      --order   ${params.order} \
      --family  ${params.family} \
      --genus   ${params.genus} \
      --country ${params.country} \
      --latmin  ${params.latmin} \
      --latmax  ${params.latmax} \
      --lonmin  ${params.lonmin} \
      --lonmax  ${params.lonmax} \
      --minyear ${params.minyear} \
      --maxyear ${params.maxyear} \
      --coordprecision          ${params.coordprecision} \
      --coorduncertainty        ${params.coorduncertainty} \
      --coorduncertaintyexclude ${params.coorduncertaintyexclude} \
      --basisofrecordinclude ${params.basisofrecordinclude} \
      --basisofrecordexclude ${params.basisofrecordexclude} \
      ${filter_specieskeys} \
      --speciestree ${speciesinphylotree} \
      ${filter_extinct} \
      --excludehuman ${params.excludehuman} \
      ${filter_polygon} \
      ${filter_wgsrpd} \
      ${filter_terrestrial} \
      ${filter_country} \
      ${filter_capitals} \
      ${filter_institutions} \
      ${filter_urban} \
      --roundcoords ${params.roundcoords} \
      --resolution  ${params.h3resolution} \
      --threads     ${task.cpus} \
      --rcode       \$(which "Shapefile_filters.R") \
      --output      "Dataset_DOIs.txt"
    """
}



// // Split dataset based on spatial constraints
// process split_by_polygons {
//
//     label "container_r"
//     queue "custom_pool"
//
//     publishDir "$params.outdir/02.Biodiverse_input", mode: 'copy'
//     // cpus 2
//
//     input:
//       path occurrences
//       path polygons
//
//     output:
//       path "SpatConstrained_*.csv",  emit: occsplit
//
//     script:
//     """
//
//     12_Split_dataset_SpatialConstraints.R \
//       --input ${occurrences} \
//       --randconstrain ${polygons} \
//       --threads ${task.cpus} \
//       --output "SpatConstrained_"
//
//     """
// }


// Create Biodiverse input files
process prep_biodiv {

    label "container_biodiverse"
    queue "custom_pool"

    publishDir "$params.outdir/02.Biodiverse_input", mode: 'copy'

    // cpus 1

    input:
      path occurrences
      path tree

    output:
      path "occ.bds",          emit: BDS
      path "tree.bts",         emit: BTS
      path "occ_analysed.bds", emit: BDA
      path "occ.bds.csv",      emit: BDOBS

    script:
    """
    ## For debugging - check which Perl are we using?
    # perl --version

    ## Prepare Biodiverse input file
    ## NB! column numbers are zero-based here
    ## Latitude = Y, Longitude = X
    00_create_bds.pl \
      --csv_file ${occurrences} \
      --out_file "occ.bds" \
      --label_column_number     '0' \
      --sampcount_column_number '6' \
      --group_column_number_x   '5' \
      --group_column_number_y   '4' \
      --cell_size_x '0' \
      --cell_size_y '0'

    ## Prepare the tree for Biodiverse
    00_create_bts.pl \
      --input_tree_file ${tree} \
      --out_file "tree.bts"

    ## Run the analyses
    02_biodiverse_analyses.pl \
      --input_bds_file "occ.bds" \
      --input_bts_file "tree.bts" \
      --calcs ${params.indices}

    """
}



// Estimate phylogenetic diversity with Biodiverse
process phylodiv {

    label "container_biodiverse"
    queue "custom_pool"

    // publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'
    // cpus 1

    input:
      path(BDA)
      val(chunkid)

    output:
      path "Biodiv_randomized_${chunkid}.bds", emit: BDArand

    script:
    """
    03_run_randomisation.pl \
      --basedata ${BDA} \
      --bd_name  ${BDA} \
      --out_file "Biodiv_randomized.bds" \
      --rand_name 'rand' \
      --iterations ${iterations_per_thread} \
      --args ${biodiverse_args} \
      seed=${chunkid}

    ## Add chunk ID into the file name
    mv "Biodiv_randomized.bds" "Biodiv_randomized_${chunkid}.bds"

    """
}



// Prepare a shapefile to spatially constrain randomizations
process prep_shapefile {

    label "container_r"
    queue "custom_pool"

    // cpus 1

    input:
      path(occurrences)
      path(polygons)

    output:
      path "shapefile*", emit: shapefile

    script:
    """
    12_Prepare_SpatialConstraints.R \
      --input ${occurrences} \
      --randconstrain ${polygons} \
      --threads ${task.cpus} \
      --output "shapefile"
    """
}


// Estimate phylogenetic diversity using spatially-constrained randomizations
process phylodiv_constrianed {

    label "container_biodiverse"
    queue "custom_pool"

    // publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'
    // cpus 1

    input:
      path(BDA)
      path(polygons)
      val(chunkid)

    output:
      path "Biodiv_randomized_${chunkid}.bds", emit: BDArand

    script:
    """
    03_run_randomisation.pl \
      --basedata ${BDA} \
      --bd_name  ${BDA} \
      --out_file "Biodiv_randomized.bds" \
      --rand_name 'rand' \
      --iterations ${iterations_per_thread} \
      --args ${biodiverse_args} \
      seed=${chunkid} \
      spatial_conditions_for_subset='sp_points_in_same_poly_shape (file => "shapefile.shp")'

    ## Add chunk ID into the file name
    mv "Biodiv_randomized.bds" "Biodiv_randomized_${chunkid}.bds"

    """
}


// Create a file with paths to all chunks with randomization results
process rand_filelist {

    // container image is required for Cloud only
    label "container_r"
    queue "custom_pool"

    input:
      path(randfiles)

      //// To avoid name collisions if analysis was done for multiple datasets
      // path(randfiles, stageAs: "?/*")

    output:
      path "randomization_results.txt", emit: RND

    shell:
    $/
    echo "${randfiles}" \
      | sed -z 's/, /\n/g; s/^\[//; s/\]//' \
      > randomization_results.txt
    /$
}
// Groovy allows an alternative syntax for string definitions
// which uses the $ as escape character in place of \ character.
// These strings are delimited with an opening $/ and and a closing /$



// Aggregate the randomization results - with Biodiverse
process aggregate_rnds_biodiv {

    label "container_biodiverse"
    queue "custom_pool"

    publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'
    // cpus 1

    input:
      path RND
      path BDArand

    output:
      path "Biodiverse.bds", emit: Biodiv

    script:
    """
    05_reintegrate_basedatas_post_rand.pl \
      --glob ${RND} \
      --output_prefix Biodiverse
    """
}


// Export Biodiverse results into CSV
process div_to_csv {

    label "container_biodiverse"
    queue "custom_pool"

    publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'

    // cpus 1

    input:
      path Biodiv

    output:
      path "RND_groups.csv",                          emit: grp
      path "RND_SPATIAL_RESULTS.csv",                 emit: spat
      path "RND_rand--SPATIAL_RESULTS.csv",           emit: spat_r
      path "RND_rand--p_rank--SPATIAL_RESULTS.csv",   emit: spat_p
      path "RND_rand--z_scores--SPATIAL_RESULTS.csv", emit: spat_z
      path "RND_rand--CANAPE--.csv",                  emit: canape, optional: true
      path "RND_HURLBERT_ES.csv",                     emit: hurl,   optional: true
      path "RND_rand--HURLBERT_ES.csv",               emit: hurl_r, optional: true
      path "RND_rand--p_rank--HURLBERT_ES.csv",       emit: hurl_p, optional: true
      path "RND_rand--z_scores--HURLBERT_ES.csv",     emit: hurl_z, optional: true

    script:
    """
    04_load_bds_and_export_results.pl \
      --input_bds_file ${Biodiv} \
      --output_csv_prefix 'RND'

    """
}



// Plot Biodiverse results
process plot_pd {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/03.Plots", mode: 'copy'

    input:
      path BDOBS     // observed indices
      path spat_z    // randomized indices
      path world

    output:
      path "*.${params.plotformat}"

    script:

    add_world_map = params.world ? "--world $world" : ""

    """
    14_Visualization.R \
      --observed   ${BDOBS} \
      --zscores    ${spat_z} \
      --threads    ${task.cpus} \
      --variables  ${params.plotvar} \
      --resolution ${params.h3resolution} \
      ${add_world_map} \
      --format "${params.plotformat}" \
      --plotz  "${params.plottype}" \
      --width  "${params.plotwidth}" \
      --height "${params.plotheight}" \
      --units  "${params.plotunits}" \
      --output "."

    """
}


// Plot PD indices (interactive map - Leaflet-based choropleth)
process plot_leaflet {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/03.Plots", mode: 'copy'

    input:
      path BDOBS     // observed indices
      path spat_z    // randomized indices
      path spat_r    // randomization-based p-values
      path NRECORDS  // number of raw records per grid cell
      path canape    // CANAPE results
      path hurl      // Hurlbert ES indices

    output:
      path "Choropleth.html"
      path "Choropleth_files",              optional: true
      path "Biodiverse_results_merged.txt", optional: true, emit: tab
      path "Diversity_estimates.gpkg",      emit: geopackage
      path "Leaflet_object.RData",          emit: leaflet

    script:

    """
    15_Leaflets.R \
      --observed    ${BDOBS} \
      --sesscores   ${spat_z} \
      --sigscores   ${spat_r} \
      --canape      ${canape} \
      --hurlbert    ${hurl} \
      --reccounts   ${NRECORDS} \
      --resolution  ${params.h3resolution} \
      --variables   ${params.leaflet_var} \
      --canapesuper ${params.leaflet_canapesuper} \
      --palette     ${params.leaflet_palette} \
      --color       ${params.leaflet_color} \
      --bins        ${params.leaflet_bins} \
      --colorses    ${params.leaflet_sescolor} \
      --redundancy  ${params.leaflet_redundancy} \
      --output      "Choropleth.html"

    """
}


// Leaflet screenshots
process leaflet_mapshots {

    label "container_r"
    queue "custom_pool"

    publishDir "$params.outdir/03.Plots/Mapshots", mode: 'copy'
    tag "$ind"

    input:
      path leaf    // leaflet object
      val  ind     // Index to plot

    output:
      path "*.png", optional: true

    script:
    """

    15_Leaflets_mapshots.R \
      --input  ${leaf} \
      --metric ${ind}  \
      --allmetrics ${params.leaflet_var} \
      --zoom       ${params.leaflet_zoom} \
      --width      ${params.leaflet_width} \
      --height     ${params.leaflet_height}

    """
}




//  The default workflow
workflow {

    // Input directory with parquet files (GBIF dump dir)
    input_ch = Channel.value(params.input)

    // Run stage-I filtering
    occ_filter(
        input_ch,
        noextinct,
        polygon,
        wgsrpd
        )

    // occ_filter.out.part_low.view()
    // occ_filter.out.part_high.view()

    // Counting the total number of records
    record_count(
        input_ch,
        noextinct,
        polygon,
        wgsrpd,
        terrestrial,
        rmcountrycentroids,
        rmcountrycapitals,
        rminstitutions,
        rmurban
        )

    // Run stage-II filtering for species with low abundance (no DBSCAN)
    outl_low(
        occ_filter.out.part_low,
        polygon,
        wgsrpd,
        terrestrial,
        rmcountrycentroids,
        rmcountrycapitals,
        rminstitutions,
        rmurban
        )

    // Channel for DBSCAN-based filtering (iterate for each species)
    species_ch = occ_filter.out.spp.splitText().map{it -> it.trim()}
    // species_ch.view()

    // Run stage-II filtering for abundant species (with DBSCAN)
    outl_high(
        species_ch,
        occ_filter.out.part_high,
        polygon,
        wgsrpd,
        terrestrial,
        rmcountrycentroids,
        rmcountrycapitals,
        rminstitutions,
        rmurban
        )

    // Collect outlier removal statistics
    // out_high_ch = outl_high.out.outliers.collect()
    // out_low_ch = outl_low.out.outlierslow.collect()
    // Channel
    //     .from('out_high_ch', 'out_low_ch')
    //     .collectFile(name: 'OutlierCounts.txt', newLine: true, keepHeader: true)
    //     .subscribe {
    //         println "${it.text}"
    //     }

    // Automatically retrive phylogenetic tree from Open Tree of Life
    if(params.phytree == null){

      // Get OTT IDs for species
      prep_ott_ids(occ_filter.out.spp_all)

      // Get the tree
      get_ott_tree(prep_ott_ids.out.spp_ott)
      phytree = get_ott_tree.out.tree

    } else {
    // Use user-supplied hylogenetic tree
      phytree = file(params.phytree)
    }


    // Use output of for Biodiverse
    flt_ch = outl_high.out.sp.mix(outl_low.out.lowabsp).collect()
    // flt_ch.view()

    // Prepare a file with paths to occurrence filtering results
    // filtered_filelist(flt_ch)

    // Merge species occurrences into a single file
    // merge_occ(filtered_filelist.out.FLT)   // using inputfile (fail if there are a lot of species)
    merge_occ(
        flt_ch,
        phytree,
        record_count.out.n_spr)

    // Prepare Biodiverse input files
    prep_biodiv(merge_occ.out.occurrences, merge_occ.out.tree)

    // Channel of randomization chunks
    rnd_ch = Channel.fromList( randomization_chunks )


    // If no spatial constraints (for randomization) are provided
    if(params.randconstrain == null){

      // Perform unconstrained randomizations
      phylodiv(prep_biodiv.out.BDA, rnd_ch)

      // Prepare a file with paths of `phylodiv` output (multiple chunks of randomizations)
      // and create a new channel from it
      rand_filelist(phylodiv.out.BDArand.collect())

      // Aggregate randomization results (with Biodiverse script)
      aggregate_rnds_biodiv(
          rand_filelist.out.RND,
          phylodiv.out.BDArand.collect())

    } else {

      // If spatial constraints are provided, split dataset in parts (for each polygon)

      // A channel with spatial polygons
      polygons = file(params.randconstrain)

      // Prepare a shapefile with polygons
      prep_shapefile(merge_occ.out.occurrences, polygons)

      // Run spatially-constrained randomizations
      phylodiv_constrianed(
        prep_biodiv.out.BDA,
        prep_shapefile.out.shapefile,
        rnd_ch)

      // Collect `phylodiv_constrianed` output (multiple chunks of randomizations)
      rand_filelist(phylodiv_constrianed.out.BDArand.collect())

      // Aggregate randomization results (with Biodiverse script)
      aggregate_rnds_biodiv(
          rand_filelist.out.RND,
          phylodiv_constrianed.out.BDArand.collect())

    } // end of randomizations


    // // Split occurrences by polygons and run randomizations independently
    // polygons = file(params.randconstrain)
    //
    // // Split dataset
    // split_by_polygons(merge_occ.out.occurrences, polygons)
    //
    // // Channel with spatially-constrained datasets
    // ch_spatconstr = split_by_polygons.out.occsplit.flatten()
    // ch_spatconstr.view()
    //
    // // Prepare data for Biodiverse
    // prep_biodiv(ch_spatconstr, merge_occ.out.tree)
    //
    // // Channel with the number of randomization chunks
    // rnd_ch_tmp = Channel.fromList( randomization_chunks )
    //
    // // Apply randomization chunks for each spatially-contrained dataset (cartesian product)
    // rnd_ch = prep_biodiv.out.BDA.combine(rnd_ch_tmp)
    // rnd_ch.view()


    // Output results as CSV
    div_to_csv(aggregate_rnds_biodiv.out.Biodiv)

    // Plot PD indices (static map)
    plot_pd(
        div_to_csv.out.spat,
        div_to_csv.out.spat_z,
        world)

    // Channels with optional input for Leaflets
   ch_canape   = div_to_csv.out.canape.flatten().collect().ifEmpty(file("no_canape"))
   ch_hurlbert = div_to_csv.out.hurl.flatten().collect().ifEmpty(file("no_hurlbert"))

    // Plot PD indices (interactive map - Leaflet-based choropleth)
    plot_leaflet(
        div_to_csv.out.spat,
        div_to_csv.out.spat_z,
        div_to_csv.out.spat_r,
        record_count.out.n_recr,
        ch_canape,
        ch_hurlbert
        )


    // Prepare a channel with diversity index names
    ind_ch = Channel.value( params.leaflet_var )
      .splitText( { it.split(',') } )
      .distinct()
      .flatten()
      .map { it.trim() }

    // Leaflet screenshots
    // leaflet_mapshots(
    //   plot_leaflet.out.leaflet,
    //   ind_ch
    //   )

    // Prepare a list of DOIs for the GBIF datasets
    if(params.deriveddataset == true){
      derived_datasets(
        input_ch,
        noextinct,
        polygon,
        wgsrpd,
        terrestrial,
        rmcountrycentroids,
        rmcountrycapitals,
        rminstitutions,
        rmurban,
        merge_occ.out.occurrences
        )
    }

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
