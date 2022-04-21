#!/usr/bin/env nextflow
/*

========================================================================================
    GBIF phylogenetic diversity pipeline
========================================================================================
    Version: v0.0.1
    License: MIT
    Github : https://github.com/vmikk/biodiverse-scripts
    Website: TBA
    Slack  : TBA
----------------------------------------------------------------------------------------
*/

// TO DO:
// - autoname outdir (with current date?)
// - add a text file with parameter summary to the results
// - update visualization results (PDP, PD, S)
// - Vizualization: add channel with PD index names?
// - include sha256 hash of the image in the container reference: https://www.nextflow.io/blog/2016/docker-and-nextflow.html
// - add local and azure profiles
// - move params into a config
// - In Nextflow, any binaries within the the bin directory of a repository are added to the path
//   --> deprecate `scripts_path`
//   make scripts executable + remove Rscript and perl ??
// - Fix `10_Filter_occurrences.R`   In value[[3L]](cond) : double expected, got “NA”
// - Add spatial constraints for the randomizations (e.g., shapefile with biomes?)
// - Dynamic computing resources for intensive tasks: https://www.nextflow.io/docs/latest/process.html#dynamic-computing-resources
// - For species name matchin to the Open Tree Taxonomy - use pre-built DB?
// - Add alternative outlier detection algorithm (e.g., OPTICS algorithm)
// - Add spatial constraints for randomizations (e.g., within biome)
// - Split the pipeline into workflows?
// - Colorize verbose messages


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.0.1'

//// Initialize parameters, set default values

params.scripts_path = "${projectDir}/bin"
params.data_path = "${projectDir}/pipeline_data"

// Filtering, stage I - "10_Filter_occurrences.R"
params.input = false
params.outdir = "${launchDir}/results"
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

// Biodiverse
params.indices = "calc_richness,calc_simpson_shannon,calc_endemism_whole,calc_pd,calc_pe,calc_phylo_rpd1,calc_phylo_rpd2,calc_phylo_rpe1,calc_phylo_rpe2,calc_phylo_mpd_mntd2"
params.randname = "rand_structured"
params.iterations = 100
params.biodiversethreads = 1

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



// Visualization
params.plotvar = "PHYLO_RPD1"
params.plotformat = "pdf"
params.world = params.data_path + "/WorldMap_NaturalEarth_Medium.RData"

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
out_plot = params.outdir + "/03.Plots"
out_logs = params.outdir + "/logs"


// Pipeline help message
def helpMsg() {
    log.info"""
    ====================================================================
    GBIF phylogenetic diversity pipeline :  Version ${version}
    ====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run vmikk/biodiverse-scripts -r ${version} --input ... --outdir ...
    
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
        --indices             Comma-seprated list of diversity and endemism indices; e.g., "calc_richness,calc_pd,calc_pe"
        --randname            Randomisation scheme type; e.g., "rand_structured"
        --iterations          Number of randomisation iterations; e.g., 1000
    NEXTFLOW-SPECIFIC:
        -qs                   Queue size (max number of processes that can be executed in parallel); e.g., 8
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
        GBIF phylogenetic diversity pipeline: Version ${version}
        ========================================================
        GBIF occurrence dump:     ${params.input}
        Output path:              ${params.outdir}
        H3 spatial resolution:    ${params.h3resolution}
        Land-mask:                ${params.terrestrial}
        Spatial outliers removal: ${params.dbscan}
        """
        .stripIndent()

if(params.dbscan == true){
    log.info "  DBSCAN epsilon:           ${params.dbscanepsilon}".stripIndent()
    log.info "  DBSCAN minpts:            ${params.dbscanminpts}".stripIndent()
}

log.info """
        Biodiverse params:
          Biodiverse indices:         ${params.indices}
          Biodiverse randomizations:  ${params.iterations}
          Biodiverse threads:         ${params.biodiversethreads}
          Biodiverse rand per thread: ${iterations_per_thread}
          Biodiverse args:            ${biodiverse_args}
        """
        .stripIndent()

log.info """
        Pipeline info:
          Pipeline version:                ${version}
          Pipeline scripts location:       ${params.scripts_path}
          Pipeline internal data location: ${params.data_path}
          Pipeline profile:                ${workflow.profile}
          Config file used:                ${workflow.configFiles}
          Container engine:                ${workflow.containerEngine}
        """
        .stripIndent()

log.info """
        Core Nextflow options:
          launchDir:              ${workflow.launchDir}
          workDir:                ${workflow.workDir}
          projectDir:             ${workflow.projectDir}
        """
        .stripIndent()

log.info "========================================================"
log.info "\n"


// Occurrence filtering, stage I
process occ_filter {

    label "container_r"
    containerOptions = { workflow.containerEngine == "docker" ? 
        "--volume ${params.input}:${params.input} --volume ${params.outdir}:${params.outdir}"
        : null }

    publishDir "$params.outdir", mode: 'copy'
    // cpus 10

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

    label "container_r"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir} --volume ${params.data_path}:${params.data_path}"
        : null }

    // publishDir "$params.outdir/01.filtered2", mode: 'copy'
    // cpus 5

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


// Outlier filtering, stage II - with DBSCAN, independently by species
process outl_high {

    label "container_r"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir} --volume ${params.data_path}:${params.data_path}"
        : null }

    // publishDir "$params.outdir/01.filtered2", mode: 'copy'
    // cpus 1

    // Add species ID to the log file
    tag "$sp"

    // Iterate for each species from `species_ch` channel
    input:
      val sp

    output:
      path "${sp}.RData", emit: sp

    script:
    """
    Rscript ${params.scripts_path}/11_Additional_filtering_and_aggregation.R \
      --input "${out_flt1}/Partition=high" \
      --specieskey ${sp} \
      --dbscan ${params.dbscan} \
      --epsilon ${params.dbscanepsilon} \
      --minpts ${params.dbscanminpts} \
      --resolution ${params.h3resolution} \
      --terrestrial ${params.terrestrial} \
      --threads ${task.cpus} \
      --output ${out_flt2}

    cp ${out_flt2}/${sp}.RData ${sp}.RData

    """
}


// Merge filtered species occurrences and prep data for Biodiverse
process merge_occ {

    label "container_r"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir} --volume ${params.phytree}:${params.phytree}"
        : null }

    publishDir "$params.outdir/02.Biodiverse_input", mode: 'copy'
    // cpus 10

    input:
      val spp

    output:
      path "H3_GridCell_Centres.csv", emit: h3coords
      path "Trimmed_occurrences.csv", emit: occurrences
      path "Trimmed_tree.nex", emit: tree

    script:
    """
    Rscript ${params.scripts_path}/12_Prepare_Biodiverse_input.R \
      --input ${out_flt2} \
      --phytree ${params.phytree} \
      --taxgroup ${params.taxgroup} \
      --threads ${task.cpus} \
      --output  "\$PWD"  # ${out_biod}

    """
}

// Create Biodiverse input files
process prep_biodiv {

    label "container_biodiverse"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir}"
        : null }

    publishDir "$params.outdir/02.Biodiverse_input", mode: 'copy'

    // cpus 1

    input:
      val occurrences
      val tree

    output:
      path "occ.bds", emit: BDS
      path "tree.bts", emit: BTS
      path "occ_analysed.bds", emit: BDA
      path "occ.bds.csv", emit: BDOBS

    script:
    """
 
    ## For debugging - check which Perl are we using?
    # perl --version

    ## Prepare Biodiverse input file
    perl ${params.scripts_path}/00_create_bds.pl \
      --csv_file ${occurrences} \
      --out_file "occ.bds" \
      --label_column_number '0' \
      --group_column_number_x '1' \
      --cell_size_x '-1'
    
    ## Prepare the tree for Biodiverse
    perl ${params.scripts_path}/00_create_bts.pl \
      --input_tree_file ${tree} \
      --out_file "tree.bts"
    
    ## Run the analyses
    perl ${params.scripts_path}/02_biodiverse_analyses.pl \
      --input_bds_file "occ.bds" \
      --input_bts_file "tree.bts" \
      --calcs ${params.indices}
 
    # cp "occ.bds.csv" "${params.outdir}/02.Biodiverse_input/Biodiverse_ObservedIndices.csv"
    # 
    # cp "occ.bds"  "${params.outdir}/02.Biodiverse_input/occ.bds"
    # cp "tree.bts" "${params.outdir}/02.Biodiverse_input/tree.bts"
    # cp "occ_analysed.bds" "${params.outdir}/02.Biodiverse_input/occ_analysed.bds"

    """
}



// Estimate phylogenetic diversity with Biodiverse
process phylodiv {

    label "container_biodiverse"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir}"
        : null }

    // publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'
    // cpus 1

    input:
      val BDA
      val chunkid

    output:
      path "Biodiv_randomized.bds", emit: BDArand

    script:
    """

    perl ${params.scripts_path}/03_run_randomisation.pl \
      --basedata ${BDA} \
      --bd_name ${BDA} \
      --out_file "Biodiv_randomized.bds" \
      --rand_name 'rand' \
      --iterations ${params.iterations} \
      --args ${biodiverse_args}

    """
}




// Create a file with paths to all chunks with randomization results
process rand_filelist {

    input:  val randfiles
    output: path "randomization_results.txt", emit: RND

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
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir}"
        : null }

    publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'
    // cpus 1
    
    input:
      val RND

    output:
      path "Biodiverse.bds", emit: Biodiv

    script:
    """
    perl ${params.scripts_path}/05_reintegrate_basedatas_post_rand.pl \
      --glob ${RND} \
      --output_prefix Biodiverse
    """
}


// Export Biodiverse results into CSV
process div_to_csv {

    label "container_biodiverse"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir}"
        : null }

    publishDir "$params.outdir/02.Biodiverse_results", mode: 'copy'

    // cpus 1

    input:
      val Biodiv

    output:
      path "RND_groups.csv", emit: RND1
      path "RND_rand--p_rank--SPATIAL_RESULTS.csv", emit: RND2
      path "RND_rand--SPATIAL_RESULTS.csv", emit: RND3
      path "RND_rand--z_scores--SPATIAL_RESULTS.csv", emit: RND4
      path "RND_SPATIAL_RESULTS.csv", emit: RND5

    script:
    """

    perl ${params.scripts_path}/04_load_bds_and_export_results.pl \
      --input_bds_file ${Biodiv} \
      --output_csv_prefix 'RND'

    """
}



// Plot Biodiverse results
process plot_pd {

    label "container_r"
    containerOptions = { workflow.containerEngine == "docker" ?
        "--volume ${params.outdir}:${params.outdir} --volume ${params.data_path}:${params.data_path}"
        : null }

    publishDir "$params.outdir/03.Plots", mode: 'copy'

    input:
      val BDOBS   // observed indices
      val RND4    // randomized indices

    output:
      path "*.${params.plotformat}"

    script:
    """
    Rscript ${params.scripts_path}/14_Visualization.R \
      --observed ${BDOBS} \
      --zscores ${RND4} \
      --threads ${task.cpus} \
      --variables ${params.plotvar} \
      --world "${params.world}" \
      --format "${params.plotformat}" \
      --output "."

      # --output "${out_plot}"

    """
}




//  The default workflow
workflow {

    // Input directory with parquet files (GBIF dump dir)
    input_ch = Channel.value(params.input)

    // Run stage-I filtering
    occ_filter(input_ch)

    // Channel with land shapefile
    // land_ch = Channel.value(params.terrestrial)

    // Run stage-II filtering for species with low abundance (no DBSCAN)
    outl_low( occ_filter.out.part_low )

    // Channel for DBSCAN-based filtering (iterate for each species)
    //species_ch = Channel.fromPath("${params.outdir}/spp.txt").splitText().map{it -> it.trim()}
    species_ch = occ_filter.out.spp.splitText().map{it -> it.trim()}
    // species_ch.view()

    // Run stage-II filtering for abundant species (with DBSCAN)
    outl_high(species_ch)

    // Use output of for Biodiverse
    flt_ch = outl_high.out.sp.mix(outl_low.out.lowabsp).collect()
    // flt_ch.view()

    // Merge species occurrences into a single file
    merge_occ(flt_ch)

    // Prepare Biodiverse input files
    prep_biodiv(merge_occ.out.occurrences, merge_occ.out.tree)


    // Channel of randomization chunks
    rnd_ch = Channel.fromList( randomization_chunks )

    // Perform Biodiverse randomizations
    phylodiv(prep_biodiv.out.BDA, rnd_ch)

    // Prepare a file with paths of `phylodiv` output (multiple chunks of randomizations)
    // and create a new channel from it
    rand_filelist(phylodiv.out.BDArand.collect())

    // Aggregate randomization results (with Biodiverse script)
    aggregate_rnds_biodiv(rand_filelist.out.RND)

    // Output results as CSV
    div_to_csv(aggregate_rnds_biodiv.out.Biodiv)

    // Plot PD indices
    plot_pd(div_to_csv.out.RND5, div_to_csv.out.RND4)
    // plot_pd(prep_biodiv.out.BDOBS, div_to_csv.out.RND4)


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
