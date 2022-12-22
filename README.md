# PhyloNext - PD (Phylogenetic Diversity) in the cloud <img src='images/PhyloNext_logo.png' align="right" height="100" />

[![Nextflow](https://img.shields.io/badge/Nextflow%20DSL2-%E2%89%A522.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-blue?style=flat&logo=singularity)](https://sylabs.io/docs/)
[![Github_Status_Badge](https://img.shields.io/badge/GitHub-0.0.2-blue.svg)](https://github.com/vmikk/PhyloNext)
[![GitHub license](https://img.shields.io/github/license/vmikk/PhyloNext)](https://github.com/vmikk/PhyloNext/blob/main/LICENSE)

The automated pipeline for the analysis of phylogenetic diversity using [GBIF occurrence data](https://www.gbif.org/occurrence/search?occurrence_status=present), species phylogenies from [Open Tree of Life](https://tree.opentreeoflife.org), and [Biodiverse software](https://shawnlaffan.github.io/biodiverse/).

## Introduction

Current pipeline brings together two critical research data infrastructures, the Global
Biodiversity Information Facility [(GBIF)](https://www.gbif.org/) and Open Tree of Life [(OToL)](https://tree.opentreeoflife.org), to make them more accessible to non-experts.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses [Docker](https://www.docker.com/) containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

The pipeline could be launched in a cloud environment (e.g., the [Microsoft Azure Cloud Computing Services](https://azure.microsoft.com/en-us/), [Amazon AWS Web Services](https://aws.amazon.com/), and [Google Cloud Computing Services](https://cloud.google.com/)).

## Pipeline summary

1. Filtering of GBIF species occurrences for various taxonomic clades and geographic areas
2. Removal of non-terrestrial records and spatial outliers (using density-based clustering)
3. Preparation of phylogenetic tree (currently, only pre-constructed phylogenetic trees are available; with the update of OToL, phylogenetic trees will be downloaded automatically using API) and name-matching with GBIF species keys
4. Spatial binning of species occurrences using Uber’s H3 system (hexagonal hierarchical spatial index)
5. Estimation of phylogenetic diversity and endemism indices using [Biodiverse program](https://shawnlaffan.github.io/biodiverse/)
6. Visualization of the obtained results

## Quick Start

An example command to run the pipilene:

```bash
nextflow run vmikk/phylonext -r main \
  --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet/" \
  --classis "Mammalia" --family  "Felidae,Canidae" \
  --country "DE,PL,CZ"  \
  --minyear 2000  \
  --dbscan true  \
  --phytree $(realpath "${HOME}/.nextflow/assets/vmikk/phylonext/test_data/phy_trees/Mammals.nwk") \
  --iterations 100  \
  -resume
```

## Documentation

The PhyloNext pipeline comes with documentation about the pipeline usage 
at [https://phylonext.github.io/](https://phylonext.github.io/).  

Main pipeline parameters and output are desribed here:
- [parameters](https://phylonext.github.io/parameters/)
- [output](https://phylonext.github.io/outputs/)

To show a help message, run `nextflow run vmikk/phylonext -r main --helpMsg`.
```
=====================================================================
PhyloNext: GBIF phylogenetic diversity pipeline :  Version 1.0
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

    --country             Country code, ISO 3166 (multiple comma-separated values allowed); e.g., "DE,PL,CZ"
    --latmin              Minimum latitude of species occurrences (decimal degrees); e.g., 5.1
    --latmax              Maximum latitude of species occurrences (decimal degrees); e.g., 15.5
    --lonmin              Minimum longitude of species occurrences (decimal degrees); e.g., 47.0
    --lonmax              Maximum longitude of species occurrences (decimal degrees); e.g., 55.5
    --minyear             Minimum year of record's occurrences; default, 1945
    --basisofrecordinclude Basis of record to include from the data; e.g., "PRESERVED_SPECIMEN"
    --basisofrecordexclude Basis of record to exclude from the data; e.g., "FOSSIL_SPECIMEN,LIVING_SPECIMEN"
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

Leaflet interactive visualization:
    --leaflet_var         Variables to plot; e.g., "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE,PE_WE,SES_PE_WE,CANAPE,Redundancy"
    --leaflet_color       Color scheme for continuous variables (default, "RdYlBu")
    --leaflet_palette     Color palette for continuous variables (default, "quantile")
    --leaflet_bins        Number of color bins for continuous variables (default, 5)
    --leaflet_redundancy  Redundancy threshold for hiding the grid cells with low number of records (default, 0 = display all grid cells)

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
```

## Credits

PhyloNext pipeline was developed by Vladimir Mikryukov and Kessy Abarenkov.

[Biodiverse program](https://shawnlaffan.github.io/biodiverse/) and Perl scripts accompanying PhyloNext were written by [Shawn Laffan](https://github.com/shawnlaffan) (Laffan et al., 2010).

Scripts for getting an induced subtree from the Open Tree of Life were developed by [Emily Jane McTavish](https://github.com/snacktavish).

We thank the following people for their extensive assistance in the development of this pipeline: Joe Miller, Shawn Laffan, Tim Robertson, Emily Jane McTavish, John Waller, Thomas Stjernegaard Jeppesen, and Matthew Blissett.

Also we are very grateful to [Manuele Simi](https://github.com/manuelesimi) and [nf-core](https://nf-co.re/) community for helpful advices on the development of this pipeline.

For more details, please see the [Acknowledgments section](https://phylonext.github.io/acknowledgements/) in the docs.

## Funding

The work is supported by a grant “PD (Phylogenetic Diversity) in the Cloud” to GBIF Supplemental funds from the GEO-Microsoft Planetary Computer Programme.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](CONTRIBUTING.md).

For further information or help, don't hesitate to file an [issue on GitHub](https://github.com/vmikk/PhyloNext/issues).

## Future plans

- Add support of [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) containers.

## Citations

If you use PhyloNext pipeline for your analysis, please cite it using the following DOI: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

Laffan SW, Lubarsky E, Rosauer DF (2010) Biodiverse, a tool for the spatial analysis of biological and related diversity. Ecography, 33: 643-647. [DOI: 10.1111/j.1600-0587.2010.06237.x](https://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2010.06237.x)

An extensive list of references for the tools used by the pipeline can be found in the [Citations](https://phylonext.github.io/citations/) section in the documentation.

