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

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/)

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run vmikk/phylonext -r main -profile test,docker
    ```
4. Start running your own analysis!

    ```console
    nextflow run vmikk/phylonext -r main \
      --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet/" \
      --classis "Mammalia" --family  "Felidae,Canidae" \
      --country "DE,PL,CZ"  \
      --minyear 2000  \
      --dbscan true  \
      --phytree  $(realpath "${HOME}/.nextflow/assets/vmikk/phylonext/test_data/phy_trees/Mammals.nwk") \
      --iterations 100  \
      --outdir "$PWD" \
      -resume
    ```

### Installation example on Ubuntu

1. Nextflow installation:

    Nextflow requires Java 11 (or later, up to 18) to be installed.
    ```
    sudo apt-get update
    sudo apt-get install default-jdk
    ```
    Install Nextflow:
    ```
    wget -qO- https://get.nextflow.io | bash
    chmod +x ./nextflow
    mkdir -p ~/bin & mv ./nextflow ~/bin/
    ```

2. Docker installation (for details see [the official Docker documentation](https://docs.docker.com/engine/install/ubuntu/)):
    ```
    sudo apt-get update
    sudo apt-get install ca-certificates curl gnupg lsb-release

    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

    echo \
      "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
      $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

    sudo apt-get update
    sudo apt-get install docker-ce docker-ce-cli containerd.io

    sudo usermod -aG docker $USER
    newgrp docker
    ```

### Singularity

As Docker is NOT supported in most HPC systems, it's possible to run the pipeline using [Singularity](https://sylabs.io/singularity/).


## Documentation
To show a help message, run `nextflow run vmikk/phylonext -r main --helpMsg`.
```
====================================================================
GBIF phylogenetic diversity pipeline :  Version 0.0.1
====================================================================

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
```

#### Passing in an input parameter file

It is possible to pass the pipeline parameters via YAML or JSON file, e.g.:
```
nextflow run vmikk/phylonext -r main -resume -params-file Mammals.yaml
```
The YAML file could contain the following:
```
input      : "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet/"
class      : "Mammalia"
family     : "Felidae,Canidae"
country    : "DE,PL,CZ"
minyear    : 2000
dbscan     : true
phytree    : "/path/to/the/phylogenetic/tree.nwk"
iterations : 100
outdir     : "${launchDir}"
```


#### The other helpful commands

```
## Download or update the pipeline
## By default, the pipeiline is stored in the '~/.nextflow/assets/vmikk/PhyloNext' directory
nextflow pull vmikk/phylonext

## Run the latest development version of the pipeline
nextflow run vmikk/phylonext -r main ...

## Run the tagged version (e.g., v0.1) of the pipeline
nextflow run vmikk/phylonext -r v0.1 ...

## Print the pipeline and system runtime information
nextflow info
nextflow info vmikk/phylonext

## Delete the local copy of the pipeline
nextflow drop vmikk/phylonext
```

#### Obtaining a local snapshot of species occurrences from GBIF

If you would like to run the pipeline locally, you may download the GBIF occurrence dump from the [Amazon AWS cloud](https://registry.opendata.aws/gbif/) using the [AWS CLI program](https://aws.amazon.com/cli/) (No AWS account required). E.g., to download `2022-05-01` dump to the `GBIF` directory in your home folder run:
```
aws s3 sync \
  s3://gbif-open-data-eu-central-1/occurrence/2022-05-01/occurrence.parquet/ \
  ~/GBIF/ \
  --no-sign-request
```

Alternative GBIF snapshot is also [hosted by the Microsoft AI for Earth program](https://github.com/microsoft/AIforEarthDataSets/blob/main/data/gbif.md). To download it using the [AzCopy command-line utility](https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10) run:
```
azcopy copy \
  "https://ai4edataeuwest.blob.core.windows.net/gbif/occurrence/2022-05-01/occurrence.parquet/*" \
  "~/GBIF/" \
  --recursive=true
```





## Funding

The work is supported by a grant “PD (Phylogenetic Diversity) in the Cloud” to GBIF Supplemental funds from the GEO-Microsoft Planetary Computer Programme.
## Future plans

- Add support of [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) containers.

