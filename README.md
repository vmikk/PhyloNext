# PD (Phylogenetic Diversity) in the cloud

[![Nextflow](https://img.shields.io/badge/Nextflow%20DSL2-%E2%89%A521.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![Github_Status_Badge](https://img.shields.io/badge/GitHub-0.0.1-blue.svg)](https://github.com/vmikk/biodiverse-scripts)
[![GitHub license](https://img.shields.io/github/license/vmikk/biodiverse-scripts)](https://github.com/vmikk/biodiverse-scripts/blob/main/LICENSE)

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
5. Estimation of phylogenetic diversity and endemism indices using [Biodeverse program](https://shawnlaffan.github.io/biodiverse/)
6. Visualization of the obtained results (to be implemented soon)


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/)

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run vmikk/biodiverse-scripts -r main -profile test
    ```
4. Start running your own analysis!

    ```console
    nextflow run vmikk/biodiverse-scripts -r main \
      --input "/mnt/GBIF/Parquet/2022-01-01/occurrence.parquet/" \
      --class "Mammalia" --family  "Felidae,Canidae" \
      --country "DE,PL,CZ"  \
      --minyear 2000  \
      --dbscan true  \
      --phytree  $(realpath "${HOME}/.nextflow/assets/vmikk/biodiverse-scripts/test_data/phy_trees/Mammals.nwk") \
      --iterations 100  \
      --outdir "$PWD" \
      -resume
    ```

### Installation example on Ubuntu

1. Nextflow installation:

    Nextflow requires Java 8 (or later, up to 17) to be installed.
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


## Documentation
To show a help message, run `nextflow run vmikk/biodiverse-scripts -r main --helpMsg`.
```
====================================================================
GBIF phylogenetic diversity pipeline :  Version 0.0.1
====================================================================

Pipeline Usage:
To run the pipeline, enter the following in the command line:
    nextflow run vmikk/biodiverse-scripts -r main --input ... --outdir ...

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

```

The other helpful commands:
```
## Download or update the pipeline
## By default, the pipeiline is stored in the '~/.nextflow/assets/vmikk/biodiverse-scripts' directory
nextflow pull vmikk/biodiverse-scripts

## Run the latest development version of the pipeline
nextflow run vmikk/biodiverse-scripts -r main ...

## Run the tagged version (e.g., v0.1) of the pipeline
nextflow run vmikk/biodiverse-scripts -r v0.1 ...

## Print the pipeline and system runtime information
nextflow info
nextflow info vmikk/biodiverse-scripts

## Delete the local copy of the pipeline
nextflow drop vmikk/biodiverse-scripts
```


If you would like to run the pipeline locally, you may download the GBIF occurrence dump from the [Amazon AWS cloud](https://registry.opendata.aws/gbif/) using the [AWS CLI program](https://aws.amazon.com/cli/) (No AWS account required). E.g., to download `2022-01-01` dump to the `GBIF` directory in your home folder run:
```
aws s3 sync \
  s3://gbif-open-data-eu-central-1/occurrence/2022-01-01/occurrence.parquet/ \
  ~/GBIF/ \
  --no-sign-request
```


