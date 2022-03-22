# PD (Phylogenetic Diversity) in the cloud
Scripts for Biodiverse pipeline
## Introduction

Current pipeline brings together two critical research data infrastructures, the Global
Biodiversity Information Facility [(GBIF)](https://www.gbif.org/) and Open Tree of Life [(OToL)](https://tree.opentreeoflife.org), to make them more accessible to non-experts.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. 

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/)


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
    --indices             Comma-seprated list of diversity and endemism indices; e.g., "calc_richness,calc_pd,calc_pe"
    --randname            Randomisation scheme type; e.g., "rand_structured"
    --iterations          Number of randomisation iterations; e.g., 1000

```
