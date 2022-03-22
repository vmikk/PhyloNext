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
