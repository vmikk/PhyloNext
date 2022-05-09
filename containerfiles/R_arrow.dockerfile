# Docker image with R packages required for GBIF occurrence filtering

FROM rocker/r-ver:4.1.2

MAINTAINER vladimir.mikryukov@ut.ee

ENV LANG C.UTF-8
ENV SHELL /bin/bash

## Install the required dependencies
RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    zip unzip \
    curl git wget less \
    build-essential \
    libgeos-dev libudunits2-dev libproj-dev libgdal-dev

## Install tidyverse packages along with arrow and data.table
RUN /rocker_scripts/install_tidyverse.sh \
  && apt-get autoremove -y \
  && apt-get autoclean -y \
  && rm -rf /var/lib/apt/lists/*

## Install additional R packages
RUN install2.r --error --skipinstalled --ncpus -1 \
    optparse \
    ape \
    plyr \
    rotl \
    rgbif \
    proj4 \
    rgeos \
    sf \
    dbscan \
    doFuture \
    && R -e 'remotes::install_github("crazycapivara/h3-r")' \
    && rm -rf /tmp/downloaded_packages/

## Install arrow
# RUN R -e 'arrow::install_arrow(minimal = FALSE)' \
#     && rm -rf /tmp/downloaded_packages/

## Run bash in the container
# ENTRYPOINT ["bash"]
