# Docker image with R packages required for GBIF occurrence filtering

FROM rocker/r-ver:4.2.2

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

## Install tidyverse packages along with arrow and data.table, + pandoc
RUN /rocker_scripts/install_tidyverse.sh \
  && /rocker_scripts/install_pandoc.sh \
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
    webshot \
    R.utils \
    tinytest \
    covr \
    && R -e 'remotes::install_github("crazycapivara/h3-r")' \
    && R -e 'remotes::install_github("chgrl/leafletR")' \
    && R -e 'remotes::install_github("r-spatial/mapview")' \
    && R -e 'webshot::install_phantomjs()' \
    && rm -rf /tmp/downloaded_packages/

## Install arrow
# RUN R -e 'arrow::install_arrow(minimal = FALSE)' \
#     && rm -rf /tmp/downloaded_packages/

## Run bash in the container
# ENTRYPOINT ["bash"]
