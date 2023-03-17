# Docker image with R packages required for GBIF occurrence filtering and mapping

FROM rocker/r-ver:4.2.2

LABEL org.opencontainers.image.authors="vladimir.mikryukov@ut.ee"

ENV LANG C.UTF-8
ENV SHELL /bin/bash

## 1. Install the basic dependencies
## 2. Install tidyverse packages along with arrow, data.table, and pandoc
## 3. Clean up
RUN apt-get update -qq \
    && apt-get -y --no-install-recommends install \
      zip unzip \
      curl git wget less \
      build-essential \
      libgeos-dev libudunits2-dev libproj-dev libgdal-dev \
      gnupg \
 && /rocker_scripts/install_tidyverse.sh \
 && /rocker_scripts/install_pandoc.sh \
 && apt-get purge --auto-remove -y gnupg \
 && apt-get autoremove -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/*

## x. Install Chrome (for webshot2 package), currently does not work well with Docker
#    && curl -sSL https://dl.google.com/linux/linux_signing_key.pub | gpg --dearmor | tee /etc/apt/trusted.gpg.d/chrome.gpg \
#    && echo "deb https://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google-chrome.list \
#    && apt-get update -qq \
#    && apt-get install -y --no-install-recommends google-chrome-stable \

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
    leaflet \
    mapview \
    webshot \
    webshot2 \
    htmlwidgets \
    R.utils \
    tinytest \
    covr \
    && R -e 'remotes::install_github("crazycapivara/h3-r")' \
    && R -e 'remotes::install_github("chgrl/leafletR")' \
    && R -e 'webshot::install_phantomjs()' \
    && rm -rf /tmp/downloaded_packages/

    # && R -e 'remotes::install_github("r-spatial/mapview")' \
    # && R -e 'remotes::install_github("rstudio/leaflet")' \
    # && R -e 'remotes::install_github("rstudio/webshot2")' \
    # && R -e 'remotes::install_github("rstudio/chromote")' \
    # && R -e 'remotes::install_github("ramnathv/htmlwidgets")' \

    ## For headless browser support: https://github.com/rstudio/chromote
    # chromote \

## Run bash in the container
# ENTRYPOINT ["bash"]
