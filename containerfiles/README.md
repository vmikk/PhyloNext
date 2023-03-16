# PhyloNext containers

Currently, all software dependencies used by PhyloNext are grouped into 3 containers:
- `rarrow` - R
- `biodiverse` - Biodiverse and Perl
- `opentree` - Open Tree of Life and Python

## `rarrow` changelog

A Docker image with R packages required for GBIF occurrence filtering and mapping.

### v.1.2.0 (Mar 15, 2023)

- removed Google Chrome
  (headless browser and mapshots does not work well inside a Docker container)

- `R` 4.2.2
- `GEOS` 3.10.2
- `GDAL` 3.4.1
- `PROJ` 8.2.1

- `ape` 5.7-1
- `arrow` 11.0.0.3
- `dbscan` 1.1.11
- `h3` 3.7.2
- `leaflet` 2.1.2
- `rgbif` 3.7.5
- `rotl` 3.0.14
- `sf` 1.0.10
- `terra` 1.7.18

Docker compressed image size ~ 706.46 MB
Docker decompressed image size ~ 2.02 GB
Singularity image size ~ 694.12 MB

Docker image digest: `sha256:23bd794c33e94c6806017a8e4f1bbfb037cda803dadbf68c055cdc65757d2a43`
Singularity unique ID: `sha256.73c7e453f4f306b1ce5b6bc607a82f407b769252b0edc7e6b599b750e83e8551`

### v.1.1.0 (Feb 01, 2023)

- added Google Chrome and `chromote` package (for `webshot2` package)

- `R` 4.2.2
- `GEOS` 3.10.2
- `GDAL` 3.4.1
- `PROJ` 8.2.1

- `ape` 5.6-2
- `arrow` 10.0.1
- `dbscan` 1.1-11
- `h3` 3.7.2
- `leaflet` 2.1.1
- `rgbif` 3.7.5
- `rotl` 3.0.14
- `sf` 1.0-9
- `terra` 1.7-3

Compressed Size ~ 862.6 MB
Decompressed Size ~ 2.42 GB

Docker image digest: `sha256:9f5662b5f9b6c450b0dd24f60e026670b72ccfda1276995f0aa476f1bd5bf23f`

### v.1.0.0 (Jan 20, 2023)

- added `pandoc`, `R.utils`, `tinytest`, `covr`

- `R` 4.2.2
- `GEOS` 3.10.2
- `GDAL` 3.4.1
- `PROJ` 8.2.1

- `ape` 5.6-2
- `arrow` 10.0.1
- `dbscan` 1.1-11
- `h3` 3.7.2
- `leaflet` 2.1.1
- `rgbif` 3.7.4
- `rotl` 3.0.14
- `sf` 1.0-9
- `terra` 1.6-47

Compressed Size ~ 756.7 MB
Decompressed Size ~ 2.15 GB

Docker image digest: `sha256:4e9e60ab9973a41808e12b2d348d794949061020cf0135d2fba49f9e5af09254`
Singularity unique ID: `sha256.35985310c3c455d02586064b5e7e09165350ea2ec4bbbfc7067b3b0b264790f5`

### v.0.0.2 (Jul 08, 2022)

- added leaflet with dependencies

- `R` 4.2.1
- `GEOS` 3.8.0
- `GDAL` 3.0.4
- `PROJ` 6.3.1

- `ape` 5.6-2
- `arrow` 8.0.0
- `dbscan` 1.1-10
- `h3` 3.7.1
- `leaflet` 2.1.1
- `rgbif` 3.7.2
- `rotl` 3.0.12
- `sf` 1.0-7
- `terra` 1.5-34

Compressed Size ~ 807.06 MB
Decompressed Size ~ 2.4 GB

### v.0.0.1 (Apr 21, 2022)

- `R` 4.1.2
- `GEOS` 3.8.0
- `GDAL` 3.0.4
- `PROJ` 6.3.1

- `ape` 5.6-2
- `arrow` 7.0.0
- `dbscan` 1.1-10
- `h3` 3.7.1
- `rgbif` 3.7.1
- `rotl` 3.0.12
- `sf` 1.0-7

Compressed Size ~ 819.36 MB
Decompressed Size ~ 2.64 GB


