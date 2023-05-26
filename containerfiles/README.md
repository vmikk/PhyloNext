# PhyloNext containers

Currently, all software dependencies used by PhyloNext are grouped into 3 containers:  
- `rarrow` - R  
- `biodiverse` - Biodiverse and Perl  
- `opentree` - Open Tree of Life and Python  

For easy access and setup, the Docker and Singularity image archives for the PhyloNext pipeline are available on Zenodo:  
- `v.1.3.0`: [DOI:10.5281/zenodo.7973798](https://zenodo.org/record/7973798)  


## `rarrow` changelog

A Docker image with R packages required for GBIF occurrence filtering and mapping.  
Docker hub: [https://hub.docker.com/r/vmikk/rarrow/tags](https://hub.docker.com/r/vmikk/rarrow/tags)  
Singularity Library: [https://cloud.sylabs.io/library/vmiks/gbif/rarrow](https://cloud.sylabs.io/library/vmiks/gbif/rarrow)  


### v.1.3.0 (May 23, 2023)

- `R` 4.3.0
- `GEOS` 3.10.2
- `GDAL` 3.4.1
- `PROJ` 8.2.1

- `ape` 5.7-1
- `arrow` 12.0.0
- `dbscan` 1.1.11
- `h3` 3.7.2
- `leaflet` 2.1.2
- `rgbif` 3.7.7
- `rotl` 3.0.14
- `sf` 1.0.12
- `terra` 1.7.29

Docker compressed image size ~ 740.81 MB  
Docker decompressed image size ~ 2.09 GB  
Singularity image size ~ 697.93 MB  

Docker image digest: `sha256:f9d8270c4ccf921f23f8748826a4892be5b18810addc6a9fb95af44bc4478630`  
Singularity unique ID: `sha256.9d06bb52f10bc1a38fe8552a058721324525d8b036cf6fc704496faa8ac686d9`  


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



## `biodiverse` changelog

A Docker image with `Biodiverse` and `Perl` for estimation of phylogenetic diversity.  
Docker hub: [https://hub.docker.com/r/vmikk/biodiverse/tags](https://hub.docker.com/r/vmikk/biodiverse/tags)  
Singularity Library: [https://cloud.sylabs.io/library/vmiks/gbif/biodiverse](https://cloud.sylabs.io/library/vmiks/gbif/biodiverse)  


### v.1.3.0 (May 23, 2023)

- Biodiverse v4.3 `76fd253` (https://github.com/vmikk/biodiverse-docker/releases/tag/v.1.3.0)  

Biodiverse changelog:  
- [Version 4.3](https://github.com/shawnlaffan/biodiverse/releases/tag/r4.3)  
- [Version 4.2](https://github.com/shawnlaffan/biodiverse/releases/tag/r4.2)  

Docker compressed image size ~ 825.52 MB  
Docker decompressed image size ~ 2.78 GB  
Singularity image size ~ 636.11 MB  

Docker image digest: `sha256:9ce19b4e4da42268841c9cd88619cde63739ed3f9870870a9c4fe76dfb7f4e51`  
Singularity unique ID: `sha256.88a4ac5ac625345a8d3d87f918213b2f57d2109b597971d990ff3165a19463b3`  


### v.1.2.0 (Mar 15, 2023)

- Biodiverse v4.1 `d70b531` (https://github.com/vmikk/biodiverse-docker/releases/tag/v.1.2.0)  

Biodiverse changelog:  
- [Version 4.1](https://github.com/shawnlaffan/biodiverse/releases/tag/r4.1)  

Docker compressed image size ~ 820.3 MB  
Docker decompressed image size ~ 2.77 GB  
Singularity image size ~ 635.50 MB  

Docker image digest: `sha256:4cf5e5800b804303f7c9a90ef201c301ce63a84590c21fde5d8eaff0ba66da44`  
Singularity unique ID: `sha256.031397ea4334a40db8d94af08a3791dff9f9209990f6990be8dc80b925880098`  

### v.1.0.0 (Dec 17, 2022)

- Biodiverse v4 `fe66697` (https://github.com/vmikk/biodiverse-docker/releases/tag/v.1.0.0)  

Biodiverse changelog:  
- [Version 4.0](https://github.com/shawnlaffan/biodiverse/releases/tag/r4.0)  
- [Version r3.99_005](https://github.com/shawnlaffan/biodiverse/releases/tag/r3.99_005)  
- [Version r3.99_004](https://github.com/shawnlaffan/biodiverse/releases/tag/r3.99_004)  

Compressed Size ~ 818.57 MB  
Decompressed Size ~ 2.76 GB  

Docker image digest: `sha256:0f6d5c64938d7bd9cd54943a8f717c15bdc008672f256f2b93bc74e0caa65f1d`  
Singularity unique ID: `sha256.3c9edf889a2339fa3d1d51ddff54a8dd7a9b969e9da4d1a651fc783f41b78326`  

### v.0.0.2 (Jun 13, 2022)

- Biodiverse v.3.99_003 `c06bbad`  
- Remove `perlbrew` as a dependency.  
https://github.com/vmikk/biodiverse-docker/releases/tag/v.0.0.2

### v.0.0.1 (Apr 22, 2022)

- Biodiverse v.3.99_002 `d20a9ce`  
https://github.com/vmikk/biodiverse-docker/releases/tag/v.0.0.1



## `OpenTree` changelog

A Docker image with `OpenTree` and `Python` for phylogenetic tree fetching from Open Tree of Life.  
Docker hub: [https://hub.docker.com/r/vmikk/opentree/tags](https://hub.docker.com/r/vmikk/opentree/tags)  
Singularity Library: [https://cloud.sylabs.io/library/vmiks/gbif/opentree](https://cloud.sylabs.io/library/vmiks/gbif/opentree)  

### v.0.0.2 (Apr 22, 2022)

Compressed Size ~ 360.47 MB  
Decompressed Size ~ 968 MB  

Docker image digest: `sha256:8b1dba9bd654e4fe5768f42e09b3444aa5a190c9bb622692eb19438171a9a85a`  
Singularity unique ID: `sha256.b41dd6c6ec3a2b6da41d594b30e16d035477ef218ebd89c907106024166003b5`  
