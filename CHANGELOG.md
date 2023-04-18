# PhyloNext: Changelog

This project tries to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).  
For version numbering, we use the following convention: `MAJOR.MINOR.PATCH`.  
Each element increases numerically (e.g., `1.9.0` -> `1.10.0` -> `1.11.0`).  


## v1.2.0 - [Apr 11, 2023]

- `Added`: spatially-constrained randomizations;  
- `Added`: spatial filtering using user-supplied polygons;  
- `Added`: CANAPE output in tabular format from Biodiverse;  
- `Added`: `leaflet_canapesuper` parameter (for 3-class or 4-class endemism types in CANAPE);  
- `Added`: geographic coordinates of H3 cells to the output table;  
- `Added`: `leaflet_mapshots` script;  
- `Changes`: previously, the spatial unit was categorical (H3 grid cell); to spatially constrain records, now we use geographic coordinates of grid cells;  
- `Changes`: record counts, derived dataset, and additional filtering now re-use the same code base;  
- `Changes`: screen output is less verbose now;  
- `Changes`: PRNG randomization seeds are explicit now;  
- `Fixed`: Number of Biodiverse randomizations per thread;  
- `Fixed`: CANAPE index requirements;  
- `Fixed`: antimeridian issue in static visualization;  
- `Fixed`: record counts and derived dataset take polygon-filtering into account;  

## v.1.0.0 - [Feb 1, 2023]

First stable version.  

- `Added`: derived datasets;  
- `Added`: basis-of-record filtering;  
- `Added`: GeoPackage output;  
- `Added`: Newick tree output;  
- `Added`: `maxage` and `phyloonly` parameters for phylogenetic tree fetching from the Open Tree of Life;  
- `Added`: `leaflet_sescolor` parameter to adjust Z-score coloring scheme;  
- `Added`: CI tests;  
- `Fixed`: urban areas shapefile;  

## v.0.0.2 - [Dec 22, 2022]

Code refactoring targeting for the execution in a cloud environment (in particular, MS Azure Batch).  
- `Added`: Leaflet-based interactive visualization;  
- `Added`: redundancy index;  
- `Added`: removal of extinct species;  
- `Added`: removal of human records (genus _Homo_);  
- `Fixed`: correct usage of `val` and `path` for input file staging;  
- `Changes`: `class` is a reserved word in Nextflow; now we are using `classis` parameter to refer taxonomic rank of Class;  

## v.0.0.1 - [Jun 1, 2022]

The initial release of PhyloNext.  
First commit - Feb 09, 2022.  

