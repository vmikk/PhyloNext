# PhyloNext: Changelog

This project tries to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).  
For version numbering, we use the following convention: `MAJOR.MINOR.PATCH`.  
Each element increases numerically (e.g., `1.9.0` -> `1.10.0` -> `1.11.0`).  


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

