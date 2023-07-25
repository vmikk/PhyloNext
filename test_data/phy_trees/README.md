# Phylogenetic trees

## Acacia

`Acacia_Mishler,2014_Latin-labels.nwk`

Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H., Laffan, S. W., & Miller, J. T. (2014). Phylogenetic measures of biodiversity and neo- and paleo-endemism in Australian Acacia. Nature Communications, 5(1), 4473. https://doi.org/10.1038/ncomms5473


`Acacia_OpenTree_bladj_Latin-labels.nwk`
`Acacia_OpenTree_bladj_OTT-labels.tre`

Trees deriverd from the Open Tree of Life and dated with `bladj`. Node age esimates are from:  

Miller, J. T., Murphy, D. J., Ho, S. Y. W., Cantrill, D. J., & Seigler, D. (2013). Comparative dating of Acacia: Combining fossils and multiple phylogenies to infer ages of clades with poor fossil records. Australian Journal of Botany, 61(6), 436–445. https://doi.org/10.1071/BT13149


## Mammals, Amphibians, and Birds

The phylogenetic trees used to benchmark the EDGE2 metric were based on published phylogenies, but were augmented with species not present in the original phylogenies. 
To factor in the uncertainty of imputation, the original authors utilized 1000 imputed trees.  
Here, only one "median" tree per taxon has been uploaded. 
This "median" tree is identified as the one with the lowest squared distance to all other trees.


- `Mammals_Gumbs,2022_Latin-labels.nwk`  
- `Amphibians_Gumbs,2022_Latin-labels.nwk`  
- `Birds_Gumbs,2022_Latin-labels.nwk`  


Reference:  
Gumbs R, Gray CL, Böhm M, Burfield IJ, Couchman OR, Faith DP, Forest F, Hoffmann M, Isaac NJB, Jetz W, Mace GM, Mooers AO, Safi K, Oenone Scott O, Steel M, Tucker CM, Pearse WD, Owen NR, Rosindell J (2023) The EDGE2 protocol: Advancing the prioritisation of Evolutionarily Distinct and Globally Endangered species for practical conservation action. **PLOS Biology** 21(2): e3001991. [DOI:10.1371/journal.pbio.3001991](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001991)  


## Open Tree of Life

The phylogenetic trees below were obtained from the Open Tree of Life on February 26, 2022. 
Use these trees with caution, as they represent the synthetic approximately dated trees. 
These trees use OTT IDs for their tip labels, so if you wish to use them in PhyloNext, 
you should specify the parameter `--phylabels OTT`.

Included in this set are trees for the following groups: 
Amphibians, Birds, Australian Fabaceae, Gymnosperms, and Mammals:  

- `Amphibians_OpenTree_OTT-labels.nwk`  
- `Birds_OpenTree_OTT-labels.nwk`  
- `Fabaceae_Australia_OpenTree_OTT-labels.nwk`  
- `Gymnosperms_OpenTree_OTT-labels.nwk`  
- `Mammals_OpenTree_OTT-labels.nwk`  
- `Reptiles_OpenTree_OTT-labels.nwk`  

You can find additional information on how these trees were procured at the following link:  
https://github.com/McTavishLab/GBIF-Biodiverse-OpenTree  

