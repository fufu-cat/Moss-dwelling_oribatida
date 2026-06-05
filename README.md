# Moss-dwelling_oribatida

This repository contains the data and R scripts used for the manuscript:

**Multidimensional Beta Diversity Dominated by Replacement Reveals Hidden Complementarity of Moss-Associated Oribatid Mites in Alpine Petrifying Springs**

## Overview

This study investigates moss-associated oribatid mite assemblages inhabiting carbonate-encrusted moss mats in alpine petrifying springs on the south-eastern Tibetan Plateau. We quantified taxonomic beta diversity (TD), functional beta diversity (FD) and phylogenetic beta diversity (PD), partitioned each dimension into replacement and richness-difference components, and examined the ecological drivers of multidimensional beta diversity across five tufa landscape types.

## File Descriptions

- **oribatida.csv**  
  Genus-level abundance matrix of moss-associated oribatid mites across 25 site-level moss-mat samples.

- **traits.csv**  
  Functional trait matrix of oribatid mite genera used for FD analysis.

- **moss.csv**  
  Moss community composition data used to represent moss community gradients.

- **env.csv**  
  Measured ecological predictor variables, including stream–substrate conditions, moss functional status and spatial factors.

- **01 Oribatida_phylo_tree_final_Grafen.tre**  
  Genus-level phylogenetic tree of oribatid mites used for PD analysis.

- **02 Oribatida_functional_tree_NJ.tre**  
  Trait-based functional dendrogram used for FD analysis.

- **4_R Scripts.R**  
  Complete R script used to reproduce the statistical analyses and visualizations presented in the manuscript.

## Usage

Download or clone this repository, open the R script in the `R Scripts` folder, and run:

```r
source("R Scripts/4_R Scripts.R")
```

Please make sure that the working directory is set to the root folder of this repository before running the script.

## Data Availability

All data and scripts required to reproduce the analyses are provided in this repository.

## License

This repository is licensed under the Apache License 2.0. See the `LICENSE` file for details.

## Contact

For questions about the data or scripts, please contact:

**Yuanyuan Zhou**  
Email: 233100170050@gznu.edu.cn

**Zhaohui Zhang**  
Email: zhaozhang9@hotmail.com
