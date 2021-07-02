# Ritmic
**Ritmic** is an R package for studying the RegulatIon of Tumor MIcroenvironment Composition.
The pipeline **Ritmic**  is ine 5 steps:
1. (optionnal) - Simulation of a D matrix, corresponding to the complex profile of samples, and that can be decomposed in T\*A, with T the pure profile of each cell type, and A the proportion of each cell type in each sample. 
2. Deconvolution of the D matrix, with the R packages *medepir* and *EDec*.
3. Identification and extraction of the tumor cell type in the data obtained
4. Application of the *PenDA* method for individual expression analysis on the tumor cell type.
5. Link between gene deregulation and micro-environment composition, graphical visualization. 


## Installation
To get the current development version from github:
```
git clone https://github.com/bcm-uga/Ritmic.git
cd Ritmic 
R
```

## Build package

```R
install.packages("devtools")
# mixtools htmltools scales yaml lazyeval plyr rlang ggplot2 gtools caTools KernSmooth penda progress
devtools::load_all(); devtools::document(); devtools::install()
```

## Build vignette 
```
setwd("vignettes")
rmarkdown::render("simulation.Rmd")
rmarkdown::render("penda_analysis.Rmd")
```
