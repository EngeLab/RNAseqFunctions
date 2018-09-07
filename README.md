RNAseqFunctions
================

Build:
<a href="https://travis-ci.org/EngeLab/sRNAseqFunctions"><img src="https://travis-ci.org/EngeLab/RNAseqFunctions.svg?branch=master"></a>

Test coverage:
[![codecov](https://codecov.io/gh/EngeLab/RNAseqFunctions/branch/master/graph/badge.svg)](https://codecov.io/gh/EngeLab/RNAseqFunctions)

### Description

A package of RNAseq QC and analysis tools for the Enge lab.

### Installation

``` r
if(!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools")
}

devtools::install_github("EngeLab/RNAseqFunctions")
```

### Vignettes

The vignette for loading and filtering counts data can be viewed in R
with the following command:

``` r
vignette("loadAndFilterFuns", package = "RNAseqFunctions")
```

The vignette for running t-SNE and plotting can be viewed in R with the
following command:

``` r
vignette("tsneAndPlotting", package = "RNAseqFunctions")
```
