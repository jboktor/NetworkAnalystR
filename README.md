
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NetworkAnalystR

<!-- badges: start -->

<!-- badges: end -->

The goal of NetworkAnalystR is to …

## Installation

You can install the released version of NetworkAnalystR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xia-lab/NetworkAnalystR")
#> Downloading GitHub repo xia-lab/NetworkAnalystR@master
#> Skipping 7 packages ahead of CRAN: S4Vectors, BiocGenerics, AnnotationDbi, DelayedArray, Biobase, IRanges, BiocParallel
#>      checking for file ‘/tmp/Rtmpp4QxUx/remotes1a2f7c99d985/xia-lab-NetworkAnalystR-032c5ad0def50d86f03c86210d8c15a1c9bd6c3b/DESCRIPTION’ ...  ✓  checking for file ‘/tmp/Rtmpp4QxUx/remotes1a2f7c99d985/xia-lab-NetworkAnalystR-032c5ad0def50d86f03c86210d8c15a1c9bd6c3b/DESCRIPTION’
#>   ─  preparing ‘NetworkAnalystR’:
#>    checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
#>   ─  excluding invalid files
#>      Subdirectory 'R' contains invalid file names:
#>      ‘_graph_eda.R’ ‘_graph_noa.R’ ‘_graph_sif.R’ ‘_ls_objects.R’
#>      ‘_parseListInput.R’ ‘_prepareExpressSeeds.R’ ‘_prepareListSeeds.R’
#>      ‘_prepareSigProteinJSON.R’ ‘_query_sqlite.R’
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>   ─  looking to see if a ‘data/datalist’ file should be added
#>   ─  building ‘NetworkAnalystR_0.0.0.9000.tar.gz’
#>      
#> 
#> Installing package into '/usr/local/lib/R/site-library'
#> (as 'lib' is unspecified)
```

## Example Data

Importing example data:

``` r
# Create objects for storing processed data
nSet <- Init.Data()

# Read in the data and fill in the dataSet list 
nSet <- ReadTabExpressData("data/test/estrogen.txt")

# Perform data annotation
nSet <- PerformDataAnnot("hsa", "array", "hgu95av2", "mean")
```

## README.Rmd

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.
