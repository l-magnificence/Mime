
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mime

<!-- badges: start -->
<!-- badges: end -->

    The Mime package provides a user-friendly solution for constructing 
    machine learning-based integration models from transcriptomic data. With the 
    widespread use of high-throughput sequencing technologies, understanding biology 
    and cancer heterogeneity has been revolutionized. Mime streamlines the process of 
    developing predictive models with high accuracy, leveraging complex datasets to 
    identify critical genes associated with disease progression, patient outcomes, 
    and therapeutic response. It offers four main applications (i) establishing 
    prognosis models using 10 machine learning algorithms, (ii) building binary 
    response models with 7 machine learning algorithms, (iii) conducting core feature 
    selection related to prognosis using 8 machine learning methods, and (iv) visualizing 
    the performance of each model.

## Installation

You can install the development version of Mime like so:

``` r
if (!requireNamespace("Mime", quietly = TRUE))
  devtools::install_github("l-magnificence/Mime")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(Mime)
## basic example code
```

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
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
