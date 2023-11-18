
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mime

<!-- badges: start -->
<!-- badges: end -->

The `Mime` package provides a user-friendly solution for constructing
machine learning-based integration models from transcriptomic data.

With the widespread use of high-throughput sequencing technologies,
understanding biology and cancer heterogeneity has been revolutionized.
Mime streamlines the process of developing predictive models with high
accuracy, leveraging complex datasets to identify critical genes
associated with disease progression, patient outcomes, and therapeutic
response.

It offers four main applications:

1.  Establishing prognosis models using 10 machine learning algorithms.
2.  Building binary response models with 7 machine learning algorithms.
3.  Conducting core feature selection related to prognosis using 8
    machine learning methods.
4.  Visualizing the performance of each model.
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/graph.jpg)
## Install

You can install the development version of Mime like so:

``` r
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' )
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE))
  devtools::install_github("souravc83/fastAdaboost")

if (!requireNamespace("Mime", quietly = TRUE))
  devtools::install_github("l-magnificence/Mime")
  
library(Mime)
```

## Quick Start

This is a basic example which shows you how to use `Mime`:

Mime need multiple cohorts containing transcriptional sequencing data with information of survival or clinical response to therapy as well as a gene set as inputs. You can dowmload the example data in `External data` file to check the format of input.
``` r
load("./Example.cohort.Rdata")
list_train_vali_Data[["Dataset1"]][1:5,1:5]
#>               ID    OS.time OS   MT-CO1   MT-CO3
#>  TCGA.DH.A66B.01 1281.65322  0 13.77340 13.67931
#>  TCGA.HT.7607.01   96.19915  1 14.96535 14.31857
#>  TCGA.DB.A64Q.01  182.37755  0 13.90659 13.65321
#>  TCGA.DU.8167.01  471.97707  0 14.90695 14.59776
#>  TCGA.HT.7610.01 1709.53901  0 15.22784 14.62756
```
The first column `ID` is the sample ID, second to third column `OS.time` and `OS` are the survival time and status of patients, other columns are the gene expression level. `Dataset1` is the training dataset, while other Datasets as validation.
``` r
load("./Example.ici.Rdata")
list_train_vali_Data[["training"]][1:5,1:5]
#>                                   ID Var      FTH1   EEF1A1      ACTB
#>                      SAMf2ce197162ce   N 10.114846 4.817746 11.230180
#>                           ERR2208915   Y  2.044180 5.038854  3.977902
#> G138701_RCCBMS-00141-T_v1_RNA_OnPrem   Y  5.406008 5.341635  5.366668
#>                      SAMe41b1e773582   N  9.215794 4.707360 11.412721
#>                      SAM5ffd7e4cd794   N  9.003710 3.908884 10.440559

```
The first column `ID` is the sample ID, second column `Var` is the theraputic response of patients (`N`: No response; `Y`: response), other columns are the gene expression level. `training` is the training dataset, while other Datasets as validation.
``` r
load("./genelist.Rdata")
#> [1] "MYC"    "CTNNB1" "JAG2"   "NOTCH1" "DLL1"   "AXIN2"  "PSEN2"  "FZD1"   "NOTCH4" "LEF1"   "AXIN1"  "NKD1"   "WNT5B" 
#>[14] "CUL1"   "JAG1"   "MAML1"  "KAT2A"  "GNAI1"  "WNT6"   "PTCH1"  "NCOR2"  "DKK4"   "HDAC2"  "DKK1"   "TCF7"   "WNT1"  
#>[27] "NUMB"   "ADAM17" "DVL2"   "PPARD"  "NCSTN"  "HDAC5"  "CCND2"  "FRAT1"  "CSNK1E" "RBPJ"   "FZD8"   "TP53"   "SKP2"  
#>[40] "HEY2"   "HEY1"   "HDAC11"

```
This gene set is associated with Wnt/β-catenin signalling from MSigDB.

We recommend training dataset with more than 100 samples and gene set with more than 50 genes.

### 1. Construct prognostic models
``` r
library(Mime)
load("./Example.cohort.Rdata")
load("./genelist.Rdata")
res= ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                     list_train_vali_Data = list_train_vali_Data,
                     unicox.filter.for.candi = T,
                     unicox_p_cutoff = 0.05,
                     candidate_genes = genelist,
                     mode = 'all',nodesize =5,seed = 5201314 )
```
``` r
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/cindex_dis_all.png)

