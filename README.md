
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

### 1. Construct predicting models for prognosis
``` r
library(Mime)
load("./Example.cohort.Rdata")
load("./genelist.Rdata")
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                     list_train_vali_Data = list_train_vali_Data,
                     unicox.filter.for.candi = T,
                     unicox_p_cutoff = 0.05,
                     candidate_genes = genelist,
                     mode = 'all',nodesize =5,seed = 5201314 )
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/cindex_dis_all.png)

``` r
cindex_dis_select(res,
                  model="StepCox[forward] + plsRcox",
                  order= names(list_train_vali_Data))
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/cindex_specific_model.png)

``` r
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + plsRcox",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/sur_km.png)

``` r
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                            auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                            auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                            auc_cal_method="KM")
```
``` r
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)

```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/auc1y_dis_all.png)

``` r
roc_vis(all.auc.1y,
        model_name = "StepCox[forward] + plsRcox",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/auc1y_roc.png)

``` r
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + plsRcox",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/auc_specific_model.png)

``` r
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = "StepCox[forward] + plsRcox",type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/meta_rs.png)

``` r
rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('LGG','GBM','Glioma'),
                                        list_input_data = list_train_vali_Data)
HR_com(rs.glioma.lgg.gbm,
       res,
       model_name="StepCox[forward] + plsRcox",
       dataset=names(list_train_vali_Data),
       type = "categorical")
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/hr_comp.png)

``` r
cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="StepCox[forward] + plsRcox",
            dataset=names(list_train_vali_Data))
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/cindex_comp.png)

``` r
auc.glioma.lgg.gbm.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data,AUC_time = 1,
                                            auc_cal_method = 'KM')
auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name="StepCox[forward] + plsRcox",
         dataset=names(list_train_vali_Data))
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/auc_comp_1y.png)

### 2. Construct predicting models for response
``` r
load("./Example.ici.Rdata")
load("./genelist.Rdata")
res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$training,
                                      list_train_vali_Data = list_train_vali_Data,
                                      candidate_genes = genelist,
                                      methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                      seed = 5201314,
                                      cores_for_parallel = 60
)
```
``` r
auc_vis_category_all(res.ici,dataset = c("training","validation"),
                     order= c("training","validation"))
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/ICI_response_auc_all.png)

``` r
plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/ICI_response_roc_all.png)

### 3. Core feature selection

``` r
load("./Example.cohort.Rdata")
load("./genelist.Rdata")
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Dataset1,
                                            candidate_genes = genelist,
                                            mode = "all",nodesize =5,seed = 5201314 )
core_feature_select(res.feature.all)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/core_feature_intersect.png)

``` r
core_feature_rank(res.feature.all, top=20)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/core_feature_intersect_rank.png)

``` r
dataset_col<-c("#3182BDFF","#E6550DFF")
corplot <- list()
for (i in c(1:2)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                               color = dataset_col[i],
                               feature1="PSEN2",
                               feature2="WNT5B",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=2)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/gene_cor.png)

``` r
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-core_feature_sur("PSEN2", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        #color=c("blue","green"),
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
```
![Screenshot](https://github.com/l-magnificence/Mime/blob/main/fig/gene_km.png)

