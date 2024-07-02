#' Calculating the Area Under the Curve of the previous signatures for immunotherapy response.
#'
#' @param list_train_vali_Data A list containing the training data and the other validation data. All the validation data have the same data form as the training data.
#' @param seed The seed you can set as any positive integer, for example, 5201314
#' @param train_data The training data with the 'ID' and 'Var' as the first two columns. Starting in the third column are the variables used to construct the model. 'Var' is the target predictor variable for constructing the model. 'Var' contains only Y or N.
#' @param cores_for_parallel The cores you can choose for parallel operation. The default is 12. The bigger the better if the configuration allows it.
#'
#' @return A list of the AUC results of each signature we collected in each data. 
#' @export
#'
#' @examples
cal_auc_previous_sig <- function(list_train_vali_Data, # lsit of the cohort, 第一列为ID,第二列为Var, Y/N, 从第三列开始是基因名称
                                 seed = 5201314,
                                 train_data,
                                 cores_for_parallel = 12) {
  library(dplyr)
  
  cat("Data processing")

  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, as.numeric)
    rownames(x) <- x$ID
    return(x)
  })

  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
    return(x)
  })


  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x <- x[!is.na(x$Var) & !is.na(x$Var), ]
    return(x)
  })

  # use the mean replace the NA
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })


    return(x)
  })



  ###########################################################################################################
  #### ================== Compare sig with other Signatures in  datasets====================####
  ###########################################################################################################

  ## ==loading signatures==##
  # IMS.Sig
  IMS.Sig <- c("CCL8", "VCAN", "CCL2", "BCAT1", "ISG15", "CD163", "AXL", "CCL13", "COL6A3", "SIGLEC1", "PDGFRB", "IL10", "STC1", "ADAM12", "OLFML2B", "FAP", "TWIST2", "INHBA")
  # TRS.Sig
  TRS.Sig <- c("CTLA4", "CXCR6", "LYST", "CD38", "GBP2", "HLA-DRB5")
  # NLRP3.Sig
  NLRP3.Sig <- c("ARRDC1-AS1", "CARD8", "GSDMD", "ATAT1", "CD36", "CPTP", "DHX33", "EIF2AK2", "GBP5", "NLRC3", "PYDC2", "SIRT2", "TLR4", "TLR6", "USP50", "APP", "CASP1", "HSP90AB1", "MEFV", "NFKB1", "NFKB2", "NLRP3", "P2RX7", "PANX1", "PSTPIP1", "PYCARD", "RELA", "SUGT1", "TXN", "TXNIP")
  # IFNG.Sig
  IFNG.Sig <- c("IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA")
  # CRMA.Sig
  CRMA.Sig <- c("CSAG1", "CSAG2", "CSAG3", "MAGEA2", "MAGEA2B", "MAGEA3", "MAGEA6", "MAGEA12")
  # T.cell.inflamed.Sig
  T.cell.inflamed.Sig <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
  # Cytotoxic.Sig
  Cytotoxic.Sig <- c("GZMA", "PRF1")
  # PDL1.Sig
  PDL1.Sig <- c("PDL1", "PDCD1")
  # LRRC15.CAF.Sig
  LRRC15.CAF.Sig <- c("MMP11", "COL11A1", "C1QTNF3", "CTHRC1", "COL12A1", "COL10A1", "COL5A2", "GJB2", "THBS2", "AEBP1", "MFAP2", "LRRC15", "PLAU", "ITGA11") # Alias for 'PRLHR' are:'GR3','GPR10','PrRPR'
  ImmmunCells.Sig <- Mime1::ImmmunCells.Sig
  TcellExc.Sig <- Mime1::TcellExc.Sig
  exc.sig <- Mime1::exc.sig

  ls_sig <- list(
    "ImmmunCells.Sig" = ImmmunCells.Sig,
    "TcellExc.Sig" = TcellExc.Sig,
    "IFNG.Sig" = IFNG.Sig,
    "T.cell.inflamed.Sig" = T.cell.inflamed.Sig,
    "Cytotoxic.Sig" = Cytotoxic.Sig,
    "NLRP3.Sig" = NLRP3.Sig,
    "LRRC15.CAF.Sig" = LRRC15.CAF.Sig,
    "PDL1.Sig" = PDL1.Sig,
    "CRMA.Sig" = CRMA.Sig,
    "TRS.Sig" = TRS.Sig
  )

  auc.list.immunotherapy <- list()


  ##### =============================1 immunotherapy response signatures==============================================#####

  #### 1.1 IFNG.Sig

  # IFNG.Sig use average gene expression as prediction socres
  # Reference: Ayers M, Lunceford J, Nebozhyn M, Murphy E, Loboda A, Kaufman DR, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–40.
  # Available from: Available from: https://www.jci.org/articles/view/91190

  cat("1.1 IFNG.Sig")

  sig <- ls_sig[["IFNG.Sig"]]
  sig

  sig <- gsub("-", ".", sig)




  AUC_INFG.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)]


    roc <- ROCit::rocit(
      score = rowMeans(new[, 2:7]),
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")

  AUC_INFG.Sig

  auc.list.immunotherapy[["IFNG.Sig"]] <- AUC_INFG.Sig

  #### 1.2 T.cell.inflamed.Sig

  # T.cell.inflamed.Sig use average gene expression as prediction socres
  # Reference: Ayers M, Lunceford J, Nebozhyn M, Murphy E, Loboda A, Kaufman DR, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–40.
  # Available from: Available from: https://www.jci.org/articles/view/91190

  cat("1.2 T.cell.inflamed.Sig")


  sig <- ls_sig[["T.cell.inflamed.Sig"]]
  sig
  sig <- gsub("-", ".", sig)


  AUC_T.cell.inflamed.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)]


    roc <- ROCit::rocit(
      score = rowMeans(new[, 2:7]),
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")

  AUC_T.cell.inflamed.Sig

  auc.list.immunotherapy[["T.cell.inflamed.Sig"]] <- AUC_T.cell.inflamed.Sig


  #### 1.3 PDL1.Sig

  # PDL1.Sig was the expression of PD-L1
  # Reference: Topalian SL, Hodi FS, Brahmer JR, Gettinger SN, Smith DC, McDermott DF, et al. Safety, Activity, and Immune Correlates of Anti–PD-1 Antibody in Cancer. N Engl J Med. 2012;366:2443–54
  # Available from: doi: 10.1056/NEJMoa1200690
  cat("1.3 PDL1.Sig")

  sig <- ls_sig[["PDL1.Sig"]]
  sig
  sig <- gsub("-", ".", sig)


  AUC_PDL1.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)]

    roc <- ROCit::rocit(
      score = new$PDCD1,
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")

  AUC_PDL1.Sig
  auc.list.immunotherapy[["PDL1.Sig"]] <- AUC_PDL1.Sig



  #### 1.4 LRRC15.CAF.Sig

  # LRRC15.CAF.Sig use eigenWeightedMean method in sparrow package to calculate the prediction scores
  # Reference: Dominguez CX, Müller S, Keerthivasan S, Koeppen H, Hung J, Gierke S, et al. Single-Cell RNA Sequencing Reveals Stromal Evolution into LRRC15 + Myofibroblasts as a Determinant of Patient Response to Cancer Immunotherapy. Cancer Discov [Internet]. 2020;10:232–53.
  # Available from: http://cancerdiscovery.aacrjournals.org/lookup/doi/10.1158/2159-8290.CD-19-0644

  library(sparrow)
  sig <- ls_sig[["LRRC15.CAF.Sig"]]
  sig
  sig <- gsub("-", ".", sig)
  cat("1.4 LRRC15.CAF.Sig")


  AUC_LRRC15.CAF.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)]

    expr <- as.data.frame(t(new[, -1]))

    scores <- sparrow::eigenWeightedMean(expr)$score

    roc <- ROCit::rocit(
      score = scores,
      class = new$Var,
      negref = "Y"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")

  AUC_LRRC15.CAF.Sig

  auc.list.immunotherapy[["LRRC15.CAF.Sig"]] <- AUC_LRRC15.CAF.Sig


  #### 1.5 NLRP3.Sig
  # NLRP3.Sig was calculated using ssGSEA
  # Reference: Ju M, Bi J, Wei Q, et al. Pan-cancer analysis of NLRP3 inflammasome with potential implications in prognosis and immunotherapy in human cancer. Brief Bioinform 2020; 00: 1–16.
  # Available from: DOI: 10.1093/bib/bbaa345
  cat("1.5 NLRP3.Sig")

  library(GSVA)
  library(GSEABase)
  sig <- NLRP3.Sig
  sig <- gsub("-", ".", sig)

  gmt <- Mime1::NLRP3.Sig.gmt
  all(duplicated(names(gmt))) # no duplicated ids

  AUC_NLRP3.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g
    new$Var <- factor(new$Var, levels = c("Y", "N"))

    expr <- as.matrix(t(new[, -c(1:2)]))

    gsva <- gsva(expr, gmt, method = "ssgsea", parallel.sz = 1)

    score <- as.numeric(gsva)

    roc <- ROCit::rocit(
      score = score,
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")


  AUC_NLRP3.Sig

  auc.list.immunotherapy[["NLRP3.Sig"]] <- AUC_NLRP3.Sig


  #### 1.6 Cytotoxic.Sig
  cat("1.6 Cytotoxic.Sig")

  # Cytotoxic.Sig use use geometric mean of gene expression as prediction socres
  # Reference:Rooney MS, Shukla SA, Wu CJ, Getz G, Hacohen N. Molecular and genetic properties of tumors associated with local immune cytolytic activity. Cell [Internet]. Elsevier Inc.; 2015;160:48–61.
  # Available from: http://dx.doi.org/10.1016/j.cell.2014.12.033

  sig <- ls_sig[["Cytotoxic.Sig"]]
  sig

  AUC_Cytotoxic.Sig <- lapply(list_train_vali_Data, function(g) {
    Cytotoxic_combat <- g[, 3:ncol(g)]
    # log2 changes should be removed for the following geometric mean calculation
    Cytotoxic_combat <- Cytotoxic_combat^2
    Cytotoxic_combat <- cbind(g[, 1:2], Cytotoxic_combat)
    rownames(Cytotoxic_combat) <- g$ID
    Cytotoxic_combat <- Cytotoxic_combat[!is.na(Cytotoxic_combat$Var), ] ## remove patient with unknown response status

    Cytotoxic_data <- Cytotoxic_combat
    new <- Cytotoxic_data[, colnames(Cytotoxic_data) %in% c("Var", sig)]


    roc <- ROCit::rocit(
      score = compositions::geometricmeanRow(new[, 2:ncol(new)]),
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")

  AUC_Cytotoxic.Sig

  auc.list.immunotherapy[["Cytotoxic.Sig"]] <- AUC_Cytotoxic.Sig


  #### 1.7 ImmuneCells.Sig

  ## ImmuneCells.Sig Model was constructed using "cancerclass" according to the original article
  # Reference:Xiong D, Wang Y, You M. A gene expression signature of TREM2hi macrophages and γδ T cells predicts immunotherapy response. Nat Commun. 2020;11:1–12.
  # Available from: http://dx.doi.org/10.1038/s41467-020-18546-x
  cat("1.7 ImmuneCells.Sig")

  sig <- ls_sig[["ImmmunCells.Sig"]]

  model.Dev <- function(training, method, sig) {
    training <- training[, colnames(training) %in% c("Var", sig)]
    # 7 models adpoted in this study as followings:
    #' nb': navie bayes
    #' svmRadialWeights': Support Vector Machines with Class Weights
    #' rf': random forest
    #' kknn': k-Nearest Neighbors
    #' adaboost':AdaBoost Classification Trees
    #' LogitBoost':Boosted Logistic Regressions
    #' cancerclass': cancerclass

    # Grid search for parameter tuning
    Grid <- list(
      nb = expand.grid(fL = c(0, 0.5, 1, 1.5, 2.0), usekernel = TRUE, adjust = c(0.5, 0.75, 1, 1.25, 1.5)),
      svmRadialWeights = expand.grid(sigma = c(0.0005, 0.001, 0.005, 0.01, 0.05), C = c(1, 3, 5, 10, 20), Weight = c(0.1, 0.5, 1, 2, 3, 5, 10)),
      rf = expand.grid(mtry = c(2, 42, 83, 124, 165, 205, 246, 287, 328, 369)),
      kknn = expand.grid(kmax = c(5, 7, 9, 11, 13), distance = 2, kernel = "optimal"),
      adaboost = expand.grid(nIter = c(50, 100, 150, 200, 250), method = c("Adaboost.M1", "Real adaboost")),
      LogitBoost = expand.grid(nIter = c(11, 21, 31, 41, 51, 61, 71, 81, 91, 101))
    )
    TuneLength <- list(
      nb = nrow(Grid[["nb"]]),
      svmRadialWeights = nrow(Grid[["svmRadialWeights"]]),
      rf = nrow(Grid[["rf"]]),
      kknn = nrow(Grid[["kknn"]]),
      adaboost = nrow(Grid[["adaboost"]]),
      LogitBoost = nrow(Grid[["LogitBoost"]])
    )
    ls_model <- lapply(method, function(m) {
      if (m == "cancerclass") { # cancerclass is not avaliable in caret
        pData <- data.frame(class = training$Var, sample = rownames(training), row.names = rownames(training))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(training[, -1])
        Sig.Exp.train <- ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)
        predictor <- fit(Sig.Exp.train, method = "welch.test")
        model.tune <- predictor
      } else {
        f <- 5 # f folds resampling
        r <- 10 # r repeats
        n <- f * r

        # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
        seeds <- vector(mode = "list", length = n + 1)
        # the number of tuning parameter
        for (i in 1:n) seeds[[i]] <- sample.int(n = 1000, TuneLength[[m]])

        # for the last model
        seeds[[n + 1]] <- sample.int(1000, 1)


        ctrl <- trainControl(
          method = "repeatedcv",
          number = f, ## 5-folds cv
          summaryFunction = twoClassSummary, # Use AUC to pick the best model
          classProbs = TRUE,
          repeats = r, ## 10-repeats cv,
          seeds = seeds
        )



        model.tune <- train(Var ~ .,
          data = training,
          method = m,
          metric = "ROC",
          trControl = ctrl,
          tuneGrid = Grid[[m]]
        )
      }
      print(m)
      return(model.tune)
    })


    names(ls_model) <- method

    return(ls_model)
  }
  if (T) {
    library(stringr)
    library(gridExtra)
    library(future)
    library(sva)
    library(e1071)
    library(pROC)
    library(ROCit)
    library(caret)
    library(doParallel)
    library(cancerclass)
    library(dplyr)
  }

  res <- model.Dev(
    training = train_data,
    method = "cancerclass",
    sig = sig
  )



  AUC_ImmmunCells.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)]
    new$Var <- factor(new$Var, levels = c("Y", "N"))
    pData <- data.frame(class = new$Var, sample = rownames(new), row.names = rownames(new))
    phenoData <- new("AnnotatedDataFrame", data = pData)
    Sig.Exp <- t(new[, -1])
    Sig.Exp.test <- ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)
    prediction <- predict(res$cancerclass, Sig.Exp.test, positive = "N", ngenes = nrow(Sig.Exp), dist = "cor")
    roc <- ROCit::rocit(
      score = as.numeric(prediction@prediction[, "z"]),
      class = new$Var,
      negref = "Y"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")
  AUC_ImmmunCells.Sig

  auc.list.immunotherapy[["ImmmunCells.Sig"]] <- AUC_ImmmunCells.Sig


  #### 1.8 TcellExc.Sig
  cat("1.8 TcellExc.Sig")

  # TcellExc.Sig use overall expression as prediction scores
  # Reference: Jerby-Arnon L, Shah P, Cuoco MS, Rodman C, Su M-J, Melms JC, et al. A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade. Cell.2018;175:984-997.e24.
  # Available from: https://linkinghub.elsevier.com/retrieve/pii/S0092867418311784

  # source("/export3/zhangw/Project_Cross/Project_Mime/data/sig/IMPRES/ImmRes_source.R") ## 'ImmRes_OE.R' was downloaded from https://github.com/livnatje/ImmuneResistance
  source(system.file("extdata", "ImmRes_source.R", package = "Mime1"))

  library(caret)
  # library(recipes,lib.loc = "/export/bioinfo-team/home/liuhw/R/x86_64-pc-linux-gnu-library/4.1")
  library(recipes)

  sig <- ls_sig[["TcellExc.Sig"]]
  sig <- gsub("-", ".", sig)
  gene.sign <- list(TcellExc.Sig = gsub("-", ".", TcellExc.Sig))



  prd <- function(df) {
    r <- list(tpm = t(df[, -1]), genes = colnames(df[, -1]))
    OE <- get.OE.bulk(
      r = r,
      gene.sign = gene.sign
    )

    roc <- ROCit::rocit(
      score = OE[, "TcellExc.Sig"],
      class = df$Var,
      negref = "Y"
    )
    auc <- roc$AUC
    print(auc)
  }


  getTcellExc.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)] #

    return(prd(new))
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")

  AUC_TcellExc.Sig <- getTcellExc.Sig
  AUC_TcellExc.Sig

  auc.list.immunotherapy[["TcellExc.Sig"]] <- AUC_TcellExc.Sig



  #### 1.9 CRMA.Sig

  # CRMA.Sig use geometric mean of gene expression as prediction socres
  # Shukla SA, Bachireddy P, Schilling B, Galonska C, Zhan Q, Bango C, et al. Cancer-Germline Antigen Expression Discriminates Clinical Outcome to CTLA-4 Blockade. Cell [Internet]. Elsevier Inc.; 2018;173:624-633.e8.
  # Available from: https://doi.org/10.1016/j.cell.2018.03.026

  cat("1.9 CRMA.Sig")

  sig <- ls_sig[["CRMA.Sig"]]
  sig
  sig <- gsub("-", ".", sig)

  # log2 changes should be removed for the following geometric mean calculation

  CRMA_combat <- lapply(list_train_vali_Data, function(x) {
    x[, 3:ncol(x)] <- 2**x[, 3:ncol(x)] - 1
    return(x)
  })


  res <- list()
  p <- list()
  AUC <- list()


  getCRMA.Sig <- lapply(CRMA_combat, function(g) {
    new <- g[, colnames(g) %in% c("Var", sig)] #

    roc <- ROCit::rocit(
      score = compositions::geometricmeanRow(new[, 2:ncol(new)]),
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(CRMA_combat)) %>%
    `colnames<-`("AUC")


  AUC_CRMA.Sig <- getCRMA.Sig
  AUC_CRMA.Sig


  auc.list.immunotherapy[["CRMA.Sig"]] <- AUC_CRMA.Sig



  #### 1.10 IMPRES.Sig

  # IMPRES.Sig uses the comparision scores between 15 gene pairs as prediction socres.
  # Reference:Auslander N, Zhang G, Lee JS, Frederick DT, Miao B, Moll T, et al. Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma. Nat Med [Internet]. Springer US; 2018;24:1545–9.
  # Available from: http://dx.doi.org/10.1038/s41591-018-0157-9
  cat("1.10 IMPRES.Sig")


  gp <- Mime1::gp_IMPRES

  gp <- as.data.frame(gp)
  gp1 <- gp
  gp2 <- gp

  for (i in c("A", "B")) {
    gp2[, i] <- str_replace(gp2[, i], "C10orf54", "VSIR") # C10orf54 and VSIR are interchangable
  }
  g1 <- unique(c(gp1$A, gp1$B))
  g2 <- unique(c(gp2$A, gp2$B))

  IMPRES <- function(pt) {
    pt <- as.data.frame(t(pt))
    IMPRES.scores <- 0
    for (kk in 1:nrow(g)) {
      A <- g[kk, "A"]
      B <- g[kk, "B"]
      try(if (pt[, A] > pt[, B]) {
        IMPRES.scores <- IMPRES.scores + 1
      })
    }
    return(IMPRES.scores)
  }


  ## number of avalible IMPRES.Sig genes for each cohort


  list_impres <- list_train_vali_Data
  for (i in names(list_impres)) {
    io.g <- colnames(list_impres[[i]])
    if ("C10orf54" %in% io.g) {
      print(paste(i, sum(g1 %in% io.g), "C10orf54"))
      list_impres[[i]]$VSIR <- list_impres[[i]]$C10orf54 # add 'VISR' in cohort with "C10orf54"
    } else if ("VSIR" %in% io.g) {
      print(paste(i, sum(g2 %in% io.g), "VSIR"))
    } else {
      print(paste(i, sum(g1 %in% io.g)))
    }
  }

  ls_IMPRES <- list()
  IMPRES_ios <- list_impres


  for (i in names(IMPRES_ios)) {
    g3 <- g2[g2 %in% colnames(IMPRES_ios[[i]])]
    io <- IMPRES_ios[[i]][, g3]
    g <- as.data.frame(gp2)
    g <- g[g$A %in% g3 & g$B %in% g3, ]
    s <- apply(io, 1, IMPRES)
    s <- data.frame(IMPRES = s)
    s$ID <- IMPRES_ios[[i]]$ID
    ls_IMPRES[[i]] <- s
  }





  getIMPRES.Sig <- function() {
    IMPRES_scores <- ls_IMPRES[grp[c(3, 7)]] %>% do.call(rbind, .) # SKCM cohort: Hugo 2016 + Van Allen 2015


    IMPRES_response <- data[data$ID %in% IMPRES_scores$ID, "response"]
    IMPRES_response <- ifelse(IMPRES_response == 0, "NR", "R") %>% factor(., levels = c("R", "NR"))

    roc <- ROCit::rocit(IMPRES_scores$IMPRES, IMPRES_response, negref = "NR")

    auc <- roc$AUC
    return(auc)
  }


  getIMPRES.Sig <- lapply(names(IMPRES_ios), function(g) {
    IMPRES_scores <- ls_IMPRES[[g]]
    g <- IMPRES_ios[[g]]




    new <- g[, colnames(g) %in% c("Var", sig)] #

    IMPRES_response <- new[, "Var"]

    roc <- ROCit::rocit(IMPRES_scores$IMPRES, IMPRES_response, negref = "N")

    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")



  AUC_IMPRES.Sig <- getIMPRES.Sig
  AUC_IMPRES.Sig


  auc.list.immunotherapy[["IMPRES.Sig"]] <- AUC_IMPRES.Sig


  #### 1.11 IPRES.Sig

  # IPRES.Sig used mean z-scores of GSVA scores calculated from 73 IPRES datasets.
  # Reference:Hugo W, Zaretsky JM, Sun L, Johnson DB, Ribas A, Lo RS, et al. Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma Article Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell [Internet]. Elsevier Inc.; 2016;1–10.
  # Available from: http://dx.doi.org/10.1016/j.cell.2016.02.065

  # Original article does not provide the get sets. We manually collected datasets from supplemtary of original article and http://www.gsea-msigdb.org
  cat("1.11 IPRES.Sig")

  ### === GSVA ===###
  library(GSVA)
  library(GSEABase)

  gmt <- Mime1::IPRES_gmt
  all(duplicated(names(gmt))) # no duplicated ids


  getIPRES.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g

    expr <- as.matrix(t(new[, -c(1:2)]))

    gsva <- gsva(expr, gmt, method = "gsva", parallel.sz = cores_for_parallel)

    score <- scale(t(gsva))

    score <- rowMeans(score)

    roc <- ROCit::rocit(
      score = score,
      class = new$Var,
      negref = "Y"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")



  AUC_IPRES.Sig <- getIPRES.Sig
  AUC_IPRES.Sig

  auc.list.immunotherapy[["IPRES.Sig"]] <- AUC_IPRES.Sig



  #### 1.12 TRS.Sig
  # TRS.Sig was calculated using GSVA
  # Reference: Yan M, Hu J, Ping Y, Xu L, Liao G, Jiang Z, Pang B, Sun S, Zhang Y, Xiao Y, Li X. Single-Cell Transcriptomic Analysis Reveals a Tumor-Reactive T Cell Signature Associated With Clinical Outcome and Immunotherapy Response In Melanoma. Front Immunol. 2021 Nov 5;12:758288.
  # Available from: doi: 10.3389/fimmu.2021.758288.
  library(GSVA)
  library(GSEABase)
  sig <- TRS.Sig
  cat("1.12 TRS.Sig")

  # gmt <- Mime1::TRS.Sig.gmt
  all(duplicated(names(gmt))) # no duplicated ids

  getTRS.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g
    colnames(new) <- gsub("\\.", "-", colnames(new))
    expr <- as.matrix(t(new[, -c(1:2)]))
    gsva <- gsva(expr, list(TRS.Sig=TRS.Sig), method = "gsva", parallel.sz = cores_for_parallel)
    score <- as.numeric(gsva)
    roc <- ROCit::rocit(
      score = score,
      class = new$Var,
      negref = "N"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")


  AUC_TRS.Sig <- getTRS.Sig

  AUC_TRS.Sig

  auc.list.immunotherapy[["TRS.Sig"]] <- AUC_TRS.Sig


  #### 1.13 IMS.Sig
  # IMS.Sig was calculated as ratio of (IMS scores / IFN-γ scores)
  # IMS scores = Average gene expression of 18 immunosupression genes
  # IFN-γ scores = Average gene expression of 6 INF-γ genes
  # Reference: Cui, C., Xu, C., Yang, W., Chi, Z., Sheng, X., Si, L., ... & Kong, Y. (2021). Ratio of the interferon-γ signature to the immunosuppression signature predicts anti-PD-1 therapy response in melanoma. NPJ genomic medicine, 6(1), 1-12.
  # Available from: https://www.nature.com/articles/s41525-021-00169-w
  sig1 <- IMS.Sig
  sig1 <- gsub("-", ".", sig1)

  sig2 <- IFNG.Sig
  sig2 <- gsub("-", ".", sig2)
  cat("1.13 IMS.Sig")


  getIMS.Sig <- lapply(list_train_vali_Data, function(g) {
    new <- g
    new1 <- new[, colnames(new) %in% c("Var", sig1)]

    new2 <- new[, colnames(new) %in% c("Var", sig2)]


    score1 <- rowMeans(new1[, -1])
    score2 <- rowMeans(new2[, -1])

    score <- score1 / score2

    roc <- ROCit::rocit(
      score = score,
      class = new$Var,
      negref = "Y"
    )
    auc <- roc$AUC
    return(auc)
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(names(list_train_vali_Data)) %>%
    `colnames<-`("AUC")



  AUC_IMS.Sig <- getIMS.Sig

  AUC_IMS.Sig



  auc.list.immunotherapy[["IMS.Sig"]] <- AUC_IMS.Sig


  return(auc.list.immunotherapy)
}
