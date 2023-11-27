# Table S3B: Denovo cell sub/type signatures.
# Table S4A: The immune resistance program.
# Figures 1D and S1F-G: tSNE plots of non-malignant cells
# Figure 5: The clinical predictive value of the immune resistance program

write.all.tables <- function() {
  write.tableS3B_cell.subtype.signatures()
  write.tableS4A_resistance.signatures()
  generate.fig1.S1FG()
  generate.fig5()
}

write.tableS3B_cell.subtype.signatures <- function(cell.sig = NULL) {
  if (is.null(cell.sig)) {
    load("../Results/CellTypes/cell.type.sig.full.RData")
  }
  cell.sig$T.CD8.NK <- NULL
  cell.sig$lymphocyte <- NULL
  names(cell.sig) <- capitalize(gsub(".", " ", names(cell.sig), fixed = T))
  names(cell.sig)[names(cell.sig) == "Mal"] <- "Malignant"
  names(cell.sig)[names(cell.sig) == "Endo "] <- "Endothelial"
  names(cell.sig) <- toupper(names(cell.sig))
  summary(cell.sig)
  m <- list.2.mat(cell.sig)
  print("Writing TableS3B: Denovo cell sub/type signatures.")
  write.csv(m,
    row.names = F,
    file = "../Output/Tables/TableS3B_denovo.cell.subtype.sig.csv"
  )
  return(m)
}

write.tableS4A_resistance.signatures <- function(res.sig = NULL, CDK7.targets = NULL) {
  if (is.null(CDK7.targets)) {
    load("../Data/PublicData/CDK7.MYC.targets.RData")
  }
  if (is.null(res.sig)) {
    load("../Results/Signatures/resistance.program.RData")
  }

  f <- function(d) {
    idx <- paste(c("exc", "exc.seed", "trt", "fnc"), d, sep = ".")
    m <- list.2.boolean.mat(res.sig[idx])
    genes <- rownames(m)
    m <- apply(m, 2, as.integer)
    rownames(m) <- genes
    colnames(m) <- paste0(c(
      "Exclusion", "Exclusion seed",
      "Post-treatment", "Functional immune resistance"
    ), " (", d, ")")
    m <- m[order(rownames(m)), ]
    m <- cbind.data.frame(GENE = rownames(m), N = rowSums(m), m)
    m <- m[order(m$N, decreasing = T), ]
    colnames(m)[1:2] <- paste0(colnames(m)[1:2], " (", d, ")")
    m$IDX <- 1:nrow(m)
    return(m)
  }
  m.up <- f("up")
  m.up$CDK7.targets <- 1 * is.element(m.up$`GENE (up)`, CDK7.targets)
  m.down <- f("down")
  m <- merge(m.up, m.down, by = "IDX", all.x = T, all.y = T)
  m$IDX <- NULL
  print("Writing TableS4A: Immune resistance program.")
  write.csv(m, file = "../Output/Tables/TableS4A_resistance.program.csv")
  return(m)
}

generate.fig1.S1FG <- function(r = NULL, type.OE = NULL) {
  if (is.null(r)) {
    r <- readRDS("../Data/scData/Mel.all.data.QC.rds")
  }
  if (is.null(type.OE)) {
    type.OE <- readRDS("../Results/CellTypes/cell.type.OE.rds")
  }
  cd8 <- r$cd["CD8A", ] > 0 | r$cd["CD8B", ] > 0
  cd4 <- r$cd["CD4", ] > 0
  markers <- cbind.data.frame(
    CD8 = ifelse(cd8 & !cd4, "Positive", "Negative"),
    CD4 = ifelse(cd4 & !cd8, "Positive", "Negative"),
    CD19 = ifelse(r$tpm["CD19", ] > 0, "Positive", "Negative"),
    CD22 = ifelse(r$tpm["CD22", ] > 0, "Positive", "Negative")
  )
  colnames(markers) <- paste(colnames(markers), "(Fig S1G)")
  pdf("../Output/Figures/Fig1D_tSNE.nonmal.pdf")
  plot.extra(r$tsne, labels = r$cell.types, main = "Figure 1D (cell type)", set.flag = T)
  plot.extra(r$tsne, labels = r$samples, main = "Figure 1D (tumor)")
  dev.off()
  idx <- c("B.cell", "CAF", "Endo.", "Macrophage", "NK", "T.CD4", "T.CD8", "T.cell")
  X <- type.OE[, idx]
  colnames(X) <- paste(gsub(".", " ", colnames(X), fixed = T), "(Fig S1F)")
  pdf("../Output/Figures/FigS1FG_tSNE.nonmal.pdf")
  apply.plot.extra(r$tsne, labels = X)
  apply.plot.extra(r$tsne, labels = markers, set.flag = T)
  dev.off()
}

generate.fig5 <- function(r.tcga = NULL, R = NULL, r.pd1 = NULL, aPD1.val = NULL) {
  if (is.null(r.tcga)) {
    r.tcga <- readRDS("../Data/PublicData/TCGA_SKCM.rds")
    R <- readRDS("../Data/PublicData/public.ICB.datasets.rds")
    r.pd1 <- readRDS("../Data/ValidationCohorts/ValidationCohort2.rds")
    r.pd1$res <- cbind.data.frame(r.pd1$res.ori, ccf = r.pd1$res.cc)
    aPD1.val <- readRDS("../Results/Predictors/ValCoh2.prf.rds")
  }
  print("Depicting the clinical value of the immune resistance program (Figure 5).")
  aPD1.val$ttestP <- rbind.data.frame(aPD1.val$ori$ttestP,
    ccf = aPD1.val$cc$ttestP
  )
  pdf("../Output/Figures/Fig5_clinicalPred.pdf")
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  km.plot3(r.tcga, r.tcga$res[, "resu"], main = "TCGA (immune resistance)", xlim = 15, X = r.tcga$tme[, "T.CD8"], direction = 1)
  km.plot3(r.tcga, r.tcga$res[, "resu.down"], main = "TCGA (immune resistance down)", xlim = 15, X = r.tcga$tme[, "T.CD8"], direction = -1)
  km.plot3(r.tcga, r.tcga$res[, "resu"] - r.tcga$tme[, "T.CD8"],
    main = "TCGA (immune resistance\nminus T cell levels)", xlim = 15, X = r.tcga$tme[, "T.CD8"], direction = 1
  )
  mtitle("Figure 5A")
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  b <- R$aCTLA4.mouse$response != "NE"
  boxplot.test(R$aCTLA4.mouse$res[b, "res"],
    ifelse(R$aCTLA4.mouse$response[b] == "CR", "1.Responder", "2.Non-\nresponder"),
    alternative = "greater", main = "Fig5B. aCTLA4 response in mice",
    dots.flag = T
  )
  boxplot.test(R$aPD1.hugo$res[, "res"],
    ifelse(R$aPD1.hugo$response == "CR", "1.Responder", "2.Non-\nresponder"),
    alternative = "greater", main = "Fig5C. aPD1 response in melanoma patients",
    dots.flag = T
  )
  mtitle("Figure 5B-C")
  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  m <- cbind.data.frame(
    sig = c("res", "ccf.res", "resF", "resF.minus.TIL"),
    names = c(
      "Immune resistance", "Immune resistance\n(cell cycle residuals)",
      "Immune resistance (refined)",
      "Immune resistance (refined)\nminus inferred T cell levels"
    )
  )
  apply(m, 1, function(x) {
    km.plot3(r.pd1, r.pd1$res[, x[1]],
      xlim = 1.5, main = x[2], qua = 0.2,
      direction = ifelse(grepl("down", x[1]), -1, 1), ylab = "PFS",
      X = r.pd1$tme[, "T.CD8"]
    )
  })
  mtitle("Figure 5D")

  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  aPD1.plot.cb(r = r.pd1, aPD1.val = aPD1.val$ori, plot.sig = c("resF", "TME.T.CD8", "resF.minus.TIL"))
  mtitle("Figure 5F")
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  aPD1.plot.response(r = r.pd1, aPD1.val = aPD1.val$ori, plot.sig = c("res", "resF", "exc"), sd.flag = F)
  mtitle("Figure 5G")
  par(mfrow = c(1, 1), oma = c(0, 0, 3, 0))
  b <- rownames(aPD1.val$ori$cox.res) != "T.CD8"
  barplot.sig.prf(
    z = -aPD1.val$ori$cox.res$con[b], l = rownames(aPD1.val$ori$cox.res)[b],
    main = "PFS", p.type = "COX", res.sig.names = aPD1.val$ori$res.sig.names
  )
  mtitle("Figure 5E")
  barplot.sig.prf(
    z = aPD1.val$ori$ttest$CR.vs.others.zscores, l = rownames(aPD1.val$ori$ttest),
    main = "CR vs. non-CR", p.type = "t-test", res.sig.names = aPD1.val$ori$res.sig.names
  )
  mtitle("Figure 5H")
  dev.off()
  return()
}
