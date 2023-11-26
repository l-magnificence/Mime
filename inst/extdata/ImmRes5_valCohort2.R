aPD1.compute.res.scores <- function(r = NULL, res.sig = NULL, cc.flag = F) {
  # This function computes the OE of the immune resistance program
  # (identified in Jerby et al. 2018) based on the TPM matrix of the program's genes.
  # The TPM data is provided through the single-cell portal:
  # https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download
  # File name: ValCo2_tpm_resistance.csv

  if (is.null(res.sig)) {
    load("../Results/Signatures/resistance.program.RData")
  }
  if (is.null(r)) {
    r <- readRDS("../Data/ValidationCohorts/ValidationCohort2.rds")
  }
  # r$tpm is the matrix provided in the protal (see above).
  r$res <- get.OE(r, sig = res.sig, bulk.flag = T)
  r$res <- t(get.residuals(t(r$res), r$n.genes))
  resistance.scores <- cmb.res.scores(r)
  if (cc.flag) {
    resistance.scores <- t(get.residuals(t(resistance.scores), r$cc.scores))
  }
  return(resistance.scores)
}

prd.aPD1 <- function(r = NULL, res.sig, fileName = "tmp") {
  if (is.null(r)) {
    load("../Results/Signatures/resistance.program.RData")
    r <- readRDS("../Data/ValidationCohorts/ValidationCohort2.rds")
    aPD1.val <- list()
    r$res <- r$res.cc
    aPD1.val$cc <- prd.aPD1(r, res.sig, fileName = "cc")
    r$res <- r$res.ori
    aPD1.val$ori <- prd.aPD1(r, res.sig, fileName = "ori")
    saveRDS(aPD1.val, "../Results/Predictors/ValCoh2.prf.rds")
    return(aPD1.val)
  }

  aPD1.val <- list()
  b1 <- is.element(r$cb, c("non-CB", "CB"))
  b2 <- r$response != "na" & !is.na(r$response) & r$response != "MR"
  b3 <- is.element(r$cbs, c("non-CB", "CB>1yr"))
  aPD1.val$ttest <- cbind.data.frame(
    CR.vs.PR.PD = apply.ttest(t(r$res[b1, ]), r$response[b1] == "CR", ranksum.flag = F),
    CR.vs.others = apply.ttest(t(r$res[b2, ]), r$response[b2] == "CR", ranksum.flag = F),
    CB.vs.nonCB = apply.ttest(t(r$res[b1, ]), r$cb[b1] == "CB", ranksum.flag = F),
    longCB.vs.nonCB = apply.ttest(t(r$res[b3, ]), r$cbs[b3] == "CB>1yr", ranksum.flag = F)
  )

  aPD1.val$auc <- cbind.data.frame(
    CR.vs.PR.PD = get.auc.mat(r$res[b1, ], r$response[b1] == "CR"),
    CR.vs.others = get.auc.mat(r$res[b2, ], r$response[b2] == "CR"),
    CB.vs.nonCB = get.auc.mat(r$res[b1, ], r$cb[b1] == "CB"),
    longCB.vs.others = get.auc.mat(r$res[b2, ], r$cbs[b2] == "CB>1yr")
  )
  aPD1.val$auc[aPD1.val$auc < 0.5] <- 1 - aPD1.val$auc[aPD1.val$auc < 0.5]

  load("../Data/scData/res.sig.names.RData")
  aPD1.val$res.sig.names <- res.sig.names
  b.denovo <- is.element(rownames(aPD1.val$ttest), res.sig.names)
  b.denovo.down <- is.element(rownames(aPD1.val$ttest), res.sig.names.down)

  idx <- unique(multi.gsub(c(".more", ".less", ".zscores"), "", colnames(aPD1.val$ttest)))
  aPD1.val$ttestP <- aPD1.val$ttest[, paste0(idx, ".less")]
  aPD1.val$ttestP[b.denovo.down, ] <- aPD1.val$ttest[b.denovo.down, paste0(idx, ".more")]
  aPD1.val$ttestP <- aPD1.val$ttestP[b.denovo, ]
  colnames(aPD1.val$ttestP) <- gsub(".less", "", colnames(aPD1.val$ttestP))

  aPD1.plot(r, aPD1.val, fileName = fileName)
  saveRDS(aPD1.val, file = paste0("../Results/Predictors/ValCoh2.prf.", fileName, ".rds"))
  return(aPD1.val)
}

set.aPD1 <- function(r = NULL, iciA.sigs, cell.sig) {
  if (is.null(r)) {
    r <- readRDS("../Data/ValidationCohorts/ValidationCohort2.rds")
    load("../Results/Resistance/all.ici.sigs.RData")
    load("../Results/CellTypes/cell.type.sig.RData")
  }
  r <- compute.samples.res.scores(
    r = r,
    res.sig = iciA.sigs,
    cell.sig = cell.sig,
    residu.flag = T,
    cc.sig = iciA.sigs[c("G1 S (Tirosh)", "G2 M (Tirosh)")],
    num.rounds = 1000
  )
  r$res.cc <- r$res
  r$sampleName <- "ValidationCohort2_aPD1"
  print("Saving Validation Cohort II..")
  saveRDS(r, file = "../Data/ValidationCohorts/ValidationCohort2.rds")
  return(r)
}

aPD1.plot <- function(r, aPD1.val, fileName = "") {
  load("../Data/scData/res.sig.names.RData")
  fileName <- paste0("~/Desktop/CellRevisions/Figures/AutoFig/Fig5_ValCoh2.prf.", fileName, ".pdf")
  r$X <- r$res[, intersect(colnames(r$res), aPD1.val$res.sig.names)]
  pdf(fileName)

  barplot.sig.prf(
    z = aPD1.val$ttest$CR.vs.others.zscores,
    l = rownames(aPD1.val$ttest),
    main = "CR vs. non-CR",
    p.type = "t-test",
    res.sig.names = res.sig.names
  )

  b <- rownames(aPD1.val$cox.res) != "T.CD8"
  barplot.sig.prf(
    z = -aPD1.val$cox.res$con[b],
    l = rownames(aPD1.val$cox.res)[b],
    main = "PFS",
    p.type = "COX",
    res.sig.names = res.sig.names
  )

  aPD1.plot.response(r = r, aPD1.val = aPD1.val, plot.sig = c("res", "resF", "TME.T.CD8", "resF.minus.TIL"), sd.flag = F)
  aPD1.plot.cb(r = r, aPD1.val = aPD1.val, plot.sig = c("res", "resF", "TME.T.CD8", "resF.minus.TIL"))
  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  lapply(c("res.up", "res.down", "resF.up", "resF.down", "TME.T.CD8", "resF.minus.TIL"), function(x) {
    km.plot3(r, r$res[, x], xlim = 2, main = x, qua = 0.2, direction = ifelse(grepl("down", x), -1, 1), ylab = "PFS", X = r$tme[, "T.CD8"])
    print(x)
    print(ifelse(grepl("down", x), -1, 1))
  })
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  lapply(c("res", "resF", "TME.T.CD8", "resF.minus.TIL"), function(x) {
    km.plot3(r, r$res[, x], xlim = 2, main = x, qua = 0.2, direction = ifelse(grepl("down", x), -1, 1), ylab = "PFS", X = r$tme[, "T.CD8"])
  })

  aPD1.plot.response(r = r, aPD1.val = aPD1.val, plot.sig = res.sig.names, sd.flag = F)
  aPD1.plot.cb(r = r, aPD1.val = aPD1.val, plot.sig = res.sig.names)
  aPD1.plot.km(r = r, plot.sig = res.sig.names, res.sig.names.down = res.sig.names.down)
  dev.off()
  return()
}

aPD1.plot.cb <- function(r, aPD1.val, plot.sig) {
  plot.sig <- intersect(plot.sig, colnames(r$res))
  b <- is.element(r$response, c("PR", "CR", "PD"))
  Xf <- rbind(r$res[b, ], r$res[r$cb == "CB", ])
  l.cb <- c(r$cbs[b], r$cb[r$cb == "CB"])
  l.cb[l.cb == "non-CB"] <- "1.non-CB"
  l.cb[l.cb == "CB"] <- "2.CB"
  l.cb[l.cb == "CB<6m"] <- "3.CB<6m"
  l.cb[l.cb == "CB<1yr"] <- "4.CB<1yr"
  l.cb[l.cb == "CB>1yr"] <- "5.CB>1yr"

  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  laply(plot.sig, function(x) {
    print(x)
    boxplot.test(Xf[, x], l.cb,
      dots.flag = T, las = 2, legend.flag = F, t.test.flag = "none", ylab = "OE",
      main = paste0(
        x, "\n",
        format.pval.private(aPD1.val$ttestP[x, c("CB.vs.nonCB", "longCB.vs.nonCB")]),
        "\nAUC = ", round(aPD1.val$auc[x, "longCB.vs.others"], 2)
      )
    )
  })
  return()
}

aPD1.plot.response <- function(r, aPD1.val, plot.sig, sd.flag = F) {
  plot.sig <- intersect(plot.sig, colnames(r$res))
  l.rs <- r$response
  l.rs[l.rs == "CR"] <- "1.CR"
  l.rs[l.rs == "PD"] <- "3.PD"
  if (sd.flag) {
    valid.rsp <- c("CR", "PR", "SD", "PD")
    l.rs[is.element(l.rs, c("PR", "SD"))] <- "2.PR/SD"
  } else {
    valid.rsp <- c("CR", "PR", "PD")
    l.rs[l.rs == "PR"] <- "2.PR"
  }
  b <- is.element(r$response, valid.rsp)
  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  laply(plot.sig, function(x) {
    boxplot.test(r$res[b, x], l.rs[b],
      dots.flag = T, las = 2, legend.flag = F, t.test.flag = "none", ylab = "OE",
      main = paste0(
        x, "\n",
        format.pval.private(aPD1.val$ttestP[x, "CR.vs.PR.PD"]),
        "\nAUC = ", round(aPD1.val$auc[x, "CR.vs.others"], 2)
      )
    )
  })
  return()
}

barplot.sig.prf <- function(z, l, main, res.sig.names, p.type = "t-test") {
  # redundant.res.sig<-add.up.down.suffix(c('exc',"trt","resu"))
  redundant.res.sig <- add.up.down.suffix(c("resu"))
  res.sig.names <- setdiff(res.sig.names, "TME.T.CD8")
  b <- !is.element(l, redundant.res.sig)
  z <- z[b]
  l <- l[b]
  idx <- order(z)
  z <- z[idx]
  l <- l[idx]
  p.cut <- (-log10(0.05))
  b.icr <- is.element(l, res.sig.names)
  p <- wilcox.test(abs(z)[b.icr], abs(z)[!b.icr], alternative = "greater")$p.value
  l <- gsub(".", " ", l, fixed = T)
  l <- capitalize(l)
  barplot(
    height = abs(z), col = ifelse(z > 0, "lightblue", "gray"),
    names.arg = l, las = 2, cex.names = 0.6, border = ifelse(b.icr, "black", "white"),
    xlab = paste("-log10(", p.type, "p-value)"), horiz = T, main = format.pval.private(p)
  )
  abline(v = (-log10(0.05)), col = "red", lty = 2)
}

aPD1.plot.km <- function(r, plot.sig, res.sig.names.down) {
  plot.sig <- intersect(plot.sig, colnames(r$res))
  l <- unique(gsub(".down", "", gsub(".up", "", plot.sig)))
  b2 <- is.element(paste0(l, ".up"), plot.sig) & is.element(paste0(l, ".down"), plot.sig)
  qua <- 0.2
  lapply(l[b2], function(x1) {
    par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
    x <- paste0(x1, ".up")
    km.plot3(r, r$res[, x], xlim = 2, main = x, qua = qua, direction = 1, ylab = "PFS", X = r$tme[, "T.CD8"])
    x <- paste0(x1, ".down")
    km.plot3(r, r$res[, x], xlim = 2, main = x, qua = qua, direction = (-1), ylab = "PFS", X = r$tme[, "T.CD8"])
    x <- x1
    km.plot3(r, r$res[, x], xlim = 2, main = x, qua = qua, direction = 1, ylab = "PFS", X = r$tme[, "T.CD8"])
  })

  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  lapply(l[!b2], function(x) {
    km.plot3(r, r$res[, x],
      main = x, X = r$tme[, "T.CD8"], xlim = 2, qua = 0.2,
      direction = ifelse(is.element(x, res.sig.names.down), -1, 1), ylab = "PFS"
    )
  })
}
