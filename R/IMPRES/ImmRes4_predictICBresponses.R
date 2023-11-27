set.TCGA <- function(r.tcga = NULL, res.sig = NULL, cell.sig = NULL) {
  if (is.null(r.tcga)) {
    r.tcga <- readRDS("../Data/PublicData/TCGA_SKCM.rds")
  }
  if (is.null(res.sig)) {
    load("../Results/Signatures/resistance.program.RData")
  }
  if (is.null(cell.sig)) {
    load("../Results/Signatures/cell.type.sig.full.RData")
  }

  r.tcga <- compute.samples.res.scores(
    r = r.tcga,
    res.sig = res.sig,
    cell.sig = cell.sig,
    cc.sig = NULL,
    residu.flag = F, num.rounds = 1000
  )
  print("Saving the TCGA melanoma cohort..")
  saveRDS(r.tcga, file = "../Data/PublicData/TCGA_SKCM.rds")
  return(r.tcga)
}

prd.TCGA.survival <- function(r) {
  results <- cbind.data.frame(
    cox = cox.mat(t(r$res), r)[, "Zscore"],
    cox.CD8.T = cox.mat(t(r$res), r, X = r$tme[, "T.CD8"])[, "Zscore"]
  )
  saveRDS(results, file = "../Results/Predictors/TCGA.OS.prf.rds")
  return(results)
}

set.public.ICB.cohorts <- function() {
  load("../Results/Signatures/cell.type.sig.full.RData")
  load("../Results/Signatures/cell.type.sig.full.mouse.RData")
  load("../Results/Resistance/all.ici.sigs.RData")
  load("../Results/Resistance/all.ici.sigs.mouse.RData")
  R <- readRDS("../Data/PublicData/public.ICB.datasets.rds")
  R <- lapply(R, function(r) {
    print(paste("Processing cohort:", r$data.name))
    if (r$data.name == "aCTLA4.mouse") {
      print("Using the mouse signatures.")
      genes <- r$genes
      r$genes <- casefold(r$genes)
      r <- compute.samples.res.scores(
        r = r, res.sig = lapply(iciA.sigs.mouse, casefold),
        cell.sig = lapply(cell.sig.mouse, casefold),
        residu.flag = F, cc.sig = NULL, num.rounds = 1000
      )
      r$genes <- genes
      return(r)
    }
    r <- compute.samples.res.scores(
      r = r, res.sig = iciA.sigs, cell.sig = cell.sig,
      residu.flag = F, cc.sig = NULL, num.rounds = 1000
    )
    return(r)
  })
  print("Saving ICB cohorts..")
  saveRDS(R, file = "../Data/PublicData/public.ICB.datasets.rds")
  return(R)
}

prd.public.ICB.response <- function(R) {
  if (is.null(R)) {
    R <- get.public.ICB.cohorts(compute.scores.flag = T)
  }

  test.sig <- function(r) {
    r$cb <- is.element(r$response, c("CR", "PR"))
    b <- is.element(r$response, c("CR", "PR", "PD"))
    bA <- !is.element(r$response, "NE")
    m <- cbind.data.frame(
      CR.vs.others = apply.ttest(t(r$res[b, ]), r$response[b] == "CR")[, "zscores"],
      CR.vs.PRPD = apply.ttest(t(r$res[bA, ]), r$response[bA] == "CR")[, "zscores"],
      CRPR.vs.others = apply.ttest(t(r$res[b, ]), r$cb[b])[, "zscores"],
      CRPR.vs.PD = apply.ttest(t(r$res[bA, ]), r$cb[bA])[, "zscores"]
    )
    return(m)
  }
  rearrange <- function(m, R1, comp) {
    m1 <- t(laply(m, function(x) x[[comp]]))
    rownames(m1) <- colnames(R1[[1]]$res)
    colnames(m1) <- names(R1)
    return(m1)
  }
  idx <- setdiff(names(R), c("mapki", "aPD1.riaz"))
  ICB.prd1 <- lapply(R[idx], test.sig)
  ICB.prd2 <- lapply(colnames(ICB.prd1[[1]]), function(x) rearrange(ICB.prd1, R[idx], x))
  names(ICB.prd2) <- colnames(ICB.prd1[[1]])
  saveRDS(ICB.prd2, file = "../Results/Predictors/public.ICB.prf.rds")
  return(ICB.prd2)
}

compute.samples.res.scores <- function(r, res.sig, cell.sig, residu.flag = F,
                                       cc.sig = NULL, num.rounds = 1000) {
  r$res.ori <- NULL
  r$res <- NULL
  r$tme <- NULL
  r$X <- NULL
  names(res.sig) <- gsub(" ", ".", names(res.sig))
  two.sided <- unique(gsub(".up", "", gsub(".down", "", names(res.sig))))
  b <- is.element(paste0(two.sided, ".up"), names(res.sig)) &
    is.element(paste0(two.sided, ".down"), names(res.sig))
  two.sided <- two.sided[b]
  r$res <- get.OE.bulk(r, res.sig, num.rounds = num.rounds)
  res2 <- r$res[, paste0(two.sided, ".up")] - r$res[, paste0(two.sided, ".down")]
  if (!is.matrix(res2)) {
    res2 <- as.matrix(res2)
  }
  colnames(res2) <- two.sided
  r$res <- cbind(res2, r$res)
  rownames(r$res) <- colnames(r$tpm)
  r$tme <- get.OE.bulk(r, cell.sig, num.rounds = num.rounds)
  r$res <- cmb.res.scores(r)
  if (residu.flag) {
    print("Computing residuales")
    r$res <- t(get.residuals(t(r$res), colSums(r$tpm > 0)))
    r$tme <- t(get.residuals(t(r$tme), colSums(r$tpm > 0)))
  }
  if (!is.null(cc.sig)) {
    print("Filtering cell-cycle effects")
    if (is.null(r$cc.scores)) {
      r$cc.scores <- get.OE.bulk(r, cc.sig, num.rounds = num.rounds)
    }
    r$res.ori <- r$res
    r$res <- t(get.residuals(t(r$res), r$cc.scores))
  }
  if (is.element("resF", colnames(r$res)) & is.element("T.CD8", colnames(r$tme))) {
    r$res <- cbind.data.frame(r$res, resF.minus.TIL = r$res[, "resF"] - r$tme[, "T.CD8"])
    if (!is.null(r$res.ori)) {
      r$res.ori <- cbind.data.frame(r$res.ori, resF.minus.TIL = r$res.ori[, "resF"] - r$tme[, "T.CD8"])
    }
  }
  return(r)
}

cmb.res.scores <- function(r, res.sig = NULL, bulk.flag = T) {
  res <- r$res
  if (!is.null(res.sig)) {
    res <- get.OE(r = r, sig = res.sig, bulk.flag = bulk.flag)
  }

  if (!all(is.element(c("trt", "exc", "exc.seed"), colnames(res)))) {
    print("No merge.")
    return(res)
  }
  res <- cbind.data.frame(
    excF = rowMeans(res[, c("exc", "exc.seed")]),
    excF.up = rowMeans(res[, c("exc.up", "exc.seed.up")]),
    excF.down = rowMeans(res[, c("exc.down", "exc.seed.down")]),
    res = rowMeans(res[, c("trt", "exc", "exc.seed")]),
    res.up = rowMeans(res[, c("trt.up", "exc.up", "exc.seed.up")]),
    res.down = rowMeans(res[, c("trt.down", "exc.down", "exc.seed.down")]),
    res
  )

  if (is.element("fnc", colnames(res))) {
    res <- cbind.data.frame(
      resF = res[, "res"] + res[, "fnc"],
      resF.up = res[, "res.up"] + res[, "fnc.up"],
      resF.down = res[, "res.down"] + res[, "fnc.down"],
      res
    )
  }

  remove.sig <- add.up.down.suffix(c("exc", "exc.seed", "fnc"))
  res <- res[, setdiff(colnames(res), remove.sig)]
  colnames(res) <- gsub("excF", "exc", colnames(res))
  return(res)
}
