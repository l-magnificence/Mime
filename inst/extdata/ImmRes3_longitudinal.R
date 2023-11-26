set.ValCo1 <- function() {
  r <- readRDS("../Data/ValidationCohorts/ValidationCohort1.rds")
  load("../Results/Signatures/cell.type.sig.full.RData")
  load("../Results/Signatures/resistance.program.RData")
  r <- compute.samples.res.scores(
    r = r, res.sig = res.sig,
    cell.sig = cell.sig[c("B.cell", "CAF", "Macrophage", "T.cell", "Mal")],
    residu.flag = F, cc.sig = NULL, num.rounds = 1000
  )
  saveRDS(r, file = "../Data/ValidationCohorts/ValidationCohort1.rds")
  return(r)
}

test.ValCo1 <- function(r) {
  if (is.null(r)) {
    r <- readRDS("../Data/ValidationCohorts/ValidationCohort1.rds")
  }
  r$residu <- t(get.residuals(t(r$res), r$tme))
  r$ici <- r$B.ici[, "on"] | r$B.ici[, "post"]
  r$tr <- r$B.target[, "on"] | r$B.target[, "post"]
  r$immu <- (r$B.immun[, "on"] | r$B.immun[, "post"]) & !r$ici
  r$ici.on <- r$B.ici[, "on"]
  r$ici.post <- r$B.ici[, "post"]
  r$tr.on <- r$B.target[, "on"]
  r$tr.post <- r$B.target[, "post"]

  results <- list(main = "Validation Cohort 1")
  results$anova <- cbind.data.frame(
    original = apply(r$res, 2, function(x) get.anova.p(x, r$patients)),
    TME.filtered = (apply(r$residu, 2, function(x) get.anova.p(x, r$patients)))
  )

  f1 <- function(y) {
    r$y <- y
    M1 <- with(r, lmer(y ~ (1 | patients) + tme + ici + tr + immu, mtcars))
    idx <- c("iciTRUE", "trTRUE", "immuTRUE")
    v <- summary(M1)$coefficients[idx, c("Estimate", "Pr(>|t|)")]
    z <- get.cor.zscores(v[, "Estimate"], v[, "Pr(>|t|)"])
    return(z)
  }
  f2 <- function(y) {
    r$y <- y
    M1 <- with(r, lmer(y ~ (1 | patients) + ici.on + ici.post + tr.on + tr.post + immu, mtcars))
    idx <- c("ici.onTRUE", "ici.postTRUE", "tr.onTRUE", "tr.postTRUE")
    v <- summary(M1)$coefficients[idx, c("Estimate", "Pr(>|t|)")]
    z <- get.cor.zscores(v[, "Estimate"], v[, "Pr(>|t|)"])
    return(z)
  }
  results$hlm.res <- t(apply(r$res, 2, f1))
  colnames(results$hlm.res) <- gsub("TRUE", "", colnames(results$hlm.res))
  results$hlm.res <- cbind.data.frame(ICI.p = 10^-abs(results$hlm[, "ici"]), results$hlm.res)

  results$hlm.tme <- t(apply(r$tme, 2, f2))
  colnames(results$hlm.tme) <- gsub("TRUE", "", colnames(results$hlm.tme))
  results$hlm.tme <- cbind.data.frame(ICI.on.p = 10^-abs(results$hlm.tme[, "ici.on"]), results$hlm.tme)

  saveRDS(results, file = "../Results/Predictors/ValCoh1.prf.rds")

  print(paste(
    "The immune resistance program is induced in resistant on/post-ICI tumors:",
    format.pval.private(results$hlm.res[c("res"), "ICI.p"])
  ))
  print(paste(
    "The refined immune resistance program is induced in resistant on/post-ICI tumors:",
    format.pval.private(results$hlm.res[c("resF"), "ICI.p"])
  ))
  print(paste("T cell infiltration on ICI :", format.pval.private(results$hlm.tme["T.cell", "ICI.on.p"])))
  print(results$hlm.res[c("resF", "res"), ])
  print(results$hlm.tme["T.cell", ])
  return(results)
}

set.matched.MAPKi.Hugo <- function() {
  r <- readRDS("../Data/PublicData/MAPKi.Hugo.Cell.2015.rds")
  load("../Results/Signatures/resistance.program.RData")
  load("../Results/Signatures/cell.type.sig.RData")
  mapki.sig <- readRDS("../Data/PublicData/public.ICR.sig.rds")[c("mapki.res.up", "mapki.res.down")]
  r <- compute.samples.res.scores(
    r = r, res.sig = res.sig, cell.sig = cell.sig,
    residu.flag = F, cc.sig = NULL, num.rounds = 1000
  )
  r$mapki <- get.OE(r = r, sig = mapki.sig, bulk.flag = T, num.rounds = 1000)
  saveRDS(r, file = "../Data/PublicData/MAPKi.Hugo.Cell.2015.rds")
  return(r)
}

test.matched.MAPKi.Hugo <- function(r) {
  if (is.null(r)) {
    r <- readRDS("../Data/PublicData/MAPKi.Hugo.Cell.2015.rds")
  }
  r$tme <- r$tme[, c("B.cell", "CAF", "Macrophage", "T.cell", "Mal")]
  f <- function(y) {
    r$y <- y
    M1 <- with(r, lmer(y ~ prog + (1 | patient) + tme, mtcars))
    v <- summary(M1)$coefficient["progTRUE", c("Estimate", "Pr(>|t|)")]
    return(v)
  }
  results <- t(apply(cbind.data.frame(r$res, r$mapki), 2, f))
  print(results[c("res", "resF", "mapki.res"), ])
  saveRDS(results, file = "../Results/Predictors/MAPKi.res.Hugo2015.prf.rds")
  return(results)
}
