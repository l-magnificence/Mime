mal.t.cell.exclusion <- function(rB = NULL, r.sc = NULL, cell.sig = NULL) {
  if (is.null(rB)) {
    r.sc <- readRDS("../Data/scData/Mel.malignant.rds")
    rB <- readRDS("../Data/PublicData/TCGA_SKCM.rds")
    load("../Results/CellTypes/cell.type.sig.RData")
    results <- mal.t.cell.exclusion(rB, r.sc, cell.sig)
    return(results)
  }

  results <- cell.cell.intr(rB, r.sc,
    cellA.markers = cell.sig$Mal,
    cellB.markers = cell.sig$T.CD8,
    cellA.name = "malignant",
    cellB.name = "T.CD8",
    fileName = "FINAL_Tcell_Exclusion", pval = 0.1
  )
  exc.sig <- results$sig.final

  save(exc.sig, file = "../Results/Resistance/Exclusion/FINAL_Tcell_Exclusion.sig.RData")
  save(exc.sig, file = "../Results/Signatures/Tcell.exclusion.sig.RData")

  results <- test.seed.sig(r.sc, results, no.itr = 1000)
  saveRDS(results, file = "../Results/Resistance/Exclusion/FINAL_Tcell_Exclusion.full.rds")

  return(results)
}

cell.cell.intr <- function(rB, r.sc, cellA.markers, cellB.markers, cellA.name = "malignant",
                           cellB.name = "T.cell", bulk.confounders = NULL, sc.confounders = NULL,
                           fileName = NULL, pval = 0.1) {
  print(paste("Characterizing", cellA.name, "cells in tumors with low", cellB.name, "infiltration."))
  cellA.markers <- sort(intersect(rB$genes, cellA.markers))
  results <- list(
    cellA.markers = cellA.markers,
    cellB.markers = cellB.markers,
    bulk.confounders = bulk.confounders,
    sc.confounders = sc.confounders
  )

  print("1. Estimating cell type B abundance in bulk gene expression")
  results$bulk.cellB.abn <- round(get.OE.bulk(rB, gene.sign = list(cellB.markers), num.rounds = 1000), 2)
  rownames(results$bulk.cellB.abn) <- rB$samples

  print("2. Looking for genes correlated with cell type B abundance in bulk gene expression")
  if (is.null(bulk.confounders)) {
    results$bulk.cor <- get.cor(t(rB$tpm), results$bulk.cellB.abn, method = "pearson")
  } else {
    bulk.confounders <- as.matrix(bulk.confounders)
    b <- !is.na(rowSums(bulk.confounders))
    print(paste("Using", sum(b), "bulk samples (with full data)."))
    bulk.confounders <- bulk.confounders[b, ]
    results$bulk.cor <- pcor.mat(t(rB$tpm[, b]), results$bulk.cellB.abn[b], bulk.confounders, method = "pearson")
  }
  results$bulk.cor <- add.onesided.p(results$bulk.cor)

  print("3. Getting the seed signatures")
  results$seed.sig <- get.top.elements(results$bulk.cor[cellA.markers, c("p.pos", "p.neg")], no.elm = 20, min.cf = pval)
  names(results$seed.sig) <- c("att", "exc")
  print(summary(results$seed.sig))

  print("4. Testing if the seed signatures are anti correlated")
  results$sc.seed.scores <- round(get.OE.sc(r.sc, results$seed.sig, num.rounds = 1000), 2)
  cor.plot(results$sc.seed.scores,
    main = "Seed overall expression",
    ylab = "Seed exclusion (up)", xlab = "Seed exclusion (down)"
  )

  print("5. Expanding the seed signatures")
  r.sc$q <- cbind(r.sc$comp, log(r.sc$comp.reads))
  if (!is.null(sc.confounders)) {
    r.sc$q <- cbind(r.sc$q, sc.confounders)
    results$sc.seed.scores.rgr <- t(get.residuals(t(results$sc.seed.scores), sc.confounders))
    cor.plot(results$sc.seed.scores.rgr,
      main = "Seed residuals",
      ylab = "Seed exclusion (up)", xlab = "Seed exclusion (down)"
    )
  }

  results$att.sc <- pcor.mat(t(r.sc$tpm), results$sc.seed.scores[, 1], r.sc$q, method = "spearman")
  results$exc.sc <- pcor.mat(t(r.sc$tpm), results$sc.seed.scores[, 2], r.sc$q, method = "spearman")
  results$att.sc <- add.onesided.p(results$att.sc)
  results$exc.sc <- add.onesided.p(results$exc.sc)

  print("6. Generating the final signatures")
  f <- function(s) {
    names(s) <- gsub("p.", "", names(s), fixed = T)
    s$exc <- intersect(s$att.neg, s$exc.pos)
    s$att <- intersect(s$att.pos, s$exc.neg)
    return(s)
  }
  b.rp <- !startsWith(r.sc$genes, "RP")
  b.seed <- is.element(r.sc$genes, unlist(results$seed.sig))
  s1 <- c(
    get.top.elements(m = results$att.sc[, c("p.pos", "p.neg")], no.elm = 200, min.cf = 0.01, main = "att"),
    get.top.elements(m = results$exc.sc[, c("p.pos", "p.neg")], no.elm = 200, min.cf = 0.01, main = "exc")
  )
  s2 <- c(
    get.top.elements(m = results$att.sc[b.rp, c("p.pos", "p.neg")], no.elm = 200, min.cf = 0.01, main = "att"),
    get.top.elements(m = results$exc.sc[b.rp, c("p.pos", "p.neg")], no.elm = 200, min.cf = 0.01, main = "exc")
  )
  s3 <- c(
    get.top.elements(m = results$att.sc[b.seed, c("p.pos", "p.neg")], no.elm = 20, min.cf = 0.01, main = "att"),
    get.top.elements(m = results$exc.sc[b.seed, c("p.pos", "p.neg")], no.elm = 20, min.cf = 0.01, main = "exc")
  )
  s1 <- f(s1)
  s2 <- f(s2)
  s3 <- f(s3)
  results$sigA <- s1
  results$sig.no.rp <- s1
  results$sig <- lapply(names(s1), function(x) sort(unique(c(s1[[x]], s2[[x]], s3[[x]]))))
  names(results$sig) <- names(s1)
  print(summary(results$sig[c("exc", "att")]))

  results$sig.final <- c(results$seed.sig[c("exc", "att")], results$sig[c("exc", "att")])
  names(results$sig.final) <- c("exc.seed.up", "exc.seed.down", "exc.up", "exc.down")
  print(summary(results$sig.final))

  if (!is.null(fileName)) {
    saveRDS(results, file = paste0("../Results/Resistance/Exclusion/", fileName, ".rds"))
  }
  return(results)
}

test.seed.sig <- function(r, results, no.itr = 100) {
  no.cellA.markers <- length(results$cellA.markers)
  seed.cor <- get.cor(results$sc.seed.scores, method = "pearson")
  rnd.seed.att <- list()
  rnd.seed.exc <- list()
  # B.rnd.cellA.markers <- get.compatible.non.overlapping.sig(r = r,sig.genes = results$cellA.markers,num.rounds = no.itr)
  for (i in 1:no.itr) {
    # rnd.cellA.markers <- r$genes[B.rnd.cellA.markers[,i]]
    rnd.cellA.markers <- sample(x = r$genes, size = no.cellA.markers)
    rnd.seed <- get.top.elements(results$bulk.cor[rnd.cellA.markers, c("p.pos", "p.neg")], no.elm = 20, min.cf = 0.1)
    rnd.seed.att[[i]] <- rnd.seed$p.pos
    rnd.seed.exc[[i]] <- rnd.seed$p.neg
  }
  scores.att <- round(get.OE.sc(r, rnd.seed.att, num.rounds = 100), 2)
  scores.att[, 1] <- results$sc.seed.scores[, "att"]
  scores.exc <- round(get.OE.sc(r, rnd.seed.exc, num.rounds = 100), 2)
  scores.exc[, 1] <- results$sc.seed.scores[, "exc"]
  results$test.seed <- list()
  results$test.seed$rand.seed.cor <- get.cor(scores.att, scores.exc, method = "pearson", match.flag = T)
  hist(results$test.seed$rand.seed.cor[, "R"], 100,
    xlab = "Correlation between semi-random seed signatures",
    main = "The seed signatures are anti-correlated\n(red line)"
  )
  abline(v = results$test.seed$rand.seed.cor[1, "R"], col = "red")
  print(get.empirical.p.value(results$test.seed$rand.seed.cor[, "R"]))

  X <- 10 * ((2^r$tpm) - 1)
  r$genes.dist <- log2(rowMeans(X, na.rm = T) + 1)
  r$genes.dist.q <- discretize(r$genes.dist, n.cat = 50)
  r$zscores <- center.matrix(r$tpm, dim = 1, sd.flag = F)
  raw.scores.att <- get.semi.random.OE(
    r = r, genes.dist.q = r$genes.dist.q,
    b.sign = is.element(r$genes, results$seed.sig$att),
    num.rounds = no.itr, full.flag = T
  )
  raw.scores.exc <- get.semi.random.OE(
    r = r, genes.dist.q = r$genes.dist.q,
    b.sign = is.element(r$genes, results$seed.sig$exc),
    num.rounds = no.itr, full.flag = T
  )
  scores.att <- t(laply(1:no.itr, function(i) {
    raw.scores.att[, i] - rowMeans(raw.scores.att[, -i])
  }))
  scores.exc <- t(laply(1:no.itr, function(i) {
    raw.scores.exc[, i] - rowMeans(raw.scores.exc[, -i])
  }))

  scores.att[, 1] <- results$sc.seed.scores[, "att"]
  scores.exc[, 1] <- results$sc.seed.scores[, "exc"]

  scores.att <- round(scores.att, 2)
  scores.exc <- round(scores.exc, 2)

  results$test.seed$rand.sig.cor <- get.cor(scores.att, scores.exc, method = "pearson", match.flag = T)
  hist(results$test.seed$rand.sig.cor[, "R"], 100,
    xlab = "Correlation between semi-random signatures",
    main = "The seed signatures are anti-correlated\n(red line)"
  )
  abline(v = results$test.seed$rand.sig.cor[1, "R"], col = "red")
  results$test.seed$emp.p <- rbind(
    rand.seed = get.empirical.p.value(results$test.seed$rand.seed.cor[, "R"]),
    rand.sig = get.empirical.p.value(results$test.seed$rand.sig.cor[, "R"])
  )
  print(results$test.seed$emp.p)
  return(results)
}

get.post.trt.sig <- function(r = NULL, subs.trt = NULL) {
  if (is.null(r)) {
    load("../Results/Resistance/PostTreatment/post.treatment.subsample.de.RData")
    r <- readRDS("../Data/scData/Mel.malignant.rds")
    trt.de <- get.post.trt.sig(r, subs.trt)
    return(trt.de)
  }
  trt.de <- find.DEG(r, main = "post.treatment", mainP = "treated")
  genes <- rownames(trt.subsmp)
  sig.subsmp <- list(
    up = genes[trt.subsmp[, "F.up"] > quantile(trt.subsmp[, "F.up"], 0.9, na.rm = T)],
    down = genes[trt.subsmp[, "F.down"] > quantile(trt.subsmp[, "F.down"], 0.9, na.rm = T)]
  )
  trt.de$sig <- list(
    trt.up = sort(intersect(sig.subsmp$up, trt.de$sig$up)),
    trt.down = sort(intersect(sig.subsmp$down, trt.de$sig$down))
  )
  trt.sig <- trt.de$sig
  print(summary(trt.sig))
  save(trt.de, file = "../Results/Resistance/PostTreatment/post.treatment.de.RData")
  save(trt.sig, file = "../Results/Resistance/PostTreatment/post.treatment.sig.RData")
  save(trt.sig, file = "../Results/Signatures/post.treatment.sig.RData")
  return(trt.de)
}

get.fnc.res.sig <- function(r = NULL) {
  if (is.null(r)) {
    r <- readRDS("../Data/scData/Mel.malignant.rds")
    fnc.res.de <- get.fnc.res.sig(r)
    return(fnc.res.de)
  }
  res.down.seed <- c("B2M", "CD58", "HLA-A", "MLANA", "SOX10", "SRP54", "TAF3", "TAP2", "TAPBP")
  res.down.seed <- intersect(res.down.seed, r$genes)
  print("Characterzing the transcriptional state of malignant cells that underexpress one of these genes:")
  print(res.down.seed)

  B <- apply(r$norm[res.down.seed, ], 1, function(x) x < quantile(x, 0.01))
  print(table(rowSums(B)))
  r$res.cells <- rowSums(B) > 0
  print(paste("Found", sum(r$res.cells), "resistant cells."))
  fnc.res.de <- find.DEG(r, main = "functional.resistance", mainP = "res.cells", no.elm = 250)
  # fnc.res.de$sig<-get.scde.ttest.sig(fnc.res.de,ttest.flag = F,rank.flag = T,no.elm = 250)
  names(fnc.res.de$sig) <- c("fnc.up", "fnc.down")
  fnc.sig <- fnc.res.de$sig
  save(fnc.res.de, file = "../Results/Resistance/Functional/functional.resistance.de.RData")
  save(fnc.sig, file = "../Results/Resistance/Functional/functional.resistance.sig.RData")
  save(fnc.sig, file = "../Results/Signatures/functional.resistance.sig.RData")
  return(fnc.res.de)
}

get.res.program <- function(exc.sig = NULL, trt.sig = NULL, fnc.sig = NULL) {
  if (is.null(exc.sig)) {
    load("../Results/Signatures/Tcell.exclusion.sig.RData")
    load("../Results/Signatures/post.treatment.sig.RData")
    load("../Results/Signatures/functional.resistance.sig.RData")
  }
  res.sig <- c(
    exc.sig, trt.sig, fnc.sig,
    list(
      resu.up = c(trt.sig$trt.up, exc.sig$exc.up),
      resu.down = c(trt.sig$trt.down, exc.sig$exc.down)
    )
  )
  res.sig <- lapply(res.sig, function(x) sort(unique(x)))
  print(summary(res.sig))
  save(res.sig, file = "../Results/Signatures/resistance.program.RData")

  load("../Data/PublicData/mouse_human_mapping.RData")
  res.sig.mouse <- lapply(res.sig, function(x) convert.to.mice(mouse.human = mouse.human, x))
  saveRDS(res.sig.mouse, "../Results/Signatures/resistance.program.mouse.rds")

  write.tableS4A_resistance.signatures(res.sig)
  iciA.sig <- save.all.ICR.sig()
  return(res.sig)
}

find.DEG <- function(r, main = "", mainP = "treated", batch = NULL, no.elm = 200) {
  if (!is.null(batch)) {
    batch <- as.factor(batch)
  }
  print(paste("Running DE:", main))
  groups <- r[[mainP]]
  print(paste0(mainP, ": Positive ", sum(groups), ", negative ", sum(!groups)))
  de <- list(main = main)
  de$ttest <- apply.ttest(r$tpm, groups, ranksum.flag = F)
  de$scde <- scde.expression.difference(r$knn, r$cd, r$varinfo$prior,
    groups = factor(!groups),
    n.randomizations = 100, n.cores = 1, verbose = 1,
    batch = batch
  )
  if (!is.null(batch)) {
    de$scdeF <- de$scde
    de$scde <- de$scdeF$batch.adjusted
  }
  b <- rowSums(r$norm, na.rm = T) == 0
  r$norm[b, ] <- r$tpm[b, ]
  de$ranksum <- apply.ttest(r$norm, groups, ranksum.flag = T)
  de$sum <- cbind.data.frame(
    scde = de$scde[, c("Z", "cZ")],
    ttest = de$ttest[, c("BH.more", "BH.less")],
    ranksum = de$ranksum[, c("BH.more", "BH.less")]
  )
  de$sig <- get.scde.ttest.sig(de,
    no.elm = no.elm, scde.c = (-1.96),
    ttest.c = 0.01, rank.c = 0.01,
    ttest.flag = F, rank.flag = T
  )
  return(de)
}

get.scde.ttest.sig <- function(de, no.elm = 200, scde.c = (-1.96),
                               ttest.c = 0.01, rank.c = 0.01,
                               ttest.flag = F, rank.flag = T) {
  sig.scde <- get.top.cor(de$sum, no.elm = no.elm, min.cf = scde.c)[c("scde.cZ.up", "scde.cZ.down")]
  sig <- sig.scde
  if (ttest.flag) {
    sig.ttest <- get.top.elements(de$sum[, c("ttest.BH.more", "ttest.BH.less")],
      no.elm = no.elm, min.cf = ttest.c
    )
    sig <- sig.ttest
  }
  if (rank.flag) {
    sig.rank <- get.top.elements(de$sum[, c("ranksum.BH.more", "ranksum.BH.less")],
      no.elm = no.elm, min.cf = rank.c
    )
    sig <- sig.rank
  }
  if (ttest.flag & rank.flag) {
    sig <- list(
      up = intersect(sig.ttest$ttest.BH.more, sig.rank$ranksum.BH.more),
      down = intersect(sig.ttest$ttest.BH.less, sig.rank$ranksum.BH.less)
    )
  }
  sig <- list(
    up = intersect(sig.scde$scde.cZ.up, sig[[1]]),
    down = intersect(sig.scde$scde.cZ.down, sig[[2]])
  )
  print(summary(sig))
  return(sig)
}

cell.cell.interactions.new.dataset <- function(bulk.tpm, sc.tpm, sc.n.reads,
                                               cellA.markers, cellB.markers,
                                               cellA.name = "malignant.cell",
                                               cellB.name = "T.cell",
                                               bulk.confounders = NULL,
                                               sc.confounders = NULL,
                                               fileName = "Malignant.T.cell.exclusion",
                                               sigFilePath = "Malignant.T.cell.exclusion.signatures") {
  rB <- list(samples = colnames(bulk.tpm), genes = rownames(bulk.tpm), tpm = bulk.tpm)
  r.sc <- list(
    cells = colnames(sc.tpm), genes = rownames(sc.tpm), tpm = sc.tpm,
    comp = colSums(sc.tpm > 0), comp.reads = sc.n.reads
  )
  results <- cell.cell.intr(
    rB = rB, r.sc = r.sc,
    cellA.markers = cellA.markers,
    cellB.markers = cellB.markers,
    cellA.name = cellA.name,
    cellB.name = cellB.name,
    bulk.confounders = bulk.confounders,
    sc.confounders = sc.confounders,
    fileName = fileName, pval = 0.1
  )
  write.csv(list.2.mat(results$sig.final), file = sigFilePath)
  return(results)
}

save.all.ICR.sig <- function() {
  load("../Results/Signatures/resistance.program.RData")
  load("../Results/Signatures/cell.type.sig.full.RData")
  ICR.sig <- readRDS("../Data/PublicData/public.ICR.sig.rds")
  load("../Data/PublicData/mouse_human_mapping.RData")

  cell.sig <- cell.sig[setdiff(names(cell.sig), c("T.CD8.NK", "lymphocyte"))]
  names(cell.sig) <- paste0("TME.", names(cell.sig))
  iciA.sigs <- c(res.sig, cell.sig, ICR.sig)
  iciA.sigs.mouse <- lapply(iciA.sigs, function(x) convert.to.mice(mouse.human = mouse.human, x))
  save(iciA.sigs, file = "../Results/Resistance/all.ici.sigs.RData")
  save(iciA.sigs.mouse, file = "../Results/Resistance/all.ici.sigs.mouse.RData")
  return(iciA.sigs)
}
