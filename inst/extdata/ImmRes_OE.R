get.OE.sc <- function(r, gene.sign = NULL, num.rounds = 1000, mat.return.flag = T) {
  set.seed(1234)
  r$genes.mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm, 1, r$genes.mean, FUN = "-")
  if (any(r$tpm < 0)) {
    print("Using counts to bin genes!!!")
    r$genes.dist <- rowMeans(r$cd > 0)
  } else {
    X <- 10 * ((2^r$tpm) - 1)
    r$genes.dist <- log2(rowMeans(X, na.rm = T) + 1)
  }
  r$genes.dist.q <- discretize(r$genes.dist, n.cat = 50)
  r$sign.scores <- matrix(data = 0, nrow = ncol(r$tpm), ncol = length(gene.sign))
  sign.names <- names(gene.sign)
  colnames(r$sign.scores) <- sign.names
  sig.names <- names(gene.sign)
  r$sign.scores.raw <- r$sign.scores
  for (i in 1:length(gene.sign)) {
    b.sign <- is.element(r$genes, gene.sign[[i]])
    if (sum(b.sign) > 1) {
      rand.scores <- get.semi.random.OE(r, r$genes.dist.q, b.sign, num.rounds = num.rounds)
      raw.scores <- colMeans(r$zscores[b.sign, ])
      final.scores <- raw.scores - rand.scores
      r$sign.scores[, i] <- final.scores
      r$sign.scores.raw[, i] <- raw.scores
    }
  }
  if (mat.return.flag) {
    sign.scores <- r$sign.scores
    return(sign.scores)
  } else {
    return(r)
  }
}

get.OE.bulk <- function(r, gene.sign = NULL, num.rounds = 1000, full.flag = F) {
  set.seed(1234)
  r$genes.mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm, 1, r$genes.mean, FUN = "-")
  r$genes.dist <- r$genes.mean
  r$genes.dist.q <- discretize(r$genes.dist, n.cat = 50)
  r$sig.scores <- matrix(data = 0, nrow = ncol(r$tpm), ncol = length(gene.sign))
  sig.names <- names(gene.sign)
  colnames(r$sig.scores) <- sig.names
  r$sig.scores.raw <- r$sig.scores
  rand.flag <- is.null(r$rand.scores) | !all(is.element(names(gene.sign), colnames(r$rand.scores)))
  if (rand.flag) {
    print("Computing also random scores.")
    r$rand.scores <- r$sig.scores
  }
  for (i in sig.names) {
    b.sign <- is.element(r$genes, gene.sign[[i]])
    if (sum(b.sign) < 2) {
      next()
    }
    if (rand.flag) {
      rand.scores <- get.semi.random.OE(r, r$genes.dist.q, b.sign, num.rounds = num.rounds)
    } else {
      rand.scores <- r$rand.scores[, i]
    }
    raw.scores <- colMeans(r$zscores[b.sign, ])
    final.scores <- raw.scores - rand.scores
    r$sig.scores[, i] <- final.scores
    r$sig.scores.raw[, i] <- raw.scores
    r$rand.scores[, i] <- rand.scores
  }
  if (full.flag) {
    return(r)
  }
  sig.scores <- r$sig.scores
  return(sig.scores)
}

get.semi.random.OE <- function(r, genes.dist.q, b.sign, num.rounds = 1000, full.flag = F) {
  # Previous name: get.random.sig.scores
  sign.q <- as.matrix(table(genes.dist.q[b.sign]))
  q <- rownames(sign.q)
  idx.all <- c()
  B <- matrix(data = F, nrow = length(genes.dist.q), ncol = num.rounds)
  Q <- matrix(data = 0, nrow = length(genes.dist.q), ncol = num.rounds)
  for (i in 1:nrow(sign.q)) {
    num.genes <- sign.q[i]
    if (num.genes > 0) {
      idx <- which(is.element(genes.dist.q, q[i]))
      for (j in 1:num.rounds) {
        idxj <- sample(idx, num.genes)
        Q[i, j] <- sum(B[idxj, j] == T)
        B[idxj, j] <- T
      }
    }
  }
  rand.scores <- apply(B, 2, function(x) colMeans(r$zscores[x, ]))
  if (full.flag) {
    return(rand.scores)
  }
  rand.scores <- rowMeans(rand.scores)
  return(rand.scores)
}

get.seurat.dataset <- function(genes, D1.data, D1.name, find.var = T) {
  print(paste("Processing", D1.name, "..."))
  D1.data <- round(D1.data, 2)
  D1.data <- D1.data[genes, !duplicated(colnames(D1.data))]
  dim(D1.data)
  D1 <- CreateSeuratObject(raw.data = D1.data)
  D1 <- NormalizeData(object = D1)
  D1 <- ScaleData(object = D1)
  if (find.var) {
    D1 <- FindVariableGenes(object = D1, do.plot = FALSE)
  }
  D1@meta.data$source <- D1.name
  return(D1)
}

get.OE <- function(r, sig, bulk.flag = F, num.rounds = 1000) {
  if (bulk.flag) {
    print("Compute bulk scores")
    scores <- get.OE.bulk(r, sig, num.rounds = num.rounds)
  } else {
    print("Compute single-cell scores")
    scores <- get.OE.sc(r, sig, num.rounds = num.rounds)
  }

  names(sig) <- gsub(" ", ".", names(sig))
  two.sided <- unique(gsub(".up", "", gsub(".down", "", names(sig))))
  b <- is.element(paste0(two.sided, ".up"), names(sig)) &
    is.element(paste0(two.sided, ".down"), names(sig))
  if (any(b)) {
    two.sided <- two.sided[b]
    scores2 <- as.matrix(scores[, paste0(two.sided, ".up")] - scores[, paste0(two.sided, ".down")])
    colnames(scores2) <- two.sided
    scores <- cbind(scores2, scores)
  }

  if (!is.null(r$cells)) {
    rownames(scores) <- r$cells
  } else {
    if (!is.null(r$samples)) {
      rownames(scores) <- r$samples
    }
  }
  return(scores)
}

find.cycling.cells <- function(r, cc.sig = NULL) {
  if (is.null(cc.sig)) {
    load("../Data/PublicData/cell.cycle.sig.RData")
  }
  r$cc.scores <- get.OE.sc(r, cc.sig)
  r$ccB <- apply(r$cc.scores, 2, function(x) {
    plot.bimodal.distribution(x,
      pos.label = "cycling",
      neg.label = "non.cycling"
    )$labelsF
  })
  r$cc <- r$ccB[, "G1_S"] == "cycling" | r$ccB[, "G2_M"] == "cycling"
  print(paste("Found", sum(r$cc), "cycling cells"))
  plot.extra(r$cc.scores, labels = r$cc, xlab = "G1/S phase", ylab = "G2/M phase")
  return(r)
}

get.OE.bulk.specific.tpm <- function(r, gene.sign = NULL, num.rounds = 1000) {
  # Previous name: get.sign.scores.bulk
  set.seed(1234)
  r$genes.mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm, 1, r$genes.mean, FUN = "-")
  r$genes.dist <- r$genes.mean
  r$genes.dist.q <- discretize(r$genes.dist, n.cat = 50)
  r$sign.scores <- matrix(data = 0, nrow = ncol(r$tpm), ncol = length(gene.sign))
  sig.names <- names(gene.sign)
  colnames(r$sign.scores) <- sig.names
  r$sign.scores.raw <- r$sign.scores
  for (i in sig.names) {
    b.sign <- is.element(r$genes, gene.sign[[i]])
    if (sum(b.sign) < 2) {
      next()
    }
    rand.scores <- r$rand.scores[, i]
    raw.scores <- colMeans(r$zscores[b.sign, ])
    final.scores <- raw.scores - rand.scores
    r$sign.scores[, i] <- final.scores
    r$sign.scores.raw[, i] <- raw.scores
  }
  sign.scores <- r$sign.scores
  return(sign.scores)
}
