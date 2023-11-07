apply.ttest <- function(m, b, ranksum.flag = F) {
  if (ranksum.flag) {
    p <- t(apply(m, 1, function(x) {
      c(
        wilcox.test(x[b], x[!b], alternative = "greater")$p.value,
        wilcox.test(x[b], x[!b], alternative = "less")$p.value
      )
    }))
  } else {
    p <- t(apply(m, 1, function(x) {
      c(
        t.test(x[b], x[!b], alternative = "greater")$p.value,
        t.test(x[b], x[!b], alternative = "less")$p.value
      )
    }))
  }
  colnames(p) <- c("more", "less")
  p <- cbind(p, get.p.zscores(p))
  colnames(p)[3] <- "zscores"
  p <- cbind.data.frame(p, BH = p.adjust.mat(p[, 1:2], method = "BH"))
  return(p)
}

ttest.class <- function(v, b, alternative = "two.sided") {
  p <- t.test(v[b], v[!b], alternative = alternative)$p.value
  return(p)
}

wilcox.class <- function(v, b, alternative = "two.sided") {
  p <- NA
  v <- v[!is.na(b)]
  b <- b[!is.na(b)]
  v1 <- v[b]
  v2 <- v[!b]
  if (sum(!is.na(v1)) == 0 | sum(!is.na(v2)) == 0) {
    return(p)
  }
  p <- wilcox.test(v1, v2, alternative = alternative)$p.value
  return(p)
}

get.p.zscores <- function(p) {
  b <- p[, 1] > 0.5
  b[is.na(b)] <- F
  zscores <- (-log10(p[, 1]))
  zscores[b] <- log10(p[b, 2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}

get.cor.zscores <- function(R, P) {
  p <- cbind(get.onesided.p.value(R, P), get.onesided.p.value(-R, P))
  z <- get.p.zscores(p)
  return(z)
}

add.onesided.p <- function(m) {
  m <- cbind.data.frame(m,
    p.pos = get.onesided.p.value(m[, "R"], m[, "P"]),
    p.neg = get.onesided.p.value(-m[, "R"], m[, "P"])
  )
  m[, c("p.pos", "p.neg")] <- p.adjust.mat(m[, c("p.pos", "p.neg")], method = "BH")
  return(m)
}

get.onesided.p.value <- function(c, p) {
  p[p == 0] <- min(p[p > 0], na.rm = T)
  p.one.side <- p
  p.one.side[] <- NA
  b <- c > 0 & !is.na(c)
  p.one.side[b] <- p[b] / 2
  b <- c <= 0 & !is.na(c)
  p.one.side[b] <- 1 - (p[b] / 2)
  return(p.one.side)
}

p.adjust.mat <- function(m, method = "BH") {
  P <- apply(m, 2, function(x) p.adjust(x, method = method))
  return(P)
}

get.top.elements <- function(m, no.elm = 100, min.cf = Inf, main = "") {
  q <- no.elm / nrow(m)
  top.l <- lapply(colnames(m), function(x) {
    v <- m[, x]
    cf <- min(quantile(x = v, probs = q, na.rm = T), min.cf)
    b <- v <= cf
    elm <- sort(rownames(m)[b])
    return(elm)
  })
  if (main != "") {
    main <- paste0(main, ".")
  }
  names(top.l) <- paste0(main, colnames(m))
  return(top.l)
}

list.2.mat <- function(l) {
  n1 <- max(laply(l, length))
  m <- t(laply(l, function(x) c(x, matrix(data = "", nrow = n1 - length(x) + 1))))
  m <- m[1:n1, ]
  colnames(m) <- names(l)
  return(m)
}

get.mat <- function(m.rows, m.cols, data = NA) {
  m <- matrix(
    data = data, nrow = length(m.rows), ncol = length(m.cols),
    dimnames = list(m.rows, m.cols)
  )
  return(m)
}

get.empirical.p.value <- function(p) {
  if (is.matrix(p)) {
    p.emp <- apply(p, 1, function(x) sum(x <= x[1], na.rm = T))
    p.emp[is.na(p[, 1])] <- NA
    no.permut <- rowSums(!is.na(p))
    p.emp <- p.emp / no.permut
  } else {
    p.rand <- p[-1]
    p.rand <- p.rand[!is.na(p.rand)]
    p.emp <- c(
      emp.p = mean(p.rand <= p[1]),
      sd.dist = (p[1] - mean(p.rand)) / sd(p.rand),
      mean.rand = mean(p.rand),
      sd.rand = sd(p.rand)
    )
  }
  return(p.emp)
}

discretize <- function(v, n.cat) {
  q1 <- quantile(v, seq(from = (1 / n.cat), to = 1, by = (1 / n.cat)))
  u <- matrix(nrow = length(v))
  for (i in 2:n.cat) {
    u[(v >= q1[i - 1]) & (v < q1[i])] <- i
  }
  return(u)
}

get.cor <- function(v1, v2 = NULL, method = "spearman",
                    use = "pairwise.complete.obs", match.flag = F,
                    alternative = "two.sided") {
  # Previous name: spearman.cor
  if (is.null(v2)) {
    v2 <- v1
  }
  if (!is.matrix(v1)) {
    v1 <- as.matrix(v1)
  }
  if (!is.matrix(v2)) {
    v2 <- as.matrix(v2)
  }
  if (is.null(colnames(v1))) {
    colnames(v1) <- 1:ncol(v1)
  }
  if (is.null(colnames(v2))) {
    colnames(v2) <- 1:ncol(v2)
  }
  if (match.flag) {
    results <- laply(colnames(v1), function(i) {
      c.i <- cor.test(v1[, i], v2[, i], method = method, use = use, alternative = alternative)
      p <- c(c.i$estimate, c.i$p.value)
      return(p)
    })
    colnames(results) <- c("R", "P")
    return(results)
  }
  m <- get.mat(m.rows = colnames(v1), m.cols = colnames(v2))
  results <- list(cor = m, p = m)
  for (i in 1:ncol(v1)) {
    f <- function(x) {
      c.i <- cor.test(v1[, i], x, method = method, use = use, alternative = alternative)
      c(c.i$estimate, c.i$p.value)
    }
    c.i <- apply(v2, 2, f)
    results$cor[i, ] <- c.i[1, ]
    results$p[i, ] <- c.i[2, ]
  }
  if (ncol(v2) == 1) {
    results <- cbind(results$cor, results$p)
    colnames(results) <- c("R", "P")
  }
  return(results)
}

pcor.mat <- function(v1, v2 = NULL, v3, method = "spearman",
                     use = "pairwise.complete.obs", match.flag = F,
                     alternative = "two.sided") {
  f1 <- function(x, y, z) {
    c.i <- pcor.test(x = x, y = y, z = z, method = method)
    p <- c(c.i$estimate, c.i$p.value)
    return(p)
  }
  f2 <- function(x, y, z) {
    p <- tryCatch(f1(x, y, z), error = function(err) {
      return(c(NA, NA))
    })
    return(p)
  }
  if (is.null(v2)) {
    v2 <- v1
  }
  if (!is.matrix(v1)) {
    v1 <- as.matrix(v1)
  }
  if (!is.matrix(v2)) {
    v2 <- as.matrix(v2)
  }
  if (!is.matrix(v3)) {
    v3 <- as.matrix(v3)
  }
  if (match.flag) {
    n <- ncol(v1)
    results <- matrix(data = NA, nrow = n, ncol = 2)
    rownames(results) <- colnames(v1)
    for (i in 1:ncol(v1)) {
      b <- !is.na(v1[, i]) & !is.na(v2[, i])
      v1i <- v1[b, i]
      v2i <- v2[b, i]
      results[i, ] <- f2(v1i, v2i, v3[b, ])
    }
    return(results)
  }
  m <- matrix(nrow = ncol(v1), ncol = ncol(v2))
  rownames(m) <- colnames(v1)
  colnames(m) <- colnames(v2)
  results <- list(cor = m, p = m)
  for (i in 1:ncol(v1)) {
    f <- function(x) {
      b <- !is.na(v1[, i]) & !is.na(x)
      c.i <- f2(v1[b, i], x[b], v3[b, ])
      return(c.i)
    }
    c.i <- apply(v2, 2, f)
    results$cor[i, ] <- c.i[1, ]
    results$p[i, ] <- c.i[2, ]
  }
  if (ncol(v2) == 1) {
    results <- cbind(results$cor, results$p)
    colnames(results) <- c("R", "P")
  }
  return(results)
}

remove.ribo <- function(l) {
  if (is.list(l)) {
    l <- lapply(l, function(x) x[!startsWith(x, "RP")])
  } else {
    l <- l[!startsWith(l, "RP")]
  }
  return(l)
}

cor.plot <- function(x, y = NULL, main = "", ylab = "", xlab = "", regression.flag = F) {
  if (is.null(y)) {
    v <- colnames(x)
    xlab <- v[1]
    ylab <- v[2]
    y <- x[, 2]
    x <- x[, 1]
  }
  v <- get.cor(x, y)
  main <- paste(main, "\nR =", format(v[1], di = 2), "P =", format(v[2], scientific = T, di = 2))
  plot(x, y, main = main, xlab = xlab, ylab = ylab, cex = 0.3)
  b <- !is.na(x) & !is.na(y)
  v <- lowess(x[b], y[b])
  lines(v, col = "red")
  if (!regression.flag) {
    return()
  }
  y.d <- y - v$y[match(x, v$x)]
  y.sd <- sd(y.d, na.rm = T)
  y.av <- mean(y.d, na.rm = T)
  labels <- matrix(data = "Moderate", nrow = length(y))
  labels[y.d > (y.av + y.sd)] <- "High"
  labels[y.d < (y.av - y.sd)] <- "Low"
  my.plot(x, y, labels = labels, main = main, xlab = xlab, ylab = ylab)
  lines(v)
  return(y.d)
}

center.matrix <- function(m, dim = 1, sd.flag = F) {
  if (dim == 1) {
    zscores <- sweep(m, 1, rowMeans(m, na.rm = T), FUN = "-")
  } else {
    zscores <- sweep(m, 2, colMeans(m, na.rm = T), FUN = "-")
  }
  if (sd.flag) {
    zscores <- sweep(zscores, dim, apply(m, dim, function(x) (sd(x, na.rm = T))), FUN = "/")
  }
  return(zscores)
}

get.top.cor <- function(m, no.elm = 100, min.cf = 0, idx = NULL, add.prefix = "") {
  m <- as.matrix(m)
  m.pos <- (-m)
  m.neg <- m
  colnames(m.pos) <- paste0(colnames(m.pos), ".up")
  colnames(m.neg) <- paste0(colnames(m.neg), ".down")
  v <- get.top.elements(cbind(m.pos, m.neg),
    no.elm = no.elm,
    min.cf = (-abs(min.cf))
  )
  names(v) <- c(colnames(m.pos), colnames(m.neg))
  if (!is.null(idx)) {
    v <- v[paste(idx, c("up", "down"), sep = ".")]
  }
  names(v) <- paste0(add.prefix, names(v))
  return(v)
}

discretize.mat.q <- function(X, q1 = 0.9) {
  qv <- t(apply(X, 2, function(x) quantile(x, q1, na.rm = T)))
  B <- sweep(X, 2, qv, FUN = "-") >= 0
  return(B)
}

remove.redundant.dashs <- function(v) {
  v <- gsub("__", "_", v, fixed = T)
  v <- laply(v, function(x) {
    if (startsWith(x, "_")) {
      x <- substr(x, 2, nchar(x))
    }
    if (endsWith(x, "_")) {
      x <- substr(x, 1, nchar(x) - 1)
    }
    return(x)
  })
  return(v)
}

set.list <- function(r, b, data.name = NULL) {
  rn <- lapply(r, set.field, b = b)
  if (!is.null(data.name)) {
    rn$data.name <- data.name
  }
  return(rn)
}

set.field <- function(v, b) {
  d <- dim(v)
  d.b <- length(b)
  if (!is.null(d)) {
    if (d[1] == d.b) {
      v <- v[b, ]
    }
    if (d[2] == d.b) {
      v <- v[, b]
    }
  } else {
    if (length(v) == d.b) {
      v <- v[b]
    }
  }
  return(v)
}

cmp.multi.refs <- function(m, b, ref.groups, ranksum.flag = T) {
  # Previous name: ranksum.t.test.multi.refs
  ref.groups.u <- unique(ref.groups[!b])
  m.up <- get.mat(rownames(m), ref.groups.u)
  m.down <- m.up
  for (x in ref.groups.u) {
    b.ref <- is.element(ref.groups, x)
    if (length(unique(b.ref[b | b.ref])) < 2) {
      next()
    }
    m.r <- apply.ttest(m[, b | b.ref], b[b | b.ref], ranksum.flag = ranksum.flag)
    m.up[, x] <- m.r[, 1]
    m.down[, x] <- m.r[, 2]
    rm(m.r)
  }
  results <- list()
  results$up <- cbind.data.frame(
    max.p = rowMaxs(m.up, na.rm = T),
    No.sup = rowSums(m.up < 0.05, na.rm = T), m.up
  )
  results$down <- cbind.data.frame(
    max.p = rowMaxs(m.down, na.rm = T),
    No.sup = rowSums(m.down < 0.05, na.rm = T), m.down
  )
  return(results)
}

list.2.boolean.mat <- function(l, ids = NULL) {
  if (is.null(ids)) {
    ids <- sort(unique(unlist(l)))
  }
  B <- t(laply(l, function(x) is.element(ids, x)))
  colnames(B) <- names(l)
  rownames(B) <- ids
  return(B)
}

convert.to.mice <- function(mouse.human, genes) {
  genes <- as.character(unique(mouse.human[is.element(mouse.human[, "human"], genes), "mouse"]))
  return(genes)
}

GO.enrichment.lapply <- function(go.env, genes, sig, valuType = 1) {
  m <- t(laply(sig, function(x) GO.enrichment(go.env, genes, x)[, valuType]))
  colnames(m) <- names(sig)
  return(m)
}

GO.enrichment <- function(go.env, genes, selected.genes) {
  b <- is.element(genes, selected.genes)
  p <- laply(go.env, function(x) get.hyper.p.value(is.element(genes, x), b))
  rownames(p) <- names(go.env)
  return(p)
}

get.hyper.p.value <- function(b1, b2, full.flag = T) {
  b.na <- is.na(b1) | is.na(b2)
  b1 <- b1[!b.na]
  b2 <- b2[!b.na]
  p <- max(1 - phyper(sum(b1 & b2) - 1, sum(b1), sum(!b1), sum(b2)), 1e-17)
  if (!full.flag) {
    return(p)
  }
  ex <- sum(b2) * (sum(b1) / length(b2))
  ob.vs.ex <- sum(b1 & b2) / ex
  v <- c(p, ob.vs.ex, sum(b1 & b2), ex)
  names(v) <- c("hyper.p.value", "ob.vs.exp", "ob", "exp")
  return(v)
}

get.hyper.p.value.mat <- function(B1, B2) {
  B1 <- set.colnames(B1)
  B2 <- set.colnames(B2)
  P <- get.mat(colnames(B1), colnames(B2))
  for (i in 1:ncol(B1)) {
    for (j in 1:ncol(B2)) {
      P[i, j] <- get.hyper.p.value(B1[, i], B2[, j], full.flag = F)
    }
  }
  return(P)
}

get.residuals <- function(X, g) {
  set.seed(1234)
  f <- function(y) {
    lm(y ~ ., data = as.data.frame(g))$residuals
  }
  residuals <- t(apply(X, 1, f))
  return(residuals)
}

set.colnames <- function(m) {
  if (is.null(colnames(m))) {
    colnames(m) <- 1:ncol(m)
  }
  return(m)
}

cox.mat <- function(m, r, X = NULL, filter.flag = T) {
  if (is.null(X)) {
    f <- function(x) {
      summary(coxph(r$survival ~ x))$coefficients[1, c("coef", "Pr(>|z|)")]
    }
  } else {
    f <- function(x) {
      summary(coxph(r$survival ~ ., cbind.data.frame(x, X)))$coefficients[1, c("coef", "Pr(>|z|)")]
    }
  }

  if (filter.flag) {
    b <- rowSums(m > 0, na.rm = T) > (ncol(m) * 0.2)
  } else {
    b <- rowSums(m > 0, na.rm = T) > 0
  }
  p <- matrix(nrow = nrow(m), ncol = 2, dimnames = list(rownames(m)))
  p[b, ] <- t(apply(m[b, ], 1, f))
  p <- cbind(
    p, get.onesided.p.value(p[, 1], p[, 2]),
    get.onesided.p.value(-p[, 1], p[, 2])
  )
  p <- cbind(p, get.p.zscores(p[, 3:4]))
  colnames(p) <- c("coef", "Pr(>|z|)", "P.high.exp.worse", "P.high.exp.better", "Zscore")
  return(p)
}

get.auc.mat <- function(P, y) {
  auc <- apply(P, 2, function(x) get.auc(x, y))
  return(auc)
}

get.auc <- function(p1, y1) {
  pr <- prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  # prf <- performance(pr, measure = "prec", x.measure = "rec")
  # plot(prf,ylim = c(0,1))
  # abline(h = mean(y1))
  return(auc)
}

multi.gsub <- function(pattern, replacement = "", x) {
  for (i in 1:length(pattern)) {
    x <- gsub(pattern = pattern[i], replacement = replacement, x = x, fixed = T)
  }
  return(x)
}

format.pval.private <- function(p) {
  if (length(p) > 1) {
    P <- laply(p, format.pval.private)
    P <- gsub("P = ", "", P)
    P <- paste0("P", 1:length(P), " = ", P)
    P <- paste(P, collapse = ", ")
    return(P)
  }
  if (abs(p) > 1) {
    p <- 10^(-abs(p))
  }
  return(paste("P =", format(p, scientific = T, di = 3)))
}

boxplot.test <- function(vals, labels, las = 1, main = "", cex = 1.1, ylab = NULL, alternative = "two.sided",
                         t.test.flag = T, dots.flag = F, legend.flag = T, ref.label = labels[1]) {
  no.labels <- length(unique(labels))
  if (no.labels == 2) {
    col <- c("cadetblue", "gray")
  }
  if (no.labels == 3) {
    col <- c("cadetblue", "gray70", "brown4")
  }
  if (no.labels > 3) {
    col <- c("cadetblue", "lightblue", "gray70", "brown4", "orange")
  }
  if (t.test.flag == T) {
    p <- ttest.class(vals, is.element(labels, ref.label), alternative = alternative)
    if (alternative == "greater") {
      auc <- get.auc(vals, is.element(labels, ref.label))
    } else {
      auc <- get.auc(vals, !is.element(labels, ref.label))
    }
    auc <- round(auc, di = 2)
  } else {
    p <- wilcox.class(vals, is.element(labels, ref.label), alternative = alternative)
  }
  # main<-set.titles(gsub("_"," ",main))
  main <- gsub("_", " ", main)
  if (is.null(ylab)) {
    ylab <- main
  }
  if (t.test.flag != "none") {
    main <- gsub("_", " ", paste0(
      main, "\n(t-test P = ", format(p, scientific = T, di = 2),
      "\nAUC = ", auc, ")"
    ))
  }
  boxplot(vals ~ labels,
    las = las, cex.axis = 1, col = col, ylab = ylab, main = main, cex.main = cex,
    cex.lab = cex
  )
  l <- as.matrix(table(labels))
  l <- paste0(rownames(l), " (n = ", l, ")")
  if (legend.flag) {
    legend("topleft", legend = l, fill = col, cex = 0.8)
  }
  if (dots.flag) {
    stripchart(vals ~ labels,
      vertical = TRUE, method = "jitter",
      add = TRUE, pch = 20, col = "gray30"
    )
  }
  return(p)
}

km.plot3 <- function(r, v, main = "", X = NULL, qua = 0.2, xlim = NULL, direction = 0,
                     legend.flag = T, ylab = "Survival probability", four.levels = F) {
  M1 <- summary(coxph(r$survival ~ v))$coefficients
  coxD <- M1[1, "coef"]
  cox.p <- M1[1, "Pr(>|z|)"]
  if (!is.null(X)) {
    Mc <- summary(coxph(r$survival ~ cbind(v, X)))$coefficients
    cox.p.c <- Mc[1, "Pr(>|z|)"]
    coxD.c <- Mc[1, "coef"]
  }
  b1 <- v >= quantile(v, 1 - qua, na.rm = T)
  b2 <- v <= quantile(v, qua, na.rm = T)
  G <- ifelse(b1, "High", "Moderate")
  col <- c("red", "blue", "darkgreen")
  G[b2] <- "Low"
  km2 <- npsurv(r$survival ~ G)
  sdf2 <- survdiff(r$survival ~ G)
  sdf2 <- (1 - pchisq(sdf2$chisq, length(sdf2$n) - 1)) / 3
  l <- paste0(c("High", "Low", "Moderate"), " (", km2$n, ")")

  if (is.null(xlim)) {
    survplot(km2, col = col, lty = c(1, 1), xlab = "Years", label.curves = F, ylab = ylab, n.risk = T)
  } else {
    survplot(km2, col = col, lty = c(1, 1), xlab = "Years", label.curves = F, xlim = c(0, xlim), ylab = ylab, n.risk = T)
  }
  if (legend.flag) {
    legend("topright",
      fill = col[c(setdiff(1:length(col), 2), 2)], cex = 0.8,
      legend = l[c(setdiff(1:length(col), 2), 2)]
    )
  }

  if (!is.null(X)) {
    if (direction == 0) {
      P <- c(cox.p, cox.p.c, sdf2)
    } else {
      P <- get.onesided.p.value(direction * c(coxD, coxD.c, coxD), c(cox.p, cox.p.c, sdf2))
    }
    P <- format(P, scientific = T, di = 2)
    main <- paste0(main, "\nP=", P[1], ", Pc=", P[2], "\nlogrank=", P[3])
  } else {
    if (direction == 0) {
      P <- c(cox.p, sdf2)
    } else {
      P <- get.onesided.p.value(direction * c(coxD, coxD), c(cox.p, sdf2))
    }
    P <- format(P, scientific = T, di = 2)
    main <- paste0(main, "\nP=", P[1], ", logrank=", P[2])
  }
  title(main, cex.main = 1)
  return(G)
}

add.up.down.suffix <- function(v) {
  v <- lapply(v, function(x) x <- c(x, paste0(x, ".up"), paste0(x, ".down")))
  v <- unlist(v, use.names = F)
  return(v)
}

get.residuals <- function(X, g) {
  f <- function(y) {
    lm(y ~ ., data = as.data.frame(g))$residuals
  }
  residuals <- t(apply(X, 1, f))
  return(residuals)
}

get.anova.p <- function(y, x, order.flag = F) {
  b <- is.infinite(y) | is.na(y) | is.na(x)
  y <- y[!b]
  x <- x[!b]
  a <- aov(y ~ x)
  p <- unlist(summary(a))["Pr(>F)1"]
  return(p)
}

plot.bimodal.distribution <- function(y, density.flag = F, xlab = "", main = "",
                                      pos.label = "pos", neg.label = "neg") {
  set.seed(1234)
  mixmdl <- normalmixEM(y)
  main <- paste(main, paste("loglikelihood =", round(mixmdl$loglik)), sep = "\n")
  plot(mixmdl, which = 2, xlab2 = xlab, main2 = main)
  lines(density(y), lty = 2, lwd = 2)
  mixmdl$rank <- rank(y)
  idx1 <- order(mixmdl$mu)[1]
  mixmdl$labels <- mixmdl$rank > round(mixmdl$lambda[idx1] * length(y))
  if (density.flag) {
    my.densityplot(y, mixmdl$labels)
  }
  if (length(mixmdl$all.loglik) == 1001) {
    mixmdl$labels[] <- NA
  }
  mixmdl$labelsF <- ifelse(mixmdl$labels, pos.label, neg.label)
  y.pos <- y[mixmdl$labels]
  y.neg <- y[!mixmdl$labels]
  mixmdl$labelsF[(y < (mean(y.pos) - sd(y.pos))) & mixmdl$labels] <- paste0(pos.label, "?")
  mixmdl$labelsF[(y > (mean(y.neg) + sd(y.neg))) & !mixmdl$labels] <- paste0(neg.label, "?")
  return(mixmdl)
}


plot.extra <- function(x, y = NULL, labels, regression.flag = F, col.v = NULL, set.flag = F, cor.flag = F,
                       pch = 16, cex = 0.5, main = "", ylab = "tSNE2", xlab = "tSNE1", cex.axis = 0.6,
                       add.N = F) {
  print(xlab)
  main <- gsub(".", ": ", capitalize(main), fixed = T)
  if (add.N) {
    labels <- add.n.of.samples(labels)
  }
  if (set.flag) {
    par(mar = c(5.1, 4.1, 4.1, 12.1), xpd = TRUE)
  }
  if (is.null(col.v)) {
    col.v <- labels.2.colors(labels)
  }
  if (is.null(y)) {
    y <- x[, 2]
    x <- x[, 1]
  }
  if (cor.flag) {
    xy.cor <- get.cor(y, x)
    main <- paste(main, "\nR =", format(xy.cor[1], di = 2), "P =", format(xy.cor[2], scientific = T, di = 2))
  }
  plot(x, y, col = col.v, pch = pch, cex = cex, main = main, ylab = ylab, xlab = xlab, cex.axis = cex.axis)

  labels <- gsub(" ", "_", labels)
  l <- (max(x, na.rm = T) - min(x, na.rm = T)) / 20
  if (length(unique(labels)) < 30) {
    if (length(pch) == length(labels)) {
      map <- unique(paste(labels, col.v, pch))
      labels.n <- as.matrix(table(labels))
      idx <- match(get.strsplit(map, " ", 1), names(labels.n))
      map[, 1] <- paste0(map[, 1], " (N = ", m[idx], ")")
      legend(
        x = max(x, na.rm = T) + l, y = max(y, na.rm = T),
        legend = get.strsplit(map, " ", 1),
        col = get.strsplit(map, " ", 2), inset = c(-0.5, 0),
        pch = as.integer(get.strsplit(map, " ", 3)), lty = 1
      )
    } else {
      map <- unique(paste(labels, col.v, pch))
      legend(
        x = max(x, na.rm = T) + l, y = max(y, na.rm = T),
        legend = get.strsplit(map, " ", 1),
        col = get.strsplit(map, " ", 2),
        bty = "n", lty = 0, # line style
        lwd = 2, cex = 0.7, pch = 0
      )
    }
  }

  if (regression.flag == 1) {
    b <- !is.na(x) & !is.na(y)
    v <- lowess(x[b], y[b])
    lines(v)
  }
  if (regression.flag == 2) {
    b <- !is.na(x) & !is.na(y)
    ulabels <- unique(labels)
    for (i in ulabels) {
      bi <- b & labels == i
      v <- lowess(x[bi], y[bi])
      lines(v)
    }
  }
  return()
}

apply.plot.extra <- function(X, labels, main = NULL, xlab = "", ylab = "",
                             pch = 16, cex.axis = 0.6, set.flag = F) {
  laply(1:ncol(labels), function(i) {
    plot.extra(X,
      labels = labels[, i],
      main = ifelse(is.null(main), colnames(labels)[i],
        paste(main, colnames(labels)[i], sep = ":")
      ),
      xlab = xlab, ylab = ylab, pch = pch, cex.axis = cex.axis,
      set.flag = set.flag
    )
    return(i)
  })
}


add.n.of.samples <- function(l) {
  num.samples <- table(l)
  idx <- match(l, names(num.samples))
  l <- paste0(l, " (N=", num.samples[idx], ")")
  return(l)
}

labels.2.colors <- function(x.class, x = NULL, number.flag = F) {
  palette("default")
  col.v <- c("black", "red", setdiff(c(
    palette(), "cadetblue", "gray", "darkgreen", "darkorange", "darkviolet", "gold3",
    "lightpink", "deeppink2", "deepskyblue", rainbow(20)
  ), c("black", "red")))
  no.classes <- length(unique(x.class))
  if (number.flag) {
    col.v <- match(x.class, sort(unique(x.class)))
  } else {
    if (is.numeric(x.class[1]) && length(x.class) > 10) {
      col.v <- color.scale(x.class, c(0, 10), 0.8, 0.8,
        color.spec = "hsv",
        xrange = c(min(x.class, na.rm = T) - 1, max(x.class, na.rm = T) + 1)
      )
    } else {
      col.v <- col.v[match(x.class, sort(unique(x.class)))]
    }
  }
  if (!is.null(x)) {
    names(col.v) <- x
  }
  return(col.v)
}

get.strsplit <- function(v, sep, idx) {
  v <- as.character(v)
  vi <- laply(strsplit(v, split = sep, fixed = T), function(x) x[idx])
  return(vi)
}

gg.densityplot <- function(x, labels, title = "", subtitle = "", xlab = "Scores", legend.name = "", caption = "") {
  theme_set(theme_classic())
  b <- labels == labels[1]
  p <- t.test(x[b], x[!b])$p.value
  subtitle <- paste0(subtitle, " (", format.pval.private(p), ")")
  mpg <- cbind.data.frame(cty = x, cyl = labels)
  g <- ggplot(mpg, aes(cty))
  g <- g + geom_density(aes(fill = factor(cyl)), alpha = 0.5) +
    labs(
      title = title,
      x = xlab,
      subtitle = subtitle,
      fill = legend.name
    )
  # g <- g + theme_classic()
  return(g)
}

center.matrix <- function(m, dim = 1, sd.flag = F) {
  if (dim == 1) {
    zscores <- sweep(m, 1, rowMeans(m, na.rm = T), FUN = "-")
  } else {
    zscores <- sweep(m, 2, colMeans(m, na.rm = T), FUN = "-")
  }
  if (sd.flag) {
    zscores <- sweep(zscores, dim, apply(m, dim, function(x) (sd(x, na.rm = T))), FUN = "/")
  }
  return(zscores)
}

average.mat.rows <- function(m, ids, f = colMeans) {
  m <- cbind.data.frame(ids = ids, m)
  m.av <- ddply(.data = m, .variables = "ids", .fun = function(x) f(x[, -1]))
  rownames(m.av) <- m.av$ids
  m.av$ids <- NULL
  return(m.av)
}

get.auc <- function(p1, y1) {
  pr <- prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  return(auc)
}
