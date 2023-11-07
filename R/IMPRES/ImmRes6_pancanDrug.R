find.sensitizing.drugs <- function(r = NULL) {
  if (is.null(r)) {
    r <- readRDS("../Data/PublicData/Sanger.Garnett.data.rds")
  }
  if (is.null(r$res)) {
    r$res <- cmb.res.scores(r = r, res.sig = res.sig, bulk.flag = T)
  }
  results <- list(
    garnett = apply.drug.response.vs.state.HLM(ic50 = r$garnett.ic50, state = r$res[, "res"], r = r),
    gong = cdk46i.efficacy.analysis(
      r = r, state = r$res[, "res"],
      q.drug.sen = 0.25, q.drug.res = 0.5, q.res.state = 0.75
    )
  )
  print("Garnett et al. 2012, top 5 immune sensitizing drugs:")
  print(results$garnett[1:5, ])
  print("Gong et al. 2017, CDK4/6i is more effective in immune resistant cell lines")
  print(format.pval.private(results$gong["discrete", "P"]))
  saveRDS(results, file = "../Results/PanCancerDrugs.rds")
  return()
}

apply.drug.response.vs.state.HLM <- function(ic50, state, r, disc.flag = T, q.drug.sen = 0.25, q.res.state = 0.75) {
  r$state <- state
  p <- t(apply(ic50, 2, function(x) {
    drug.response.vs.state.HLM(x, r,
      disc.flag = disc.flag,
      q.drug.sen = q.drug.sen,
      q.res.state = q.res.state
    )
  }))
  p <- cbind.data.frame(Z = get.cor.zscores(p[, "Estimate"], p[, "Pr(>|z|)"]), p)
  p <- p[order(p[, "Z"], decreasing = T), ]
  return(p)
}

drug.response.vs.state.HLM <- function(y, r0, disc.flag = T, q.drug.sen = 0.25, q.res.state = 0.75) {
  r0$y <- y
  b <- !is.na(y)
  r0 <- list(state = r0$state[b], y = r0$y[b], types = r0$types[b])
  r0$y <- (r0$y < quantile(r0$y, q.drug.sen))
  if (disc.flag) {
    r0$state <- r0$state >= quantile(r0$state, q.res.state)
  }
  f <- function(r0) {
    M1 <- with(r0, glmer(y ~ state + (1 | types),
      family = binomial(link = "logit"),
      control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 20000))
    ))
    m.coef <- summary(M1)$coef
    idx <- match(c("Estimate", "Pr(>|z|)"), colnames(m.coef))
    if (!is.element("state", rownames(m.coef))) {
      c1 <- m.coef["stateTRUE", idx]
    } else {
      c1 <- m.coef["state", idx]
    }
    if (!is.null(summary(M1)$optinfo$conv$lme4$messages)) {
      c1[] <- NA
    }
    # print(c1)
    return(c1)
  }
  c1 <- tryCatch(
    {
      f(r0)
    },
    error = function(err) {
      return(c(NA, NA))
    }
  )
  return(c1)
}

cdk46i.efficacy.analysis <- function(r, state, q.drug.sen = 0.1, q.drug.res = 0.5, q.res.state = 0.9) {
  # Here we will use the efficacies of two CDK4/6 inhibitors (palbo. and abem.) to determine if a given
  # cell line is sensitive to CDK4/6 inhibition
  f <- function(x) {
    v <- ifelse(x < quantile(x, q.drug.sen, na.rm = T), "sen",
      ifelse(x > quantile(x, q.drug.res, na.rm = T), "res", "?")
    )
    return(v)
  }
  r$state <- state
  l <- cbind.data.frame(
    palbo = paste0("palbo.", f(r$cdk46.ic50$palbo)),
    abema = paste0("abem.", f(r$cdk46.ic50$abema))
  )
  l$final <- paste(l[, 1], l[, 2], sep = "_")
  l$final[l$final == "palbo.sen_abem.sen"] <- "CDK4/6i sensitive"
  l$final[l$final == "palbo.res_abem.res"] <- "CDK4/6i resistant"
  b <- cbind(!is.na(r$cdk46.ic50), grepl("CDK", l$final))
  r$cdk <- l$final
  r0 <- set.list(r, b[, 3])
  r0$stateB <- r0$state > quantile(r0$state, q.res.state)
  r0$y <- grepl("sen", r0$cdk)
  M1 <- with(r0, glmer(y ~ stateB + (1 | types),
    family = binomial(link = "logit"),
    control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 20000))
  ))
  M2 <- with(r0, glmer(y ~ state + (1 | types),
    family = binomial(link = "logit"),
    control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 20000))
  ))
  Z <- rbind(
    discrete = summary(M1)$coef["stateBTRUE", c("Estimate", "Pr(>|z|)")],
    cont = summary(M2)$coef["state", c("Estimate", "Pr(>|z|)")]
  )
  Z <- cbind(Z, P = get.onesided.p.value(Z[, "Estimate"], Z[, "Pr(>|z|)"]))
  return(Z)
}
