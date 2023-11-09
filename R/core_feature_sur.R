#' Survival curve of patients based on specific gene
#'
#' Creates a survival curve of patients according to expression of specific gene among different datasets
#'
#' @param gene Specific gene name
#' @param InputMatrix A matrix including OS, OS.time and gene expression
#' @param dataset Identical name of dataset
#' @param cutoff Cutoff for patient's group, which can be mean or quantile value
#' @param conf.int If TRUE, plot with confidence intervals otherwise FALSE 
#' @param median.line Drawing a horizontal/vertical line at median survival, it can be c("none", "hv", "h", "v"). If NULL, plot without line (default)
#' @param color If NULL, color values are set to the default colors. Otherwise, you can specify two color values for high and low group
#' @param xlab X axis title
#' @param pval.coord The position of p value
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' survplot <- vector("list",2) 
#' for (i in c(1:2)) {
#'   print(survplot[[i]]<-core_feature_sur("PSEN2", 
#'                                         InputMatrix=list_train_vali_Data[[i]],
#'                                         dataset = names(list_train_vali_Data)[i],
#'                                         #color=c("blue","green"),
#'                                         median.line = "hv",
#'                                         cutoff = 0.5,
#'                                         conf.int = T,
#'                                         xlab="Day",pval.coord=c(1000,0.9)))
#' }
#' aplot::plot_list(gglist=survplot,ncol=2)
#' 
core_feature_sur <- function(gene, # gene name
                             InputMatrix, # InputMatrix with os os.time and gene expression need contain gene
                             dataset, # the cohort name of InputMatrix
                             cutoff = NULL, # cutoff for group: mean or quantile value
                             conf.int = F, # show confidence intervals T or F
                             median.line = NULL, # drawing a horizontal/vertical line at median survival. c("none", "hv", "h", "v")
                             color = NULL, # two color value for high and low group
                             xlab = NULL, # x axis title
                             pval.coord = NULL # p value position
) {
  library(survival)
  library(survminer)

  if (is.null(color) == T) {
    color <- c("#868686", "#B24745") ## default color value
  } else {
    color <- color
  }

  if (is.null(median.line) == T) {
    median.line <- "none" ## default not show
  } else {
    median.line <- median.line
  }

  if (is.null(xlab) == T) {
    xlab <- "Time"
  } else {
    xlab <- xlab
  }

  if (is.null(pval.coord) == T) {
    pval.coord <- c(1, 0.25)
  } else {
    pval.coord <- pval.coord
  }


  tmp <- InputMatrix[, c("OS.time", "OS", gene)]
  colnames(tmp)[3] <- "gene"
  tmp <- tmp[!is.na(tmp$gene), ]
  mySurv <- Surv(tmp$OS.time, tmp$OS)

  if (is.null(cutoff) == T) {
    cutoff <- 0.5
    value <- quantile(tmp$gene, probs = c(cutoff))
  } else {
    if (cutoff == "mean") {
      value <- mean(tmp$gene)
    } else {
      cutoff <- cutoff
      value <- quantile(tmp$gene, probs = c(cutoff))
    }
  }

  tmp$Group <- ifelse(tmp$gene > value, "High", "Low")
  Group <- tmp$Group
  Group <- factor(Group, levels = c("Low", "High"))
  tmp$Group <- factor(Group, levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = tmp)

  # calculate HR and 95%CI
  data.survdiff <- survdiff(mySurv ~ Group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")

  p1 <- ggsurvplot(fit,
                   data = tmp,
                   conf.int = conf.int,
                   censor = F,
                   palette = color,
                   legend.title = paste0("Expression in ", dataset),
                   font.legend = 11,
                   surv.median.line = median.line,
                   legend.labs = c(
                     paste0("Low", "(n=", fit$n[1], ")"),
                     paste0("High", "(n=", fit$n[2], ")")
                   ),
                   pval = paste(
                     pval = ifelse(p.val < 0.001, "p < 0.001",
                                   paste("p = ", round(p.val, 3), sep = "")
                     ),
                     HR, CI, sep = "\n"
                   ),
                   xlab = xlab,
                   title = gene,
                   pval.coord = pval.coord
  )
  p1 <- p1$plot + theme(
    plot.title = element_text(hjust = 0.5)
  )

  print(p1)
}
