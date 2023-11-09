#' Title
#'
#' @param object
#' @param object2
#' @param color
#' @param model_name
#' @param dataset_col
#' @param dataset
#' @param type
#'
#' @return
#' @export
#'
#' @examples
HR_com <- function(object, # output of cal_auc_pre.prog.sig
                   object2, # output of cal_AUC_ml_res
                   color = NULL, # three color value for HR
                   model_name, ## if input object2 is from mode=all, then define model as specific model name while SOD
                   dataset_col = NULL, # color value for cohort
                   dataset, # input datasets name
                   type # c('categorical'),c('continuous) Univariate cox regression analysis, suggest categorical group by median
) {
  library(survival)
  library(survminer)
  library(ggplot2)

  if (is.null(color) == T) {
    color <- c("#0084A7", "white", "#E05D00") ## default color value
  } else {
    color <- color
  }

  if (is.null(dataset_col) == T) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    dataset_col <- dataset_col
  }

  hr <- data.frame()
  for (i in names(object)) {
    for (n in dataset) {
      rs <- object[[i]][[n]]

      if (type == "categorical") {
        # OS
        group <- ifelse(rs$RS > median(rs$RS), "high", "low")

        if (length(unique(group)) == 1) {
          p.val <- 1
          HR <- 1

          hr <- rbind(hr, data.frame(
            "P" = p.val,
            "HR" = HR,
            "ID" = n,
            "Signature" = i
          ))
        } else {
          group <- factor(group, levels = c("low", "high"))

          my.surv <- Surv(rs$OS.time, rs$OS)
          univ_cox <- coxph(my.surv ~ group, data = rs)
          univ_coxSummary <- summary(univ_cox)
          p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
          HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

          hr <- rbind(hr, data.frame(
            "P" = p.val,
            "HR" = HR,
            "ID" = n,
            "Signature" = i
          ))
        }
      } else {
        my.surv <- Surv(rs$OS.time, rs$OS)
        univ_cox <- coxph(my.surv ~ RS, data = rs)
        univ_coxSummary <- summary(univ_cox)
        p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
        HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

        hr <- rbind(hr, data.frame(
          "P" = p.val,
          "HR" = HR,
          "ID" = n,
          "Signature" = i
        ))
      }
    }
  }
  hr$text <- "black"

  hr2 <- data.frame()
  if (model_name != "SOD") {
    object3 <- object2[["riskscore"]][[model_name]]
  } else {
    object3 <- object2[["riskscore"]][[1]]
  }
  for (n in dataset) {
    rs <- object3[[n]]

    if (type == "categorical") {
      # OS
      group <- ifelse(rs$RS > median(rs$RS), "high", "low")

      if (length(unique(group)) == 1) {
        p.val <- 1
        HR <- 1

        hr2 <- rbind(hr2, data.frame(
          "P" = p.val,
          "HR" = HR,
          "ID" = n,
          "Signature" = "Our model"
        ))
      } else {
        group <- factor(group, levels = c("low", "high"))

        my.surv <- Surv(rs$OS.time, rs$OS)
        univ_cox <- coxph(my.surv ~ group, data = rs)
        univ_coxSummary <- summary(univ_cox)
        p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
        HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

        hr2 <- rbind(hr2, data.frame(
          "P" = p.val,
          "HR" = HR,
          "ID" = n,
          "Signature" = "Our model"
        ))
      }
    } else {
      my.surv <- Surv(rs$OS.time, rs$OS)
      univ_cox <- coxph(my.surv ~ RS, data = rs)
      univ_coxSummary <- summary(univ_cox)
      p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
      HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

      hr2 <- rbind(hr2, data.frame(
        "P" = p.val,
        "HR" = HR,
        "ID" = n,
        "Signature" = "Our model"
      ))
    }
  }
  hr2$text <- "red"

  hr <- rbind(hr, hr2)

  hr$pstar <- ifelse(hr$P < 0.05,
                     ifelse(hr$P < 0.01, "**", "*"),
                     ""
  )
  tmp <- hr
  tmp$Signature <- factor(hr$Signature, levels = c("Our model", unique(names(object))))
  tmp$ID <- factor(tmp$ID, levels = dataset)

  p <- ggplot(tmp, aes(Signature, ID)) +
    geom_tile(aes(fill = HR), colour = "white", size = 1) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = 1) +
    geom_text(aes(label = pstar), col = "black", size = 5) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(colour = dataset_col),
      axis.text.x = element_text(colour = c("red", rep("black", length(names(object)))), angle = 45, hjust = 1)
    ) +
    labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "HR"))
  print(p)
}
