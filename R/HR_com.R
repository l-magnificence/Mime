#' Compare the HR of our model with previously published models
#'
#' Creates a distribution plot of hazard ratio among different mdoels 
#'
#' @param object Output of function cal_auc_pre.prog.sig
#' @param object2 Output of function cal_AUC_ml_res
#' @param color If NULL, color values are set to the default colors. Otherwise, you can specify three color value for HR to plot
#' @param model_name Model name in object2 used to compare. If object2 is from all mode in cal_AUC_ml_res, then define model_name as specific model name while define model_name as "SOD"
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for datasets
#' @param dataset A vector of names for all datasets 
#' @param type Input "categorical" or "continuous" used by Univariate cox regression analysis, suggesting categorical which group by median
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' HR_com(rs.glioma.lgg.gbm,
#'     res,
#'     model_name="RSF + survival-SVM",
#'     dataset=names(list_train_vali_Data),
#'     type = "categorical")
#'
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
