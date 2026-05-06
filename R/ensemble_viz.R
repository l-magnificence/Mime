#' Plot Kaplan-Meier survival curves for ensemble risk groups
#'
#' @param res.by.ML.Dev.Prog.Sig The results of ML.Dev.Prog.Sig (must contain Ensemble_riskscore)
#' @param list_train_vali_Data The list of cohort data used in ML.Dev.Prog.Sig
#' @param seed Random seed for reproducibility
#' @return A ggplot object (or list of ggplot objects, one per cohort)
#' @export
ensemble_km_plot <- function(res.by.ML.Dev.Prog.Sig,
                              list_train_vali_Data,
                              seed = 42) {
  library(survival)
  library(survminer)
  library(ggplot2)

  ensemble_rs <- res.by.ML.Dev.Prog.Sig$Ensemble_riskscore
  if (is.null(ensemble_rs)) {
    warning("No Ensemble_riskscore found in results. Run ML.Dev.Prog.Sig with ensemble_top_n > 0.")
    return(NULL)
  }

  plot_list <- list()
  for (cohort_name in names(ensemble_rs)) {
    df <- data.frame(
      ID = list_train_vali_Data[[cohort_name]]$ID,
      OS.time = as.numeric(list_train_vali_Data[[cohort_name]]$OS.time),
      OS = as.numeric(list_train_vali_Data[[cohort_name]]$OS),
      RS = ensemble_rs[[cohort_name]]
    )
    df$RiskGroup <- ifelse(df$RS >= median(df$RS), "High Risk", "Low Risk")

    fit <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = df)
    p <- survminer::ggsurvplot(
      fit, data = df,
      pval = TRUE, risk.table = TRUE,
      palette = c("#E7B800", "#2E9FDF"),
      title = paste0("Ensemble Model — ", cohort_name),
      xlab = "Time (Days)",
      ylab = "Overall Survival Probability"
    )
    plot_list[[cohort_name]] <- p
  }

  return(plot_list)
}

#' Plot C-index comparison bar chart including ensemble
#'
#' @param res.by.ML.Dev.Prog.Sig The results of ML.Dev.Prog.Sig
#' @param cohort_name Name of cohort to plot (default: first cohort)
#' @return A ggplot object
#' @export
ensemble_cindex_compare <- function(res.by.ML.Dev.Prog.Sig,
                                     cohort_name = NULL) {
  library(ggplot2)
  library(dplyr)

  result <- res.by.ML.Dev.Prog.Sig$Cindex.res
  if (is.null(cohort_name)) {
    cohort_name <- unique(result$ID)[1]
  }
  df <- result[result$ID == cohort_name, ]
  df <- df[order(-df$Cindex), ]

  # Highlight ensemble model
  df$IsEnsemble <- grepl("Ensemble", df$Model)

  p <- ggplot(df, aes(x = reorder(Model, Cindex), y = Cindex, fill = IsEnsemble)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "steelblue"), guide = "none") +
    labs(title = paste("C-index Comparison —", cohort_name), x = "Model", y = "C-index") +
    theme_minimal()

  return(p)
}
