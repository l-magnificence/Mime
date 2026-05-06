#' Compute extra evaluation metrics for prognosis models
#'
#' @param res.by.ML.Dev.Prog.Sig The results of function ML.Dev.Prog.Sig
#' @param train_data The training data
#' @param inputmatrix.list A list of validation/training data frames
#' @param auc_time_vec Numeric vector of time points in years for AUC (default c(1,3,5))
#' @return A list with IBS, UnoC, Dindex results per model
#' @export
cal_extra_metrics <- function(res.by.ML.Dev.Prog.Sig,
                               train_data,
                               inputmatrix.list,
                               auc_time_vec = c(1, 3, 5)) {
  library(survival)

  result <- list(IBS = list(), UnoC = list(), Dindex = list())

  sig_genes <- res.by.ML.Dev.Prog.Sig$Sig.genes
  ml_names <- names(res.by.ML.Dev.Prog.Sig$ml.res)
  cohort_names <- names(inputmatrix.list)

  # Preprocess input data (same as other downstream functions)
  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    x
  })
  colnames(train_data) <- gsub("-", ".", colnames(train_data))

  for (model_name in ml_names) {
    if (is.null(res.by.ML.Dev.Prog.Sig$ml.res[[model_name]])) next

    # Integrated Brier Score
    result$IBS[[model_name]] <- tryCatch({
      library(pec)
      # Build Cox model from risk scores for pec input
      rs_list <- res.by.ML.Dev.Prog.Sig$riskscore[[model_name]]
      if (!is.null(rs_list) && !is.null(rs_list[[1]])) {
        train_rs <- rs_list[[1]]  # first cohort = training
        df <- data.frame(
          OS.time = as.numeric(train_rs$OS.time),
          OS = as.numeric(train_rs$OS),
          RS = train_rs$RS
        )
        cox_fit <- coxph(Surv(OS.time, OS) ~ RS, data = df)
        pec_res <- pec::pec(
          object = list("Mime_Model" = cox_fit),
          formula = Surv(OS.time, OS) ~ 1,
          data = df,
          exact = FALSE,
          verbose = FALSE
        )
        list(IBS = pec_res$AppErr$Mime_Model, time = pec_res$time)
      } else {
        NULL
      }
    }, error = function(e) {
      message(paste0("[SKIP] IBS for ", model_name, ": ", e$message))
      NULL
    })

    # Uno's C-statistic
    result$UnoC[[model_name]] <- tryCatch({
      library(survC1)
      rs_list <- res.by.ML.Dev.Prog.Sig$riskscore[[model_name]]
      if (!is.null(rs_list) && !is.null(rs_list[[1]])) {
        train_rs <- rs_list[[1]]
        df <- data.frame(
          time = as.numeric(train_rs$OS.time),
          status = as.numeric(train_rs$OS),
          marker = train_rs$RS
        )
        tau <- max(df$time) * 0.9
        survC1::Est.Cval(df, tau = tau, nofit = TRUE)
      } else {
        NULL
      }
    }, error = function(e) {
      message(paste0("[SKIP] UnoC for ", model_name, ": ", e$message))
      NULL
    })

    # D-index
    result$Dindex[[model_name]] <- tryCatch({
      rs_list <- res.by.ML.Dev.Prog.Sig$riskscore[[model_name]]
      if (!is.null(rs_list) && !is.null(rs_list[[1]])) {
        train_rs <- rs_list[[1]]
        df <- data.frame(
          OS.time = as.numeric(train_rs$OS.time),
          OS = as.numeric(train_rs$OS),
          RS = train_rs$RS
        )
        cox_fit <- coxph(Surv(OS.time, OS) ~ RS, data = df)
        # D-index = exp(coef) from Cox model
        list(
          Dindex = as.numeric(exp(coef(cox_fit))),
          lower = as.numeric(exp(confint(cox_fit))[1]),
          upper = as.numeric(exp(confint(cox_fit))[2]),
          pvalue = summary(cox_fit)$coefficients[, "Pr(>|z|)"]
        )
      } else {
        NULL
      }
    }, error = function(e) {
      message(paste0("[SKIP] Dindex for ", model_name, ": ", e$message))
      NULL
    })
  }

  return(result)
}
