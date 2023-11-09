#' Meta analysis of the result from univariate cox regression
#'
#' Creates a summary table of meta analysis results 
#'
#' @param object Output of function cal_unicox_meta_ml_res
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for datasets
#' @param dataset A vector of names for all datasets 
#'
#' @return a summary table
#' @export
#'
#' @examples
#' meta_unicox_vis(metamodel,
#'                 dataset = names(list_train_vali_Data))
#' 
meta_unicox_vis <- function(object, # output of cal_unicox_meta_ml_res
                            dataset_col = NULL, # color value for cohorts
                            dataset # input datasets name
) {
  library(forestploter)
  library(grid)

  if (is.null(dataset_col) == T) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    dataset_col <- dataset_col
  }

  metamodel <- object
  dt <- metamodel[["data"]]
  dt$Weight_random <- paste(round(100 * metamodel$w.random / sum(metamodel$w.random), 2), "%", sep = "")
  dt$Weight_fixed <- paste(round(100 * metamodel$w.fixed / sum(metamodel$w.fixed), 2), "%", sep = "")

  dt[nrow(dt) + 1, ] <- NA
  dt[nrow(dt), 6] <- "Random effect model"
  dt[nrow(dt), 1] <- exp(metamodel$TE.random)
  dt[nrow(dt), 3] <- exp(metamodel$lower.random)
  dt[nrow(dt), 4] <- exp(metamodel$upper.random)
  dt[nrow(dt), 2] <- metamodel$pval.random
  dt[nrow(dt), 10] <- "100%"
  dt[nrow(dt), 11] <- "--"

  dt[nrow(dt) + 1, ] <- NA
  dt[nrow(dt), 6] <- "Fixed effect model"
  dt[nrow(dt), 1] <- exp(metamodel$TE.fixed)
  dt[nrow(dt), 3] <- exp(metamodel$lower.fixed)
  dt[nrow(dt), 4] <- exp(metamodel$upper.fixed)
  dt[nrow(dt), 2] <- metamodel$pval.fixed
  dt[nrow(dt), 10] <- "--"
  dt[nrow(dt), 11] <- "100%"

  dt$se <- (log(dt$HCI) - log(dt$HR)) / 1.96
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  dt$`HR (95% CI)` <- ifelse(is.na(dt$Group), "",
                             sprintf(
                               "%.2f (%.2f - %.2f)",
                               dt$HR, dt$LCI, dt$HCI
                             )
  )
  dt$p <- ifelse(dt$pvalue < 0.001, "P<0.001", sprintf("%.3f", dt$pvalue))
  colnames(dt)[c(6:8, 10:11, 15)] <- c("Cohorts", "TE", "SE(TE)", "Weight(random)", "Weight(fixed)", "P value")
  dt$TE <- round(dt$TE, 2)
  dt$`SE(TE)` <- round(dt$`SE(TE)`, 2)
  dt[c((nrow(dt) - 1):nrow(dt)), c(7, 8)] <- ""
  rownames(dt) <- dt$Cohorts
  dt <- dt[c(dataset, "Random effect model", "Fixed effect model"), ]

  tm <- forest_theme(
    core = list(bg_params = list(
      fill = c(dataset_col[1:length(dataset)], "grey", "grey"),
      alpha = 0.5
    )),
    base_size = 10,
    # Confidence interval point shape, line type/color/width
    ci_pch = 16,
    ci_col = "#762a83",
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = 0.2, # Set an T end at the end of CI
    # Reference line width/type/color
    refline_lwd = 1,
    refline_lty = "dashed",
    refline_col = "grey20",
    # Vertical line width/type/color
    vertline_lwd = 1,
    vertline_lty = "dashed",
    vertline_col = "grey20",
    # Change summary color for filling and borders
    summary_fill = "#4575b4",
    summary_col = "#4575b4",
    # Footnote font size/face/color
    footnote_cex = 1,
    footnote_fontface = "italic",
    footnote_col = "red"
  )

  p <- forestploter::forest(dt[, c(6:8, 10:11, 13, 14, 15)],
                            est = dt$HR,
                            lower = dt$LCI,
                            upper = dt$HCI,
                            sizes = dt$se,
                            is_summary = c(rep(FALSE, nrow(dt) - 2), TRUE, TRUE),
                            ci_column = 6,
                            ref_line = 1,
                            arrow_lab = c("Better", "Worse"),
                            x_trans = "log2",
                            xlim = c(0, ceiling(max(dt$HCI))),
                            ticks_at = c(0.5, 2^seq(0, floor(log2((ceiling(max(dt$HCI)) - 1))), by = 1)),
                            footnote = " Univariate Cox Regression",
                            theme = tm
  )

  p <- add_text(p,
                text = "Meta analysis of univariate Cox regression",
                part = "header",
                row = 0,
                col = 4:6,
                just = c("center"),
                gp = gpar(fontface = "bold")
  )

  p <- add_border(p,
                  part = "header",
                  row = c(0, 1),
                  gp = gpar(lwd = 1)
  )

  p <- insert_text(p,
                   text = "Meta analysis",
                   row = length(dataset) + 1,
                   just = "left",
                   gp = gpar(cex = 0.8, col = "blue", fontface = "italic")
  )
  print(p)
}
