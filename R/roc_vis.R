#' Title
#'
#' @param object
#' @param model_name
#' @param dataset
#' @param order
#' @param dataset_col
#' @param anno_position
#' @param year
#'
#' @return
#' @export
#'
#' @examples
roc_vis <- function(object, # output of cal_AUC_ml_res  mode = "all",'single' or 'double'
                    model_name, # input one identical model name in output of cal_AUC_ml_res
                    dataset, # input identical name of cohort in output of cal_AUC_ml_res
                    order = NULL, # cohort order plot and input identical name of cohort in output of cal_AUC_ml_res
                    dataset_col = NULL, # color value for dataset
                    anno_position = NULL, # auc value position
                    year # year auc
) {
  library(ggplot2)

  if (is.null(anno_position) == T) {
    anno_position <- c(0.53, 0.35)
  } else {
    anno_position <- anno_position
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

  tmp <- data.frame()
  for (i in 1:length(dataset)) {
    tmp2 <- object[[model_name]][[dataset[i]]]
    tmp2$dataset <- dataset[i]
    tmp <- rbind(tmp, tmp2)
  }

  if (is.null(order) == T) {
    order <- dataset
    tmp$dataset <- factor(tmp$dataset, levels = order) ## default order
  } else {
    order <- order
    tmp$dataset <- factor(tmp$dataset, levels = order)
  }

  p1 <- ggplot(tmp, aes(FP, TP, color = dataset)) +
    geom_line(size = 0.5, alpha = 1) +
    scale_color_manual(values = dataset_col, name = "Cohort") +
    # geom_area(aes(fill = dataset),alpha = 0.5,stat = "align",position = "identity")+
    labs(
      title = paste(year, "Year AUC predicted by", model_name),
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      plot.title = element_text(hjust = 0.5),
      legend.position = "",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    scale_x_continuous(expand = c(0.01, 0)) + # right most position will be 1 + (1-0) * 0.01 + 0 = 10
    scale_y_continuous(expand = c(0.01, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 0.5)

  for (i in 1:length(dataset)) {
    ano_y <- anno_position[2] - (i - 1) * 0.05
    p1 <- p1 +
      annotate("segment",
               x = anno_position[1], xend = anno_position[1] + 0.05, y = ano_y, yend = ano_y,
               colour = dataset_col[i]
      ) +
      annotate("text",
               x = anno_position[1] + 0.07, y = ano_y, hjust = 0,
               label = paste("AUC in", order[i], ":", round(unique(tmp[tmp$dataset == order[i], "AUC"]), 3))
      )
  }
  print(p1)
}
