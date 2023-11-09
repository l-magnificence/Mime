#' Title
#'
#' @param obj
#' @param dataset
#' @param color
#' @param feature1
#' @param feature2
#' @param method
#'
#' @return
#' @export
#'
#' @examples
cor_plot <- function(obj, # expression data frame
                     dataset, # the cohort name of InputMatrix
                     color = NULL, # color value
                     feature1, # signature selected
                     feature2,
                     method) { # pearson or spearman
  library(ggplot2)
  library(ggpubr)

  if (is.null(color) == T) {
    color <- "#1a9641" ## default color value
  } else {
    color <- color
  }

  scatter <- ggplot(obj, aes(obj[, feature1], obj[, feature2])) +
    geom_point(color = "grey", alpha = 0.8) +
    geom_smooth(method = "lm", color = color) +
    stat_cor(method = method, aes(x = obj[, feature1], y = obj[, feature2], color = "#fdae61")) +
    theme_bw() +
    theme(
      legend.position = "none", plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    labs(x = feature1, y = feature2, title = dataset)

  return(scatter)
}
