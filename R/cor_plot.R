#' Correlation between two genes
#'
#' Creates a scatter plot between two genes
#'
#' @param obj A data frame of gene expression
#' @param dataset Identical name of dataset
#' @param color If NULL, color values are set to the default colors. Otherwise, you can specify one color value for line
#' @param feature1 A specific gene
#' @param feature2 Another specific gene
#' @param method Methods for correlation analysis, it can be "pearson" or "spearman"
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' dataset_col<-c("#3182BDFF","#E6550DFF")
#' corplot <- list()
#' for (i in c(1:2)) {
#'   print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
#'                                dataset=names(list_train_vali_Data)[i],
#'                                color = dataset_col[i],
#'                                feature1="PSEN2",
#'                                feature2="WNT5B",
#'                                method="pearson"))
#' }
#' aplot::plot_list(gglist=corplot,ncol=2)
#' 
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
