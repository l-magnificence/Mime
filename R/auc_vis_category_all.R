#' Area Under Curve distribution plot when predicting response
#'
#' Creates a distribution plot of AUC among different methods when predicting response
#'
#' @param object Output of function ML.Dev.Pred.Category.Sig
#' @param dataset A vector of names for all datasets in object
#' @param order If NULL, order is set to the default order. Otherwise, you can specify the order of datasets to plot 
#' @param method_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for methods in object
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' auc_vis_category_all(res.ici,dataset = c("training","validation"),
#'                      order= c("training","validation"))
#' 
auc_vis_category_all <- function(object, # output of ML.Dev.Pred.Category.Sig
                                 dataset, # input identical name of cohort in output of ML.Dev.Pred.Category.Sig
                                 order = NULL, # cohort order plot and input identical name of cohort in output of ML.Dev.Pred.Category.Sig
                                 method_col = NULL # color value for method
) {
  library(ggplot2)

  if (is.null(method_col) == T) {
    method_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    method_col <- method_col
  }

  tmp <- data.frame()
  for (i in 1:length(dataset)) {
    auc <- object[["auc"]][[dataset[i]]]
    auc$method <- rownames(auc)
    auc$dataset <- dataset[i]
    rownames(auc) <- NULL
    tmp <- rbind(tmp, auc)
  }
  tmp$ROC <- as.numeric(sprintf("%.2f", tmp$ROC))
  tmp <- tmp[order(tmp$ROC, decreasing = T), ]
  tmp$method <- factor(tmp$method, levels = unique(tmp$method))

  if (is.null(order) == T) {
    order <- dataset
    tmp$dataset <- factor(tmp$dataset, levels = order) ## default order
  } else {
    order <- order
    tmp$dataset <- factor(tmp$dataset, levels = order)
  }

  p1 <- ggplot(data = tmp, aes(x = dataset, y = ROC, fill = method)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = method_col, name = "Method") +
    geom_text(aes(label = ROC), position = position_dodge(width = 0.9), vjust = 1.5, hjust = 0.5, size = 3) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      # legend.position = "",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    ) +
    # scale_y_continuous(position = "right")+
    labs(y = "AUC", x = "", title = "")
  print(p1)
}
