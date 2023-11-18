#' Select critical genes from an gene list
#'
#' Creates an upset plot of genes filtered by different methods 
#'
#' @param object Output of function ML.Corefeature.Prog.Screen whose mode can be "all" or "all_without_SVM"
#' @param sets A vector of names for methods in object. If NULL, show all methods (default)
#' @param mb.ratio Ratio between matrix plot and main bar plot. If NULL, ratio is set to c(0.6, 0.4) (default)
#' @param col Color values for sets.bar.color,main.bar.color,matrix.color, shade.color. If NULL, color values are set to the default colors
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' core_feature_select(res.feature.all)
#' 
core_feature_select <- function(object, # output of ML.Corefeature.Prog.Screen, mode= all or all_without_SVM
                                sets = NULL, # names of method
                                mb.ratio = NULL, # Ratio between matrix plot and main bar plot
                                col = NULL # color value for sets.bar.color,main.bar.color,matrix.color, shade.color
) {
  library(UpSetR)

  if (is.null(col) == T) {
    col <- c("#E18727", "#B09C85", "#ADB6B6", "#B09C85") ## default 4 color values
  } else {
    col <- col
  }

  if (is.null(mb.ratio) == T) {
    mb.ratio <- c(0.6, 0.4) ## default ratio
  } else {
    mb.ratio <- mb.ratio
  }

  core_feature_list <- list()
  for (i in unique(object$method)) {
    core_feature_list[[i]] <- object[object$method == i, "selected.fea"]
  }

  if (is.null(sets) == T) {
    sets <- names(core_feature_list) ## default order
  } else {
    sets <- sets
  }

  p1 <- upset(fromList(core_feature_list),
              # nsets=2,
              sets = sets,
              order.by = "freq", nintersects = NA,
              mb.ratio = mb.ratio,
              keep.order = T,
              mainbar.y.label = "Shared gene number",
              sets.x.label = "Total gene number",
              point.size = 2,
              line.size = 1,
              sets.bar.color = col[1],
              main.bar.color = col[2],
              matrix.color = col[3],
              shade.color = col[4]
  )
  print(p1)
}
