#' Compare the c-index of our model with previously published models
#'
#' Creates a distribution plot of C-index among different mdoels 
#'
#' @param object Output of function cal_RS_pre.prog.sig
#' @param object2 Output of function ML.Dev.Prog.Sig
#' @param model_name Model name in object2 used to compare. If object2 is from all mode in ML.Dev.Prog.Sig, then define model_name as specific model name while define model_name as "SOD"
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for datasets
#' @param dataset A vector of names for all datasets 
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' cindex_comp(cc.glioma.lgg.gbm,
#'             res,
#'             model_name="StepCox[both] + plsRcox",
#'             dataset=names(list_train_vali_Data))
#' 
cindex_comp <- function(object, # output of cal_RS_pre.prog.sig
                        object2, # output of ML.Dev.Prog.Sig
                        model_name, ## if input object2 is from mode=all, then define model as specific model name while SOD
                        dataset_col = NULL, # color value for cohort
                        dataset # input datasets name
) {
  library(ggplot2)
  library(aplot)

  if (is.null(dataset_col) == T) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    dataset_col <- dataset_col
  }

  cindex_comp_other <- data.frame()
  for (i in names(object)) {
    tmp <- object[[i]]
    tmp$Signature <- i
    cindex_comp_other <- rbind(cindex_comp_other, tmp)
  }

  cindex_comp_other$text <- "black"

  if (model_name != "SOD") {
    cindex <- object2[["Cindex.res"]]
    cindex <- cindex[cindex$Model == model_name, ]
    colnames(cindex)[3] <- "Signature"
    cindex$Signature <- "Our model"
  } else {
    cindex <- object2[["Cindex.res"]]
  }
  cindex$text <- "red"

  cindex_comp_other <- rbind(cindex_comp_other, cindex)

  plot_list <- list()
  for (t in 1:length(dataset)) {
    cindex_comp_other_select <- cindex_comp_other[cindex_comp_other$ID == dataset[t], ]
    cindex_comp_other_select <- cindex_comp_other_select[order(cindex_comp_other_select$Cindex, decreasing = F), ]
    cindex_comp_other_select$Signature <- factor(cindex_comp_other_select$Signature, levels = unique(cindex_comp_other_select$Signature))

    plot_list[[dataset[t]]] <-
      ggplot(cindex_comp_other_select, aes(x = Signature, y = Cindex)) +
      geom_segment(aes(x = Signature, xend = Signature, y = 0, yend = Cindex, color = ID)) +
      geom_point(aes(color = ID)) +
      scale_color_manual(values = dataset_col[t], name = "Cohort") +
      theme(
        panel.grid = element_blank(),
        axis.text.y.left = (element_text(color = cindex_comp_other_select$text)),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        legend.position = "",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white")
      ) +
      # scale_y_continuous(position = "right")+
      labs(y = "C-index", x = "", title = "") +
      coord_flip()
  }

  p1 <- aplot::plot_list(gglist = plot_list, ncol = length(dataset))
  print(p1)
}
