#' Compare the AUC of our model with previously published models
#'
#' Creates a distribution plot of AUC among different mdoels 
#'
#' @param object Output of function cal_auc_pre.prog.sig
#' @param object2 Output of function cal_AUC_ml_res
#' @param model_name Model name in object2 used to compare
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for datasets
#' @param dataset A vector of names for all datasets 
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' auc_comp(auc.glioma.lgg.gbm.1,
#'          all.auc.1y,
#'          model_name="StepCox[both] + plsRcox",
#'          dataset=names(list_train_vali_Data))
#' 
auc_comp <- function(object, # output of cal_auc_pre.prog.sig
                     object2, # output of cal_AUC_ml_res
                     model_name, ## input specific model name
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


  auc_comp_other <- data.frame()
  for (i in names(object)) {
    tmp <- object[[i]]
    tmp2 <- data.frame()
    for (n in dataset) {
      tmp2 <- rbind(tmp2, data.frame(
        "AUC" = unique(tmp[[n]][["AUC"]]),
        "HR" = unique(tmp[[n]][["HR"]])
      ))
    }
    tmp2$ID <- dataset
    tmp2$Signature <- i
    auc_comp_other <- rbind(auc_comp_other, tmp2)
  }

  auc_comp_other$text <- "black"

  object3 <- object2[[model_name]]
  auc_all <- data.frame()
  for (n in dataset) {
    auc_all <- rbind(auc_all, data.frame(
      "AUC" = unique(object3[[n]][["AUC"]]),
      "HR" = unique(object3[[n]][["HR"]]),
      "ID" = n
    ))
  }
  auc_all$Signature <- "Our Model"
  auc_all$text <- "red"

  auc_comp_other <- rbind(auc_comp_other, auc_all)

  plot_list <- list()
  for (t in 1:length(dataset)) {
    auc_comp_other_select <- auc_comp_other[auc_comp_other$ID == dataset[t], ]
    auc_comp_other_select <- auc_comp_other_select[order(auc_comp_other_select$AUC, decreasing = F), ]
    auc_comp_other_select$Signature <- factor(auc_comp_other_select$Signature, levels = unique(auc_comp_other_select$Signature))

    plot_list[[dataset[t]]] <-
      ggplot(auc_comp_other_select, aes(x = Signature, y = AUC)) +
      geom_segment(aes(x = Signature, xend = Signature, y = 0, yend = AUC, color = ID)) +
      geom_point(aes(color = ID)) +
      scale_color_manual(values = dataset_col[t], name = "Cohort") +
      theme(
        panel.grid = element_blank(),
        axis.text.y.left = (element_text(color = auc_comp_other_select$text)),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        legend.position = "",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white")
      ) +
      # scale_y_continuous(position = "right")+
      labs(y = "AUC", x = "", title = "") +
      coord_flip()
  }

  p1 <- aplot::plot_list(gglist = plot_list, ncol = length(dataset))
  print(p1)
}
