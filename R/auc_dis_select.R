#' Area Under Curve distribution plot of specific selected method when predicting prognosis
#'
#' Creates a distribution plot of AUC among different datasets for specific selected method when predicting prognosis
#'
#' @param object A list of output of function cal_AUC_ml_res whose mode can be "all", "single" or "double"
#' @param model_name Model name to plot
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for cohorts
#' @param dataset A vector of names for all datasets
#' @param order If NULL, order is set to the default order. Otherwise, you can specify the order of datasets to plot 
#' @param year Specific year of auc in object list, such as c(1,3,...,9)
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
#'             model_name="StepCox[both] + plsRcox",
#'             dataset = names(list_train_vali_Data),
#'             order= names(list_train_vali_Data),
#'             year=c(1,3,5))
#'
auc_dis_select <- function(object, # a list of output of cal_AUC_ml_res  mode = "all",'single' or 'double'
                           model_name, # input specific model name
                           dataset_col = NULL, # color value for cohort
                           dataset, # input datasets name
                           order = NULL, # cohort order plot and input identical name of cohort in output of cal_AUC_ml_res
                           year # consistent number with object length, such as c(1,3,...,9)
) {
  library(ggplot2)

  if (is.null(dataset_col) == T) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    dataset_col <- dataset_col
  }

  auc_select <- data.frame()
  for (j in 1:length(object)) {
    object2 <- object[[j]]
    auc_all <- data.frame()
    for (i in names(object2)) {
      tmp <- object2[[i]]
      tmp2 <- data.frame()
      for (n in dataset) {
        tmp2 <- rbind(tmp2, data.frame(
          "AUC" = unique(tmp[[n]][["AUC"]]),
          "HR" = unique(tmp[[n]][["HR"]])
        ))
      }
      tmp2$ID <- dataset
      tmp2$Model <- i
      auc_all <- rbind(auc_all, tmp2)
    }
    auc_all$HR <- ifelse(auc_all$HR > 1, "HR>1", "HR<1")
    auc_all$AUC <- as.numeric(sprintf("%.2f", auc_all$AUC))

    auc_all_select <- auc_all[auc_all$Model == model_name, ]
    auc_all_select$Year <- paste(year[j], "Year", sep = "-")
    auc_select <- rbind(auc_select, auc_all_select)
  }


  if (is.null(order) == T) {
    order <- dataset
    auc_select$ID <- factor(auc_select$ID, levels = order) ## default order
  } else {
    order <- order
    auc_select$ID <- factor(auc_select$ID, levels = order)
  }


  p1 <- ggplot(data = auc_select, aes(x = ID, y = AUC, fill = ID)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = dataset_col, name = "Cohort") +
    geom_text(aes(label = AUC), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 3) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      legend.position = "",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    ) +
    # scale_y_continuous(position = "right")+
    labs(y = "AUC", x = "", title = model_name) +
    coord_flip() +
    facet_grid(. ~ Year)
  print(p1)
}
