#' Area Under Curve distribution plot when predicting prognosis
#'
#' Creates a distribution plot of auc among different methods when predicting prognosis
#'
#' @param object Output of function cal_AUC_ml_res whose mode is 'all'
#' @param color If NULL, color values are set to the default colors. Otherwise, you can specify three color values for cindex, two color values for mean cindex and two color values for HR 
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for datasets
#' @param dataset A vector of names for all datasets   
#' @param validate_set A vector of names for validate datasets  
#' @param order If NULL, order is set to the default order. Otherwise, you can specify the order of datasets to plot 
#' @param year Specific year of auc, it can be 1, 3, and 5
#' @param width Width of right plot. If NULL, width is set to the default width
#' @param height Height of right plot. If NULL, height is set to the default height
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' auc_dis_all(all.auc.1y,
#'          dataset = names(list_train_vali_Data),
#'          validate_set=names(list_train_vali_Data)[-1],
#'          order= names(list_train_vali_Data),
#'          width = 0.2,
#'          year=1)
#'
auc_dis_all <- function(object, # output of cal_AUC_ml_res mode = 'all'
                        color = NULL, # three color value for cindex, two color value for mean cindex and two color for HR
                        dataset_col = NULL, # color value for cohort
                        dataset, # input datasets name
                        validate_set, # input validate datasets name
                        order = NULL, # cohort order plot and input identical name of cohort in output of cal_AUC_ml_res
                        year, # year auc
                        width = NULL, # width of right plot
                        height = NULL # height of right plot
) {
  library(ggplot2)
  library(aplot)

  if (is.null(width) == T) {
    width <- 0.35 ## default 0.35
  } else {
    width <- width
  }

  if (is.null(height) == T) {
    height <- 0.01 ## default 0.01
  } else {
    height <- height
  }

  if (is.null(color) == T) {
    color <- c("#0084A7", "#F5FACD", "#E05D00", "#79AF97", "#8491B4", "#1B1919", "#808180") ## default color value
  } else {
    color <- color
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

  auc_all <- data.frame()
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
    tmp2$Model <- i
    auc_all <- rbind(auc_all, tmp2)
  }
  auc_all$HR <- ifelse(auc_all$HR > 1, "HR>1", "HR<1")
  auc_all$AUC <- as.numeric(sprintf("%.2f", auc_all$AUC))

  mean_AUC <- aggregate(
    x = auc_all$AUC,
    by = list(auc_all$Model),
    FUN = mean
  )
  colnames(mean_AUC) <- c("Model", "mean")
  mean_AUC$mean <- as.numeric(sprintf("%.3f", mean_AUC$mean))
  mean_AUC$Value <- "Mean AUC in all cohorts"

  AUC_vaidate <- auc_all[auc_all$ID %in% validate_set, ]
  mean_AUC_vaidate <- aggregate(
    x = AUC_vaidate$AUC,
    by = list(AUC_vaidate$Model),
    FUN = mean
  )
  colnames(mean_AUC_vaidate) <- c("Model", "mean")
  mean_AUC_vaidate$mean <- as.numeric(sprintf("%.3f", mean_AUC_vaidate$mean))
  mean_AUC_vaidate <- mean_AUC_vaidate[order(mean_AUC_vaidate$mean, decreasing = F), ]
  mean_AUC_vaidate$Value <- "Mean AUC in validate cohorts"

  labels <- as.data.frame(unique(auc_all$ID))
  colnames(labels) <- "Cohort"

  if (is.null(order) == T) {
    labels <- labels ## default order
  } else {
    order <- order
    labels$Cohort <- factor(labels$Cohort, levels = order)
    auc_all$ID <- factor(auc_all$ID, levels = order)
  }

  auc_all$Model <- factor(auc_all$Model, levels = mean_AUC_vaidate$Model)
  mean_AUC$Model <- factor(mean_AUC$Model, levels = mean_AUC_vaidate$Model)
  mean_AUC_vaidate$Model <- factor(mean_AUC_vaidate$Model, levels = mean_AUC_vaidate$Model)

  p1 <- ggplot(auc_all, aes(x = ID, y = Model)) +
    geom_tile(aes(fill = AUC), color = "white", size = 0.5) +
    geom_text(aes(label = AUC, color = HR), vjust = 0.5, size = 2) +
    scale_color_manual(values = color[6:7]) +
    scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = median(auc_all$AUC), name = "AUC") +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_blank()
    )

  p2 <- ggplot(labels, aes(Cohort, y = 1)) +
    geom_tile(aes(fill = Cohort), color = "white", size = 0.5) +
    scale_fill_manual(values = dataset_col, name = "Cohort") +
    labs(title = paste(year, "Year AUC predicted by all models")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  p3 <- ggplot(data = mean_AUC, aes(x = Model, y = mean, fill = Value)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = color[4]) +
    geom_text(aes(label = mean), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 2) +
    theme(
      axis.title = element_text(size = 8),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      # panel.border = element_rect(colour = "black", fill=NA,size = 0.2),
      panel.grid = element_blank(),
      # legend.position = "",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    ) +
    # scale_y_continuous(#position = "right",
    #   expand = c(0, 0.001))+
    labs(y = "") +
    coord_flip()

  p4 <- ggplot(data = mean_AUC_vaidate, aes(x = Model, y = mean, fill = Value)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = color[5], name = "") +
    geom_text(aes(label = mean), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 2) +
    theme(
      axis.title = element_text(size = 8),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      # panel.border = element_rect(colour = "black", fill=NA,size = 0.2),
      panel.grid = element_blank(),
      # legend.position = "",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    ) +
    # scale_y_continuous(#position = "right",
    #   expand = c(0, 0.001))+
    labs(y = "") +
    coord_flip()

  print(p1 %>%
          insert_top(p2, height = height) %>%
          insert_right(p3, width = width) %>%
          insert_right(p4, width = width))
}
