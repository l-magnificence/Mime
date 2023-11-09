#' C-index distribution plot
#'
#' Creates a distribution plot of C-index among different methods
#'
#' @param object Output of function ML.Dev.Prog.Sig whose mode is 'all'
#' @param color If NULL, color values are set to the default colors. Otherwise, you can specify three color values for cindex and two color values for mean cindex 
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for cohorts
#' @param validate_set A vector of names for validate datasets   
#' @param order If NULL, order is set to the default order. Otherwise, you can specify the order of cohorts to plot
#' @param width Width of right plot. If NULL, width is set to the default width
#' @param height Height of right plot. If NULL, height is set to the default height
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.2)
#'
cindex_dis_all <- function(object, # output of ML.Dev.Prog.Sig mode = 'all'
                           color = NULL, # three color value for cindex and two color value for mean cindex
                           dataset_col = NULL, # color value for cohort
                           validate_set, # input validate datasets name
                           order = NULL, # cohort order plot and input identical name of cohort in output of ML.Dev.Prog.Sig
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
    color <- c("#0084A7", "#F5FACD", "#E05D00", "#79AF97", "#8491B4") ## default color value
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

  cindex <- object[["Cindex.res"]]
  cindex$Cindex <- as.numeric(sprintf("%.2f", cindex$Cindex))

  mean_cindex <- aggregate(
    x = cindex$Cindex,
    by = list(cindex$Model),
    FUN = mean
  )
  colnames(mean_cindex) <- c("Model", "mean")
  mean_cindex$mean <- as.numeric(sprintf("%.3f", mean_cindex$mean))
  mean_cindex$Value <- "Mean C-index in all cohorts"

  cindex_vaidate <- cindex[cindex$ID %in% validate_set, ]
  mean_cindex_validate <- aggregate(
    x = cindex_vaidate$Cindex,
    by = list(cindex_vaidate$Model),
    FUN = mean
  )
  colnames(mean_cindex_validate) <- c("Model", "mean")
  mean_cindex_validate$mean <- as.numeric(sprintf("%.3f", mean_cindex_validate$mean))
  mean_cindex_validate <- mean_cindex_validate[order(mean_cindex_validate$mean, decreasing = F), ]
  mean_cindex_validate$Value <- "Mean C-index in validate cohorts"

  labels <- as.data.frame(unique(cindex$ID))
  colnames(labels) <- "Cohort"

  if (is.null(order) == T) {
    labels <- labels ## default order
  } else {
    order <- order
    labels$Cohort <- factor(labels$Cohort, levels = order)
    cindex$ID <- factor(cindex$ID, levels = order)
  }

  cindex$Model <- factor(cindex$Model, levels = mean_cindex_validate$Model)
  mean_cindex$Model <- factor(mean_cindex$Model, levels = mean_cindex_validate$Model)
  mean_cindex_validate$Model <- factor(mean_cindex_validate$Model, levels = mean_cindex_validate$Model)

  p1 <- ggplot(cindex, aes(x = ID, y = Model)) +
    geom_tile(aes(fill = Cindex), color = "white", size = 0.5) +
    geom_text(aes(label = Cindex), vjust = 0.5, color = "black", size = 2) +
    scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = median(cindex$Cindex), name = "C-index") +
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
    theme_void()

  p3 <- ggplot(data = mean_cindex, aes(x = Model, y = mean, fill = Value)) +
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

  p4 <- ggplot(data = mean_cindex_validate, aes(x = Model, y = mean, fill = Value)) +
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
