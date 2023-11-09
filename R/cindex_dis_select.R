#' Title
#'
#' @param object
#' @param dataset_col
#' @param model
#' @param order
#'
#' @return
#' @export
#'
#' @examples
cindex_dis_select <- function(object, # output of ML.Dev.Prog.Sig mode = "all",'single' or 'double'
                              dataset_col = NULL, # color value for cohort
                              model, ## if input object is from mode=all, then define model as specific model name while SOD
                              order = NULL # cohort order plot and input identical name of cohort in output of ML.Dev.Prog.Sig
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


  if (model != "SOD") {
    cindex <- object[["Cindex.res"]]
    cindex <- cindex[cindex$Model == model, ]
    cindex$Cindex <- as.numeric(sprintf("%.2f", cindex$Cindex))

    if (is.null(order) == T) {
      cindex <- cindex ## default order
    } else {
      order <- order
      cindex$ID <- factor(cindex$ID, levels = order)
    }
    p1 <- ggplot(data = cindex, aes(x = ID, y = Cindex, fill = ID)) +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = dataset_col, name = "Cohort") +
      geom_text(aes(label = Cindex), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 3) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        legend.position = "",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white")
      ) +
      # scale_y_continuous(position = "right")+
      labs(y = "C-index", x = "", title = model) +
      coord_flip()
    print(p1)
  } else {
    cindex <- object[["Cindex.res"]]
    cindex$Cindex <- as.numeric(sprintf("%.2f", cindex$Cindex))

    if (is.null(order) == T) {
      cindex <- cindex ## default order
    } else {
      order <- order
      cindex$ID <- factor(cindex$ID, levels = order)
    }

    p1 <- ggplot(data = cindex, aes(x = ID, y = Cindex, fill = ID)) +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = dataset_col, name = "Cohort") +
      geom_text(aes(label = Cindex), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 3) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        legend.position = "",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white")
      ) +
      # scale_y_continuous(position = "right")+
      labs(y = "C-index", x = "", title = cindex$Model) +
      coord_flip()
    print(p1)
  }
}
