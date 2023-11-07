#' Title
#'
#' @param object
#' @param color
#' @param dataset_col
#' @param validate_set
#' @param order
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
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

rs_sur <- function(object, # output of ML.Dev.Prog.Sig  mode = "all",'single' or 'double'
                   model_name, # input identical model name in output of ML.Dev.Prog.Sig
                   dataset, # input identical name of cohort in output of ML.Dev.Prog.Sig
                   cutoff = NULL, # cutoff for group: mean or quantile value
                   conf.int = F, # show confidence intervals T or F
                   median.line = NULL, # drawing a horizontal/vertical line at median survival. c("none", "hv", "h", "v")
                   color = NULL, # two color value for high and low group
                   xlab = NULL, # x axis title
                   pval.coord = NULL # p value position
) {
  library(survival)
  library(survminer)

  if (is.null(color) == T) {
    color <- c("#868686", "#B24745") ## default color value
  } else {
    color <- color
  }

  if (is.null(median.line) == T) {
    median.line <- "none" ## default not show
  } else {
    median.line <- median.line
  }

  if (is.null(xlab) == T) {
    xlab <- "Time"
  } else {
    xlab <- xlab
  }

  if (is.null(pval.coord) == T) {
    pval.coord <- c(1, 0.25)
  } else {
    pval.coord <- pval.coord
  }

  tmp <- object[["riskscore"]][[model_name]][[dataset]]
  mySurv <- Surv(tmp$OS.time, tmp$OS)

  if (is.null(cutoff) == T) {
    cutoff <- 0.5
    value <- quantile(tmp$RS, probs = c(cutoff))
  } else {
    if (cutoff == "mean") {
      value <- mean(tmp$RS)
    } else {
      cutoff <- cutoff
      value <- quantile(tmp$RS, probs = c(cutoff))
    }
  }

  tmp$Group <- ifelse(tmp$RS > value, "High", "Low")
  Group <- tmp$Group
  Group <- factor(Group, levels = c("Low", "High"))
  tmp$Group <- factor(Group, levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = tmp)

  # calculate HR and 95%CI
  data.survdiff <- survdiff(mySurv ~ Group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")

  p1 <- ggsurvplot(fit,
    data = tmp,
    conf.int = conf.int,
    censor = F,
    palette = color,
    legend.title = paste0("Riskscore in ", dataset),
    font.legend = 11,
    surv.median.line = median.line,
    legend.labs = c(
      paste0("Low", "(n=", fit$n[1], ")"),
      paste0("High", "(n=", fit$n[2], ")")
    ),
    pval = paste(
      pval = ifelse(p.val < 0.001, "p < 0.001",
        paste("p = ", round(p.val, 3), sep = "")
      ),
      HR, CI, sep = "\n"
    ),
    xlab = xlab,
    title = model_name,
    pval.coord = pval.coord
  )
  p1 <- p1$plot + theme(
    plot.title = element_text(hjust = 0.5)
  )
  print(p1)
}

roc_vis <- function(object, # output of cal_AUC_ml_res  mode = "all",'single' or 'double'
                    model_name, # input one identical model name in output of cal_AUC_ml_res
                    dataset, # input identical name of cohort in output of cal_AUC_ml_res
                    order = NULL, # cohort order plot and input identical name of cohort in output of cal_AUC_ml_res
                    dataset_col = NULL, # color value for dataset
                    anno_position = NULL, # auc value position
                    year # year auc
) {
  library(ggplot2)

  if (is.null(anno_position) == T) {
    anno_position <- c(0.53, 0.35)
  } else {
    anno_position <- anno_position
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

  tmp <- data.frame()
  for (i in 1:length(dataset)) {
    tmp2 <- object[[model_name]][[dataset[i]]]
    tmp2$dataset <- dataset[i]
    tmp <- rbind(tmp, tmp2)
  }

  if (is.null(order) == T) {
    order <- dataset
    tmp$dataset <- factor(tmp$dataset, levels = order) ## default order
  } else {
    order <- order
    tmp$dataset <- factor(tmp$dataset, levels = order)
  }

  p1 <- ggplot(tmp, aes(FP, TP, color = dataset)) +
    geom_line(size = 0.5, alpha = 1) +
    scale_color_manual(values = dataset_col, name = "Cohort") +
    # geom_area(aes(fill = dataset),alpha = 0.5,stat = "align",position = "identity")+
    labs(
      title = paste(year, "Year AUC predicted by", model_name),
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      plot.title = element_text(hjust = 0.5),
      legend.position = "",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    scale_x_continuous(expand = c(0.01, 0)) + # right most position will be 1 + (1-0) * 0.01 + 0 = 10
    scale_y_continuous(expand = c(0.01, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 0.5)

  for (i in 1:length(dataset)) {
    ano_y <- anno_position[2] - (i - 1) * 0.05
    p1 <- p1 +
      annotate("segment",
        x = anno_position[1], xend = anno_position[1] + 0.05, y = ano_y, yend = ano_y,
        colour = dataset_col[i]
      ) +
      annotate("text",
        x = anno_position[1] + 0.07, y = ano_y, hjust = 0,
        label = paste("AUC in", order[i], ":", round(unique(tmp[tmp$dataset == order[i], "AUC"]), 3))
      )
  }
  print(p1)
}

roc_vis_category <- function(object, # output of ML.Dev.Pred.Category.Sig
                             model_name, # input one identical model name in output of ML.Dev.Pred.Category.Sig
                             dataset, # input identical name of cohort in output of ML.Dev.Pred.Category.Sig
                             order = NULL, # cohort order plot and input identical name of cohort in output of ML.Dev.Pred.Category.Sig
                             dataset_col = NULL, # color value for  dataset
                             anno_position = NULL # auc value position
) {
  library(ggplot2)

  if (is.null(anno_position) == T) {
    anno_position <- c(0.53, 0.35)
  } else {
    anno_position <- anno_position
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

  if (model_name != "cancerclass") {
    tmp <- data.frame()
    for (i in 1:length(dataset)) {
      roc <- cbind(
        object[["roc"]][[dataset[i]]][[model_name]][["TPR"]], object[["roc"]][[dataset[i]]][[model_name]][["FPR"]],
        object[["roc"]][[dataset[i]]][[model_name]][["AUC"]]
      )
      roc <- as.data.frame(roc)
      colnames(roc) <- c("TP", "FP", "AUC")
      roc$dataset <- dataset[i]
      tmp <- rbind(tmp, roc)
    }
  } else {
    tmp <- data.frame()
    for (i in 1:length(dataset)) {
      roc <- cbind(
        object[["roc"]][[dataset[i]]][[model_name]][["sensitivities"]],
        1 - object[["roc"]][[dataset[i]]][[model_name]][["specificities"]],
        object[["roc"]][[dataset[i]]][[model_name]][["auc"]]
      )
      roc <- as.data.frame(roc)
      colnames(roc) <- c("TP", "FP", "AUC")
      roc$dataset <- dataset[i]
      tmp <- rbind(tmp, roc)
    }
  }

  if (is.null(order) == T) {
    order <- dataset
    tmp$dataset <- factor(tmp$dataset, levels = order) ## default order
  } else {
    order <- order
    tmp$dataset <- factor(tmp$dataset, levels = order)
  }

  p1 <- ggplot(tmp, aes(FP, TP, color = dataset)) +
    geom_line(size = 0.5, alpha = 1) +
    scale_color_manual(values = dataset_col, name = "Cohort") +
    # geom_area(aes(fill = dataset),alpha = 0.5,stat = "align",position = "identity")+
    labs(
      title = paste("AUC predicted by", model_name),
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      plot.title = element_text(hjust = 0.5),
      legend.position = "",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    scale_x_continuous(expand = c(0.01, 0)) + # right most position will be 1 + (1-0) * 0.01 + 0 = 10
    scale_y_continuous(expand = c(0.01, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 0.5)

  for (i in 1:length(dataset)) {
    ano_y <- anno_position[2] - (i - 1) * 0.05
    p1 <- p1 +
      annotate("segment",
        x = anno_position[1], xend = anno_position[1] + 0.05, y = ano_y, yend = ano_y,
        colour = dataset_col[i]
      ) +
      annotate("text",
        x = anno_position[1] + 0.07, y = ano_y, hjust = 0,
        label = paste("AUC in", order[i], ":", round(unique(tmp[tmp$dataset == order[i], "AUC"]), 3))
      )
  }
  print(p1)
}

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


HR_com <- function(object, # output of cal_auc_pre.prog.sig
                   object2, # output of cal_AUC_ml_res
                   color = NULL, # three color value for HR
                   model_name, ## if input object2 is from mode=all, then define model as specific model name while SOD
                   dataset_col = NULL, # color value for cohort
                   dataset, # input datasets name
                   type # c('categorical'),c('continuous) Univariate cox regression analysis, suggest categorical group by median
) {
  library(survival)
  library(survminer)
  library(ggplot2)

  if (is.null(color) == T) {
    color <- c("#0084A7", "white", "#E05D00") ## default color value
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

  hr <- data.frame()
  for (i in names(object)) {
    for (n in dataset) {
      rs <- object[[i]][[n]]

      if (type == "categorical") {
        # OS
        group <- ifelse(rs$RS > median(rs$RS), "high", "low")

        if (length(unique(group)) == 1) {
          p.val <- 1
          HR <- 1

          hr <- rbind(hr, data.frame(
            "P" = p.val,
            "HR" = HR,
            "ID" = n,
            "Signature" = i
          ))
        } else {
          group <- factor(group, levels = c("low", "high"))

          my.surv <- Surv(rs$OS.time, rs$OS)
          univ_cox <- coxph(my.surv ~ group, data = rs)
          univ_coxSummary <- summary(univ_cox)
          p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
          HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

          hr <- rbind(hr, data.frame(
            "P" = p.val,
            "HR" = HR,
            "ID" = n,
            "Signature" = i
          ))
        }
      } else {
        my.surv <- Surv(rs$OS.time, rs$OS)
        univ_cox <- coxph(my.surv ~ RS, data = rs)
        univ_coxSummary <- summary(univ_cox)
        p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
        HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

        hr <- rbind(hr, data.frame(
          "P" = p.val,
          "HR" = HR,
          "ID" = n,
          "Signature" = i
        ))
      }
    }
  }
  hr$text <- "black"

  hr2 <- data.frame()
  if (model_name != "SOD") {
    object3 <- object2[["riskscore"]][[model_name]]
  } else {
    object3 <- object2[["riskscore"]][[1]]
  }
  for (n in dataset) {
    rs <- object3[[n]]

    if (type == "categorical") {
      # OS
      group <- ifelse(rs$RS > median(rs$RS), "high", "low")

      if (length(unique(group)) == 1) {
        p.val <- 1
        HR <- 1

        hr2 <- rbind(hr2, data.frame(
          "P" = p.val,
          "HR" = HR,
          "ID" = n,
          "Signature" = "Our model"
        ))
      } else {
        group <- factor(group, levels = c("low", "high"))

        my.surv <- Surv(rs$OS.time, rs$OS)
        univ_cox <- coxph(my.surv ~ group, data = rs)
        univ_coxSummary <- summary(univ_cox)
        p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
        HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

        hr2 <- rbind(hr2, data.frame(
          "P" = p.val,
          "HR" = HR,
          "ID" = n,
          "Signature" = "Our model"
        ))
      }
    } else {
      my.surv <- Surv(rs$OS.time, rs$OS)
      univ_cox <- coxph(my.surv ~ RS, data = rs)
      univ_coxSummary <- summary(univ_cox)
      p.val <- as.numeric(univ_coxSummary$coefficients[, "Pr(>|z|)"])[1]
      HR <- as.numeric(univ_coxSummary$coefficients[, "exp(coef)"])[1]

      hr2 <- rbind(hr2, data.frame(
        "P" = p.val,
        "HR" = HR,
        "ID" = n,
        "Signature" = "Our model"
      ))
    }
  }
  hr2$text <- "red"

  hr <- rbind(hr, hr2)

  hr$pstar <- ifelse(hr$P < 0.05,
    ifelse(hr$P < 0.01, "**", "*"),
    ""
  )
  tmp <- hr
  tmp$Signature <- factor(hr$Signature, levels = c("Our model", unique(names(object))))
  tmp$ID <- factor(tmp$ID, levels = dataset)

  p <- ggplot(tmp, aes(Signature, ID)) +
    geom_tile(aes(fill = HR), colour = "white", size = 1) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = 1) +
    geom_text(aes(label = pstar), col = "black", size = 5) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(colour = dataset_col),
      axis.text.x = element_text(colour = c("red", rep("black", length(names(object)))), angle = 45, hjust = 1)
    ) +
    labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "HR"))
  print(p)
}

meta_unicox_vis <- function(object, # output of cal_unicox_meta_ml_res
                            dataset_col = NULL, # color value for cohorts
                            dataset # input datasets name
) {
  library(forestploter)
  library(grid)

  if (is.null(dataset_col) == T) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    dataset_col <- dataset_col
  }

  metamodel <- object
  dt <- metamodel[["data"]]
  dt$Weight_random <- paste(round(100 * metamodel$w.random / sum(metamodel$w.random), 2), "%", sep = "")
  dt$Weight_fixed <- paste(round(100 * metamodel$w.fixed / sum(metamodel$w.fixed), 2), "%", sep = "")

  dt[nrow(dt) + 1, ] <- NA
  dt[nrow(dt), 6] <- "Random effect model"
  dt[nrow(dt), 1] <- exp(metamodel$TE.random)
  dt[nrow(dt), 3] <- exp(metamodel$lower.random)
  dt[nrow(dt), 4] <- exp(metamodel$upper.random)
  dt[nrow(dt), 2] <- metamodel$pval.random
  dt[nrow(dt), 10] <- "100%"
  dt[nrow(dt), 11] <- "--"

  dt[nrow(dt) + 1, ] <- NA
  dt[nrow(dt), 6] <- "Fixed effect model"
  dt[nrow(dt), 1] <- exp(metamodel$TE.fixed)
  dt[nrow(dt), 3] <- exp(metamodel$lower.fixed)
  dt[nrow(dt), 4] <- exp(metamodel$upper.fixed)
  dt[nrow(dt), 2] <- metamodel$pval.fixed
  dt[nrow(dt), 10] <- "--"
  dt[nrow(dt), 11] <- "100%"

  dt$se <- (log(dt$HCI) - log(dt$HR)) / 1.96
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  dt$`HR (95% CI)` <- ifelse(is.na(dt$Group), "",
    sprintf(
      "%.2f (%.2f - %.2f)",
      dt$HR, dt$LCI, dt$HCI
    )
  )
  dt$p <- ifelse(dt$pvalue < 0.001, "P<0.001", sprintf("%.3f", dt$pvalue))
  colnames(dt)[c(6:8, 10:11, 15)] <- c("Cohorts", "TE", "SE(TE)", "Weight(random)", "Weight(fixed)", "P value")
  dt$TE <- round(dt$TE, 2)
  dt$`SE(TE)` <- round(dt$`SE(TE)`, 2)
  dt[c((nrow(dt) - 1):nrow(dt)), c(7, 8)] <- ""
  rownames(dt) <- dt$Cohorts
  dt <- dt[c(dataset, "Random effect model", "Fixed effect model"), ]

  tm <- forest_theme(
    core = list(bg_params = list(
      fill = c(dataset_col[1:length(dataset)], "grey", "grey"),
      alpha = 0.5
    )),
    base_size = 10,
    # Confidence interval point shape, line type/color/width
    ci_pch = 16,
    ci_col = "#762a83",
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = 0.2, # Set an T end at the end of CI
    # Reference line width/type/color
    refline_lwd = 1,
    refline_lty = "dashed",
    refline_col = "grey20",
    # Vertical line width/type/color
    vertline_lwd = 1,
    vertline_lty = "dashed",
    vertline_col = "grey20",
    # Change summary color for filling and borders
    summary_fill = "#4575b4",
    summary_col = "#4575b4",
    # Footnote font size/face/color
    footnote_cex = 1,
    footnote_fontface = "italic",
    footnote_col = "red"
  )

  p <- forestploter::forest(dt[, c(6:8, 10:11, 13, 14, 15)],
    est = dt$HR,
    lower = dt$LCI,
    upper = dt$HCI,
    sizes = dt$se,
    is_summary = c(rep(FALSE, nrow(dt) - 2), TRUE, TRUE),
    ci_column = 6,
    ref_line = 1,
    arrow_lab = c("Better", "Worse"),
    x_trans = "log2",
    xlim = c(0, ceiling(max(dt$HCI))),
    ticks_at = c(0.5, 2^seq(0, floor(log2((ceiling(max(dt$HCI)) - 1))), by = 1)),
    footnote = " Univariate Cox Regression",
    theme = tm
  )

  p <- add_text(p,
    text = "Meta analysis of univariate Cox regression",
    part = "header",
    row = 0,
    col = 4:6,
    just = c("center"),
    gp = gpar(fontface = "bold")
  )

  p <- add_border(p,
    part = "header",
    row = c(0, 1),
    gp = gpar(lwd = 1)
  )

  p <- insert_text(p,
    text = "Meta analysis",
    row = length(dataset) + 1,
    just = "left",
    gp = gpar(cex = 0.8, col = "blue", fontface = "italic")
  )
  print(p)
}

core_feature_rank <- function(object, # output of ML.Corefeature.Prog.Screen
                              col = NULL, # color value for segment and point
                              top = NULL # top number gene to show
) {
  library(ggplot2)

  if (is.null(col) == T) {
    col <- c("#B09C85", "#E18727")
  } else {
    col <- col
  }

  if (is.null(top) == T) {
    top <- 50
  } else {
    top <- top
  }

  tmp <- as.data.frame(table(object$selected.fea))
  tmp$Var1 <- gsub("\\.", "-", tmp$Var1)
  tmp <- tmp[order(tmp$Freq, decreasing = F), ]
  tmp$Var1 <- factor(tmp$Var1, levels = tmp$Var1)
  colnames(tmp)[2] <- "Frequence"

  p1 <- ggplot(tmp[(nrow(tmp) - top + 1):nrow(tmp), ], aes(x = Var1, y = Frequence)) +
    geom_segment(aes(x = Var1, xend = Var1, y = 0, yend = Frequence), color = col[1]) +
    geom_point(aes(size = Frequence), color = col[2], alpha = 0.5) +
    DOSE::theme_dose(10) +
    # scale_color_manual(values=dataset_col[t],name="Cohort")+
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      # legend.position = "",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    ) +
    # scale_y_continuous(position = "right")+
    labs(y = "Frequence of screening", x = "", title = paste("Top", top, "selected genes")) +
    coord_flip()

  print(p1)
}

core_feature_sur <- function(gene, # gene name
                             InputMatrix, # InputMatrix with os os.time and gene expression need contain gene
                             dataset, # the cohort name of InputMatrix
                             cutoff = NULL, # cutoff for group: mean or quantile value
                             conf.int = F, # show confidence intervals T or F
                             median.line = NULL, # drawing a horizontal/vertical line at median survival. c("none", "hv", "h", "v")
                             color = NULL, # two color value for high and low group
                             xlab = NULL, # x axis title
                             pval.coord = NULL # p value position
) {
  library(survival)
  library(survminer)

  if (is.null(color) == T) {
    color <- c("#868686", "#B24745") ## default color value
  } else {
    color <- color
  }

  if (is.null(median.line) == T) {
    median.line <- "none" ## default not show
  } else {
    median.line <- median.line
  }

  if (is.null(xlab) == T) {
    xlab <- "Time"
  } else {
    xlab <- xlab
  }

  if (is.null(pval.coord) == T) {
    pval.coord <- c(1, 0.25)
  } else {
    pval.coord <- pval.coord
  }


  tmp <- InputMatrix[, c("OS.time", "OS", gene)]
  colnames(tmp)[3] <- "gene"
  tmp <- tmp[!is.na(tmp$gene), ]
  mySurv <- Surv(tmp$OS.time, tmp$OS)

  if (is.null(cutoff) == T) {
    cutoff <- 0.5
    value <- quantile(tmp$gene, probs = c(cutoff))
  } else {
    if (cutoff == "mean") {
      value <- mean(tmp$gene)
    } else {
      cutoff <- cutoff
      value <- quantile(tmp$gene, probs = c(cutoff))
    }
  }

  tmp$Group <- ifelse(tmp$gene > value, "High", "Low")
  Group <- tmp$Group
  Group <- factor(Group, levels = c("Low", "High"))
  tmp$Group <- factor(Group, levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = tmp)

  # calculate HR and 95%CI
  data.survdiff <- survdiff(mySurv ~ Group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")

  p1 <- ggsurvplot(fit,
    data = tmp,
    conf.int = conf.int,
    censor = F,
    palette = color,
    legend.title = paste0("Expression in ", dataset),
    font.legend = 11,
    surv.median.line = median.line,
    legend.labs = c(
      paste0("Low", "(n=", fit$n[1], ")"),
      paste0("High", "(n=", fit$n[2], ")")
    ),
    pval = paste(
      pval = ifelse(p.val < 0.001, "p < 0.001",
        paste("p = ", round(p.val, 3), sep = "")
      ),
      HR, CI, sep = "\n"
    ),
    xlab = xlab,
    title = gene,
    pval.coord = pval.coord
  )
  p1 <- p1$plot + theme(
    plot.title = element_text(hjust = 0.5)
  )

  print(p1)
}

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
