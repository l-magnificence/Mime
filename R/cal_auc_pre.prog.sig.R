#' Calculate AUC scores of the previous signarues in all data
#'
#' @param use_your_own_collected_sig Whether to use your own collected signatures. T or F.
#' @param collected_sig_table If use_your_own_collected_sig set as T, you should provide a data frame containing the information of the signatures. The column names of the data frame are "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol". For example, 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'. The 'model' consists of the first name of the first author and the PMID. The 'PMID' is from the paper. The 'Cancer' uses abbreviations in the format of the TCGA. 'Author' is the first name of the first author of the paper. 'Coef' is the coefficient of the variable in the signature. 'symbol' is the variable in the signature. If use_your_own_collected_sig is F, you don't need to provide this data. 
#' @param type.sig If the use_your_own_collected_sig is F, here we collected some signatures of the LGG, GBM, and Glioma. You can choose c('Glioma', 'LGG', 'GBM'), c('Glioma'),c('GBM'), c('Glioma', 'LGG'), and some other combination for the signatures you want.
#' @param list_input_data A list of the cohorts. Column names are 'ID', 'OS.time', 'OS', and the other variables. OS.time (Day). OS(1/0).
#' @param AUC_time  c(1,2,3,4,5,6,7,······), for 1 year, 2 years, 3 years......We recommend using the shortest survival time among all queues.
#' @param auc_cal_method 'KM', 'NNE'. The default is 'KM'.
#'
#' @return A list of the AUC results of each previous signature in each cohort you provide.
#' @export
#'
#' @examples
cal_auc_pre.prog.sig <- function(use_your_own_collected_sig, # 是否使用您自己收集的signature， T or F
                                 collected_sig_table, # 列名分别为
                                 # "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol"
                                 # 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'
                                 type.sig, ### prognostic signature 的类型，c('Glioma','LGG','GBM')， c('Glioma'),c('Glioma','LGG')
                                 list_input_data, # list of the cohorts(ID,OS.time, OS····)经过了log2（x+1）转化
                                 AUC_time = NULL, ### 时间 年份
                                 auc_cal_method = "KM" # KM, NNE 默认为KM
) {
  library(tidyverse)
  library(survival)
  library(survivalROC)


  if (use_your_own_collected_sig) {
    sig.input <- collected_sig_table
  } else {
    pre.prog.sig <- Mime1::pre.prog.sig

    if (all(type.sig %in% names(pre.prog.sig))) {
      if (length(type.sig) == 1) {
        sig.input <- pre.prog.sig[[type.sig[1]]]
      } else {
        sig.input <- pre.prog.sig[[type.sig[1]]]

        for (i in 2:length(type.sig)) {
          sig.input <- rbind(sig.input, pre.prog.sig[[type.sig[i]]])
        }
      }
    } else {
      stop("please provide correct type.sig")
    }
  }





  sig.input$Coef <- as.numeric(sig.input$Coef)
  sig.input$symbol <- gsub("-", ".", sig.input$symbol)



  # Replace '-' in column names with '.'
  list_input_data <- lapply(list_input_data, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    return(x)
  })

  common_feature <- c("ID", "OS.time", "OS", unique(sig.input$symbol))

  # for (i in names(list_input_data)) {
  #   common_feature = intersect(common_feature, colnames(list_input_data[[i]]))
  # }

  returnIDtoRS <- function(rs.table.list, rawtableID) {
    for (i in names(rs.table.list)) {
      rs.table.list[[i]]$ID <- rawtableID[[i]]$ID
      rs.table.list[[i]] <- rs.table.list[[i]] %>% dplyr::select("ID", everything())
    }

    return(rs.table.list)
  }


  returnRStoROC <- function(rs.table.list, AUC_time) {
    roc.rs <- lapply(rs.table.list, function(x) {
      mySurv <- Surv(x$OS.time, x$OS)




      x$Group <- ifelse(x$RS > median(x$RS), "High", "Low")


      if (length(unique(x$Group)) > 1) {
        x$Group <- factor(x$Group, levels = c("Low", "High"))
      } else {
        x$Group <- ifelse(x$RS > mean(x$RS), "High", "Low")
      }

      if (length(unique(x$Group)) > 1) {
        x$Group <- factor(x$Group, levels = c("Low", "High"))
        data.survdiff <- survdiff(mySurv ~ x$Group)
        HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])

        if (HR > 1) {
          if (auc_cal_method == "NNE") {
            risk.survivalROC <- survivalROC(
              Stime = x$OS.time,
              status = x$OS,
              marker = x$RS,
              predict.time = 365 * AUC_time,
              method = "NNE", span = 0.25 * nrow(x)^(-0.20)
            )
          } else if (auc_cal_method == "KM") {
            risk.survivalROC <- survivalROC(
              Stime = x$OS.time,
              status = x$OS,
              marker = x$RS,
              predict.time = 365 * AUC_time,
              method = "KM"
            )
          } else {
            print("Please provide the correct parameters for method")
          }
        } else {
          if (auc_cal_method == "NNE") {
            risk.survivalROC <- survivalROC(
              Stime = x$OS.time,
              status = x$OS,
              marker = -x$RS,
              predict.time = 365 * AUC_time,
              method = "NNE", span = 0.25 * nrow(x)^(-0.20)
            )
          } else if (auc_cal_method == "KM") {
            risk.survivalROC <- survivalROC(
              Stime = x$OS.time,
              status = x$OS,
              marker = -x$RS,
              predict.time = 365 * AUC_time,
              method = "KM"
            )
          } else {
            print("Please provide the correct parameters for method")
          }
        }

        roc_1 <- cbind(round(risk.survivalROC$TP, 3), round(risk.survivalROC$FP, 3))


        roc_1 <- as.data.frame(roc_1)
        colnames(roc_1) <- c("TP", "FP")
        roc_1$AUC <- risk.survivalROC$AUC
        roc_1$HR <- HR
        return(roc_1)
      } else {
        roc_1 <- data.frame(
          TP = rep(0, nrow(x)),
          FP = rep(1, nrow(x)),
          AUC = rep(0.5, nrow(x)),
          HR = rep(1, nrow(x))
        )

        return(roc_1)
      }
    })






    return(roc.rs)
  }


  list_input_data <- lapply(list_input_data, function(x) {
    not.gene <- common_feature[which(!common_feature %in% colnames(x))]
    cons.mat <- as.data.frame(matrix(rep(0, length(not.gene) * nrow(x)), nrow = nrow(x), ncol = length(not.gene)))
    colnames(cons.mat) <- not.gene


    x <- cbind(x, cons.mat)

    x <- x[, common_feature]
    return(x)
  })




  list_input_data <- lapply(list_input_data, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, as.numeric)
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x[, c(2:3)] <- apply(x[, c(2:3)], 2, as.numeric)
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x <- x[!is.na(x$OS.time) & !is.na(x$OS), ]
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x <- x[x$OS.time > 0, ]
    return(x)
  })
  # use the mean replace the NA
  list_input_data <- lapply(list_input_data, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })


    return(x)
  })




  less.os.time <- max(list_input_data[[1]]$OS.time)

  for (i in names(list_input_data)) {
    print(i)
    if (max(list_input_data[[i]]$OS.time) < less.os.time) {
      less.os.time <- max(list_input_data[[i]]$OS.time)

      print(less.os.time)
    } else {
      less.os.time <- less.os.time
    }
  }

  print("Please wait a few minutes.")
  if (less.os.time > 365 * AUC_time) {
    model.name <- unique(sig.input$model)

    val_dd_list <- lapply(list_input_data, function(x) {
      x[, c("OS.time", "OS", common_feature)]
    })


    roc.table <- lapply(model.name, function(z) {
      coef.tab <- sig.input[sig.input$model == z, c("Coef", "symbol")]

      val_dd_list2 <- lapply(val_dd_list, function(x) {
        x[, c("OS.time", "OS", coef.tab$symbol)]
      })


      rs <- lapply(val_dd_list2, function(x) {
        cbind(x[, 1:2],
          RS = apply(as.data.frame(x[, -c(1:2)]), 1, function(x) {
            x %*% coef.tab$Coef
          })
        )
      })

      roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)




      return(roc.test)
    })
    names(roc.table) <- model.name


    return(roc.table)
  } else {
    print(paste0("The shortest overall survival time in the queue you provided is ", less.os.time, " days"))
    print("Please set a reasonable AUC_time that is less than this length of time!")
  }
}
