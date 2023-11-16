#' Calculate AUC scores of Machine Learning Models in all data
#'
#' @param res.by.ML.Dev.Prog.Sig the results of function ML.Dev.Prog.Sig
#' @param train_data the training data using in ML.Dev.Prog.Sig
#' @param inputmatrix.list A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
#' @param mode Choose MF models: all, single, double
#' @param AUC_time  c(1,2,3,4,5,6,7,······), for 1 year, 2 years, 3 years......We recommend using the shortest survival time among all queues.
#' @param single_ml If the mode is set to "single", you must fill in the following models: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso").
#' @param double_ml1 If the mode is set to "double", you need to fill in the modeling methods here: c('RSF', "StepCox","CoxBoost","Lasso").
#' @param double_ml2 If the mode is set to "double", you need to fill in the modeling methods here: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
#' @param alpha_for_Enet 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
#' @param direction_for_stepcox  c("both", "backward", "forward")
#' @param auc_cal_method KM, NNE
#'
#' @return AUC Calculation Result
#' @export
#'
#' @examples
#' all.auc.1y <- cal_AUC_ml_res(
#'   res.by.ML.Dev.Prog.Sig = res,
#'   train_data = list_train_vali_Data[["TCGA"]],
#'   inputmatrix.list = list_train_vali_Data, mode = "all",
#'   AUC_time = 1,
#'   auc_cal_method = "KM"
#' )
#'
cal_AUC_ml_res <- function(res.by.ML.Dev.Prog.Sig = NULL, # ML.Dev.Prog.Sig, 函数计算结果
                           train_data, # ML.Dev.Prog.Sig 中的训练集
                           inputmatrix.list, # A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
                           mode = NULL, # all, single, double
                           AUC_time = 1, # c(1,2,3,4,5,6,7,······),1年， 2年， 3年····。这里建议使用所有队列中最短的生存时间。
                           single_ml = NULL, # 如果 mode 为single 则必须要填 c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           double_ml1 = NULL, #  如果mode 为double 则这里需要填写建模的 方法。c('RSF', "StepCox","CoxBoost","Lasso")
                           double_ml2 = NULL, # c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           alpha_for_Enet = NULL, # 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
                           direction_for_stepcox = NULL, #  c("both", "backward", "forward")
                           auc_cal_method = "KM" # KM, NNE 默认为KM
) {
  #### loading the packages ########
  if (T) {
    library(survival)
    library(randomForestSRC)
    library(glmnet)
    library(plsRcox)
    library(superpc)
    library(gbm)
    library(CoxBoost)
    library(survivalsvm)
    library(dplyr)
    library(tibble)
    library(BART)
    library(miscTools)
    library(compareC)
    library(ggplot2)
    library(ggsci)
    library(tidyr)
    library(ggbreak)
    library(mixOmics)
    library(data.table)
    library(survivalROC)



    Sys.setenv(LANGUAGE = "en") # 显示英文报错信息
    options(stringsAsFactors = FALSE) # 禁止chr转成factor
  }

  message("--- Data preprocessing ---")

  ####  Data preprocessing #####
  # Replace '-' in column names with '.'
  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    return(x)
  })


  colnames(train_data) <- gsub("-", ".", colnames(train_data))


  sig.gene <- res.by.ML.Dev.Prog.Sig$Sig.genes


  # Matching candidate genes to genes in each cohort
  common_feature <- c("ID", "OS.time", "OS", sig.gene)
  common_feature <- gsub("-", ".", common_feature)

  for (i in names(inputmatrix.list)) {
    common_feature <- intersect(common_feature, colnames(inputmatrix.list[[i]]))
  }


  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x <- x[, common_feature]
    return(x)
  })

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, as.numeric)
    return(x)
  })

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
    return(x)
  })

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, c(2:3)] <- apply(x[, c(2:3)], 2, as.numeric)
    return(x)
  })

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x <- x[!is.na(x$OS.time) & !is.na(x$OS), ]
    return(x)
  })

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x <- x[x$OS.time > 0, ]
    return(x)
  })
  # use the mean replace the NA
  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })


    return(x)
  })



  less.os.time <- max(inputmatrix.list[[1]]$OS.time)

  for (i in names(inputmatrix.list)) {
    print(i)
    if (max(inputmatrix.list[[i]]$OS.time) < less.os.time) {
      less.os.time <- max(inputmatrix.list[[i]]$OS.time)

      print(less.os.time)
    } else {
      less.os.time <- less.os.time
    }
  }


  if (less.os.time > 365 * AUC_time) {
    #### calculating the risk score of the signatures #########

    AUC.rs <- list()

    if (mode %in% c("all", "single", "double") & identical(colnames(inputmatrix.list[[1]])[1:3], c("ID", "OS.time", "OS"))) {
      feature.a <- res.by.ML.Dev.Prog.Sig$Sig.genes
      feature.a <- gsub("-", ".", feature.a)

      num.a <- length(feature.a)

      val_dd_list <- lapply(inputmatrix.list, function(x) {
        x[, c("OS.time", "OS", feature.a)]
      })
      val_dd_list2 <- val_dd_list


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
        })






        return(roc.rs)
      }



      est_dd2 <- train_data[, c("OS.time", "OS", feature.a)]
      data <- list(
        x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
        censoring.status = est_dd2$OS,
        featurenames = colnames(est_dd2)[-c(1, 2)]
      )


      if (mode == "all") {
        feature.ab <- feature.a
        for (i in names(inputmatrix.list)) {
          feature.ab <- intersect(feature.ab, colnames(inputmatrix.list[[i]]))
        }

        num.b <- length(feature.ab)



        if (num.b == num.a) {
          ml.names <- names(res.by.ML.Dev.Prog.Sig$ml.res)



          for (i in ml.names[c(which(ml.names %in% c("RSF")), grep("+ RSF", ml.names, fixed = T))]) {
            print(i)


            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- fit[["xvar.names"]]


            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)




            AUC.rs[[i]] <- roc.test
          }


          for (i in ml.names[c(
            which(ml.names %in% c("survival - SVM")), grep("+ survival-SVM", ml.names, fixed = T)
          )]) {
            print(i)


            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- fit[["var.names"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[i]] <- roc.test
          }



          for (i in ml.names[c(which(ml.names %in% c("CoxBoost")), grep("+ CoxBoost", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- fit[["xnames"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 2], type = "lp")))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }


          for (i in ml.names[c(grep("Enet[", ml.names, fixed = T), grep("+ Enet", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[i]] <- roc.test
          }


          for (i in ml.names[c(which(ml.names %in% c("GBM")), grep("+ GBM", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]][[1]]
            best <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]][[2]]

            feature.ac <- fit[["var.names"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }

          for (i in ml.names[c(which(ml.names %in% c("Lasso")), grep("+ Lasso", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }



          for (i in ml.names[c(which(ml.names %in% c("plsRcox")), grep("+ plsRcox", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- colnames(fit[["dataX"]])
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }




          for (i in ml.names[which(ml.names %in% c("Ridge"))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]][[1]]
            cv.fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]][[2]]

            feature.ac <- fit[["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = cv.fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }

          for (i in ml.names[c(grep("+ Ridge", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[i]] <- roc.test
          }



          for (i in ml.names[c(which(ml.names %in% c("StepCox[both]", "StepCox[backward]", "StepCox[forward]")), grep("+ StepCox", ml.names, fixed = T))]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]]

            feature.ac <- names(fit[["coefficients"]])
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })

            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }


          for (i in ml.names[grep("SuperPC", ml.names, fixed = T)]) {
            print(i)

            fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]][[1]]
            cv.fit <- res.by.ML.Dev.Prog.Sig[["ml.res"]][[i]][[2]]


            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ab)]
            })

            rs <- lapply(val_dd_list, function(w) {
              test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 2)])
              ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
              rr <- as.numeric(ff$v.pred)
              rr2 <- cbind(w[, 1:2], RS = rr)
              return(rr2)
            })

            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[i]] <- roc.test
          }

          return(AUC.rs)
        } else {
          miss.gene <- feature.a[which(!feature.a %in% feature.ab)]

          print(paste0("Sorry, There are some genes in not matched,
                     meaning that some of the genes appear in the model
                     but do not exist inside certain cohorts."))
        }
      } else if (mode == "single") {
        feature.ab <- feature.a
        for (i in names(inputmatrix.list)) {
          feature.ab <- intersect(feature.ab, colnames(inputmatrix.list[[i]]))
        }

        num.b <- length(feature.ab)

        if (num.b == num.a) {
          if (single_ml == "RSF") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            print(i)

            feature.ac <- fit[["xvar.names"]]


            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (single_ml == "Enet") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (single_ml == "StepCox") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- names(fit[["coefficients"]])
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })

            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (single_ml == "CoxBoost") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- fit[["xnames"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 2], type = "lp")))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)
            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (single_ml == "plsRcox") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- colnames(fit[["dataX"]])
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (single_ml == "superpc") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$SuperPC$fit
            cv.fit <- res.by.ML.Dev.Prog.Sig$ml.res$SuperPC$cv.fit
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ab)]
            })

            rs <- lapply(val_dd_list, function(w) {
              test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 2)])
              ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
              rr <- as.numeric(ff$v.pred)
              rr2 <- cbind(w[, 1:2], RS = rr)
              return(rr2)
            })

            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)
            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (single_ml == "GBM") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$GBM$fit
            best <- res.by.ML.Dev.Prog.Sig$ml.res$GBM$best
            feature.ac <- fit[["var.names"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (single_ml == "survivalsvm") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- fit[["var.names"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (single_ml == "Ridge") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$Ridge$fit
            cv.fit <- res.by.ML.Dev.Prog.Sig$ml.res$Ridge$cv.fit

            feature.ac <- fit[["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = cv.fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (single_ml == "Lasso") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$Lasso
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)
            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }
        } else {
          miss.gene <- feature.a[which(!feature.a %in% feature.ab)]

          print(paste0("Sorry, There are some genes in not matched,
                     meaning that some of the genes appear in the model
                     but do not exist inside certain cohorts."))
        }
      } else if (mode == "double") {
        feature.ab <- feature.a
        for (i in names(inputmatrix.list)) {
          feature.ab <- intersect(feature.ab, colnames(inputmatrix.list[[i]]))
        }

        num.b <- length(feature.ab)

        if (num.b == num.a) {
          if (double_ml2 == "RSF") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            print(i)


            feature.ac <- fit[["xvar.names"]]


            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (double_ml2 == "Enet") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (double_ml2 == "StepCox") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- names(fit[["coefficients"]])
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })

            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (double_ml2 == "CoxBoost") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- fit[["xnames"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 2], type = "lp")))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)
            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (double_ml2 == "plsRcox") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- colnames(fit[["dataX"]])
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (double_ml2 == "superpc") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$SuperPC$fit
            cv.fit <- res.by.ML.Dev.Prog.Sig$ml.res$SuperPC$cv.fit
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ab)]
            })

            rs <- lapply(val_dd_list, function(w) {
              test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 2)])
              ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
              rr <- as.numeric(ff$v.pred)
              rr2 <- cbind(w[, 1:2], RS = rr)
              return(rr2)
            })

            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)
            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }


          if (double_ml2 == "GBM") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$GBM$fit
            best <- res.by.ML.Dev.Prog.Sig$ml.res$GBM$best
            feature.ac <- fit[["var.names"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (double_ml2 == "survivalsvm") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res[[1]]
            feature.ac <- fit[["var.names"]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)


            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (double_ml2 == "Ridge") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$Ridge$fit
            cv.fit <- res.by.ML.Dev.Prog.Sig$ml.res$Ridge$cv.fit

            feature.ac <- fit[["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })


            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = cv.fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)

            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }

          if (double_ml2 == "Lasso") {
            fit <- res.by.ML.Dev.Prog.Sig$ml.res$Lasso
            feature.ac <- fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]]
            val_dd_list2 <- lapply(val_dd_list, function(x) {
              x[, c("OS.time", "OS", feature.ac)]
            })



            rs <- lapply(val_dd_list2, function(x) {
              cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
            })
            roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)
            AUC.rs[[names(res.by.ML.Dev.Prog.Sig$ml.res)]] <- roc.test

            return(AUC.rs)
          }
        } else {
          miss.gene <- feature.a[which(!feature.a %in% feature.ab)]

          print(paste0("Sorry, There are some genes in not matched,
                     meaning that some of the genes appear in the model
                     but do not exist inside certain cohorts."))
        }
      }
    } else {
      print("Please refer to the sample data and process to improve the relevant parameters and preprocessing files")
    }
  } else {
    print(paste0("The shortest overall survival time in the queue you provided is ", less.os.time, " days"))
    print("Please set a reasonable AUC_time that is less than this length of time!")
  }
}
