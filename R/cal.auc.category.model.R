#' Calculating the area under the curve after developing the category predictive model
#'
#' @param res.by.ML.Dev.Pred.Category.Sig  Output of function ML.Dev.Pred.Category.Sig
#' @param cohort.for.cal A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#'
#' @return A data frame containing the AUC of each predictive model.
#' @export
#'
#' @examples
cal.auc.category.model <- function(res.by.ML.Dev.Pred.Category.Sig, ### 函数计算结果
                                   cohort.for.cal # 队列要求第一列为ID,第二列为分类变量Var, 值为Y或者N, 从第三列开始为基因，表达矩阵经过了log2(x+1)处理
) {
  if (T) {
    library(stringr)
    library(gridExtra)
    library(future)
    library(sva)
    library(e1071)
    library(pROC)
    library(ROCit)
    library(caret)
    library(doParallel)
    library(cancerclass)
    library(dplyr)
  }

  sig <- res.by.ML.Dev.Pred.Category.Sig$sig.gene

  rownames(cohort.for.cal) <- cohort.for.cal$ID
  colnames(cohort.for.cal) <- gsub("-", ".", colnames(cohort.for.cal))
  sig <- gsub("-", ".", sig)


  if (all(is.element(sig, colnames(cohort.for.cal))) & identical(colnames(cohort.for.cal)[1:2], c("ID", "Var"))) {
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", "Y"))

    ls_model <- res.by.ML.Dev.Pred.Category.Sig$model
    auc <- lapply(ls_model, function(model.tune) {
      if (class(model.tune) == "predictor") {
        pData <- data.frame(class = validation$Var, sample = rownames(validation), row.names = rownames(validation))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(validation[, -1])
        Sig.Exp.test <- ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)
        prediction <- predict(model.tune, Sig.Exp.test, "N", ngenes = nrow(Sig.Exp), dist = "cor")
        roc <- roc(
          response = prediction@prediction[, "class_membership"],
          predictor = as.numeric(prediction@prediction[, "z"])
        )
        roc_result <- coords(roc, "best")
        auc <- data.frame(ROC = as.numeric(roc$auc), Sens = roc_result$sensitivity[1], Spec = roc_result$specificity[1])
      } else {
        prob <- predict(model.tune, validation[, -1], type = "prob")
        pre <- predict(model.tune, validation[, -1])
        test_set <- data.frame(obs = validation$Var, N = prob[, "N"], Y = prob[, "Y"], pred = pre)
        auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
      }

      return(auc)
    }) %>% base::do.call(rbind, .)

    rownames(auc) <- names(ls_model)

    return(auc)
  } else {
    print("Please provide the correct parameters and the cohorts with correct format")
  }
}
