#' Performing the Univariable Cox Regression Analysis for the Optimal Model
#'
#' @param res.by.ML.Dev.Prog.Sig The output result of ML.Dev.Prog.Sig
#' @param optimal.model  One model name you choose. For example, 'RSF'.
#' @param type c('categorical') or c('continuous). The type for the univariable Cox Regression. The 'categorical' means that the patients are divided into two subgroups based on the median of the risk score. This means conducting  univariable Cox regression with the categorical subgroup. The 'continuous' means calculating the univariable Cox regression with the continuous risk score.      
#'
#' @return A list of the univariable Cox regression result of the model.
#' @export
#'
#' @examples
cal_unicox_ml_res <- function(res.by.ML.Dev.Prog.Sig = NULL, # ML.Dev.Prog.Sig, 函数计算结果，内置包含有riskscore table
                              optimal.model, # 选择要计算的model
                              type = "categorical" # c('categorical'),c('continuous) 单因素回归 采用的变量，连续变量和分类变量，如果是分类变量，则按照中位值区分High and Low # 默认为分类变量
) {
  library(tidyverse)
  library(survival)
  ### 选定一个model 计算risk score 的单因素回归分析

  optimal.model <- optimal.model

  mm <- res.by.ML.Dev.Prog.Sig[["riskscore"]][[optimal.model]]


  unicox.cal <- lapply(mm, function(x) {
    tmp <- x

    if (type == "categorical") {
      # 分类变量
      tmp$Group <- ifelse(tmp$RS > median(tmp$RS), "High", "Low")

      tmp$Group <- factor(tmp$Group, levels = c("Low", "High"))
      cox <- coxph(Surv(OS.time, OS) ~ Group, data = tmp)
    } else if (type == "continuous") {
      # 连续变量

      cox <- coxph(Surv(OS.time, OS) ~ RS, data = tmp)
    } else {
      print("Please provided the correct parameters")
    }




    coxSummary <- summary(cox)

    unicox <- c(
      as.numeric(coxSummary$coefficients[, "exp(coef)"])[1],
      as.numeric(coxSummary$coefficients[, "Pr(>|z|)"])[1],
      as.numeric(coxSummary$conf.int[, 3][1]),
      as.numeric(coxSummary$conf.int[, 4][1])
    )

    return(unicox)
  })

  unicox.rs.res <- unicox.cal %>% do.call(rbind, .)
  colnames(unicox.rs.res) <- c("HR", "pvalue", "LCI", "HCI")


  return(unicox.rs.res)
}
