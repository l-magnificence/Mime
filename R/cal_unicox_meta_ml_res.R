#' Perform the meta-analysis based on univariable Cox regression result.
#'
#' @param input A data frame containing the univariable Cox regression result. The row names are the cohort names. The column names are "HR", "pvalue", "LCI", and "HCI".  LCI and HCI are matrices containing, respectively, the lower and higher bounds of the 95% confidence interval of each correlation coefficient. 
#'
#' @return The result of the meta-analysis.
#' @export
#'
#' @examples
cal_unicox_meta_ml_res <- function(input # 列名为 HR，pvalue ，LCI， HCI的一个队列，行名为队列名称，在这里主要是unicox的结果
) {
  if (identical(colnames(input), c("HR", "pvalue", "LCI", "HCI"))) {
    library(meta)
    input <- input %>%
      as.data.frame() %>%
      mutate(`Hazard Ratio(95%CI)` = paste(HR, "(", LCI, "-", HCI, ")", sep = ""))
    input$Group <- rownames(input)
    # 对不服从正态分布的HR对数转换
    lnhr <- log(input$HR)
    lnlci <- log(input$LCI)
    lnhci <- log(input$HCI)
    selnhr <- (lnhci - lnlci) / (2 * 1.96)
    metamodel <- metagen(
      TE = lnhr, seTE = selnhr,
      sm = "HR", # 待合并的效应量
      data = input, # 输入数据
      studlab = Group
    )
    return(metamodel)
  } else {
    print("Please provide the matrix with correct format!")
  }
}
