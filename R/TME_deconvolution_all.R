#' Title
#'
#' @param inputmatrix.list
#' @param deconvolution_methods
#' @param microarray_names
#'
#' @return
#' @export
#'
#' @examples
TME_deconvolution_all <- function(inputmatrix.list, # A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
                                  deconvolution_methods = "ALL", # Deconvolution Methods in IOBR::deconvo_tme
                                  microarray_names = "none" # Please tell us which datasets are microarray, use the names of elements in inputmatrix.list. such as c("CGGA.array", "GSE108474", "GSE16011", "GSE43289", "GSE7696") if none, enter "none".
) {
  #### loading the packages ########
  if (T) {
    library(IOBR)
    library(dplyr)
    library(magrittr)
    library(data.table)
    library(readr)
    Sys.setenv(LANGUAGE = "en") # 显示英文报错信息
    options(stringsAsFactors = FALSE) # 禁止chr转成factor
  }

  message("--- Data preprocessing ---")
  # Deconvolution by IOBR
  TME_deconvolution <- function(test.matrix,
                                deconvolution_methods = "ALL",
                                arrays = F) {
    if (deconvolution_methods == "ALL") {
      cibersort <- deconvo_tme(eset = test.matrix, method = "cibersort", arrays = arrays, perm = 100)
      epic <- deconvo_tme(eset = test.matrix, method = "epic", arrays = arrays)
      mcp <- deconvo_tme(eset = test.matrix, method = "mcpcounter")
      xcell <- deconvo_tme(eset = test.matrix, method = "xcell", arrays = arrays)
      estimate <- deconvo_tme(eset = test.matrix, method = "estimate")
      estimate$ID <- gsub("-", ".", estimate$ID)
      timer <- deconvo_tme(eset = test.matrix, method = "timer", group_list = rep("stad", dim(test.matrix)[2]))
      quantiseq <- deconvo_tme(eset = test.matrix, tumor = TRUE, arrays = arrays, scale_mrna = TRUE, method = "quantiseq")
      ips <- deconvo_tme(eset = test.matrix, method = "ips", plot = FALSE)

      tme_combine <- cibersort %>%
        inner_join(., mcp, by = "ID") %>%
        inner_join(., xcell, by = "ID") %>%
        inner_join(., epic, by = "ID") %>%
        inner_join(., estimate, by = "ID") %>%
        inner_join(., timer, by = "ID") %>%
        inner_join(., quantiseq, by = "ID") %>%
        inner_join(., ips, by = "ID")

      resultList <- list("tme_combine" = tme_combine)
      return(resultList)
    } else if (deconvolution_methods %in% tme_deconvolution_methods) {
      if ("cibersort" == deconvolution_methods) {
        cibersort <- deconvo_tme(eset = test.matrix, method = "cibersort", arrays = arrays, perm = 100)
        resultList <- list("cibersort" = cibersort)
        return(resultList)
      }
      if ("epic" == deconvolution_methods) {
        epic <- deconvo_tme(eset = test.matrix, method = "epic", arrays = arrays)
        resultList <- list("epic" = epic)
        return(resultList)
      }
      if ("mcpcounter" == deconvolution_methods) {
        mcp <- deconvo_tme(eset = test.matrix, method = "mcpcounter")
        resultList <- list("mcpcounter" = mcp)
        return(resultList)
      }
      if ("xcell" == deconvolution_methods) {
        xcell <- deconvo_tme(eset = test.matrix, method = "xcell", arrays = arrays)
        resultList <- list("xcell" = xcell)
        return(resultList)
      }
      if ("estimate" == deconvolution_methods) {
        estimate <- deconvo_tme(eset = test.matrix, method = "estimate")
        estimate$ID <- gsub("-", ".", estimate$ID)
        resultList <- list("estimate" = estimate)
        return(resultList)
      }
      if ("timer" == deconvolution_methods) {
        timer <- deconvo_tme(eset = test.matrix, method = "timer", group_list = rep("stad", dim(test.matrix)[2]))
        resultList <- list("timer" = timer)
        return(resultList)
      }
      if ("quantiseq" == deconvolution_methods) {
        quantiseq <- deconvo_tme(eset = test.matrix, tumor = TRUE, arrays = arrays, scale_mrna = TRUE, method = "quantiseq")
        resultList <- list("quantiseq" = quantiseq)
        return(resultList)
      }
      if ("ips" == deconvolution_methods) {
        ips <- deconvo_tme(eset = test.matrix, method = "ips", plot = FALSE)
        resultList <- list("ips" = ips)
        return(resultList)
      }
    } else {
      print("Please provide the correct parameters for deconvolution method")
    }
  }

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })
    return(x)
  })
  tme_decon_list <- list()
  # cat("Datasets names:", paste(names(inputmatrix.list), collapse = ", "), "\n")
  # cat("Please tell us which datasets are microarray, if none, skip the step\n")  # 提供一个示例


  selected_columns <- microarray_names

  for (i in 1:length(inputmatrix.list)) {
    train_data <- inputmatrix.list[[i]]
    test.matrix <- train_data %>%
      magrittr::set_rownames(.$ID) %>%
      select(-c(ID, OS.time, OS)) %>%
      t()

    tryCatch(
      {
        if (selected_columns == "none") {
          resultList <- TME_deconvolution(test.matrix = test.matrix, deconvolution_methods = deconvolution_methods, arrays = F)
        } else if (all(selected_columns %in% names(inputmatrix.list))) {
          # 将用户输入的名字拆分成一个字符向量
          if (names(inputmatrix.list)[i] %in% selected_columns) {
            # 确保用户输入的名字都存在
            resultList <- TME_deconvolution(test.matrix = test.matrix, deconvolution_methods = deconvolution_methods, arrays = T)
          } else {
            resultList <- TME_deconvolution(test.matrix = test.matrix, deconvolution_methods = deconvolution_methods, arrays = F)
          }
        } else {
          cat("Invalid dataset name(s). Please try again.")
        }
        print("Success")
      },
      error = function(e) {
        print(paste("An error occurred:", conditionMessage(e)))
        resultList <- NULL
      },
      finally = {
        resultList_1 <- list(resultList)
        names(resultList_1) <- names(inputmatrix.list)[i]
        tme_decon_list <- append(tme_decon_list, resultList_1)
      }
    )
  }
  return(tme_decon_list)
}
