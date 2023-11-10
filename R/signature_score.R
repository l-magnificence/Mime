#' Calculating Signature Score using IOBR
#'
#' @param res.by.ML.Dev.Prog.Sig The output result of ML.Dev.Prog.Sig
#' @param inputmatrix.list A gene expression dataframe after log2(x+1) scaled. The first three of the column names are, in order, ID,OS.time, OS. Columns starting with the fourth are gene symbols. OS.time is a numeric variable in days. OS is a numeric variable containing 0, 1. 0: Alive, 1: Dead.
#' @param signature signature using in IOBR, such as signature::signature_tme, signature::signature_metabolism ,signature::signature_tumor, signature::signature_collection
#' @param estima_method estimating signature scores methods: "pca" "ssgsea" "zscore" "integration"
#' @param diff.test Whether to stratify patients into high and low-risk groups based on the median of risk scores, and obtain the differential expression of signatures between the two groups.
#' @param correlation.test Whether to perform correlation analysis on the specified signature scores.
#' @param cor.method Method using in correlation.test, such as "pearson", "kendall", "spearman"
#' @param cor.siganatue any signatures in the signature collection, such as CD_8_T_effector and DDR in IOBR::signature_collection
#'
#' @return a list containing signature scores in each cohorts and ML methods
#' @export
#'
#' @examples
#' test.sig <- signature_score(res, list_train_vali_Data)

signature_score <- function(res.by.ML.Dev.Prog.Sig, # ML.Dev.Prog.Sig, 函数计算结果
                            inputmatrix.list, # A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
                            signature = IOBR::signature_collection, # signature_tme signature_metabolism signature_tumor signature_collection
                            estima_method = "zscore", # "pca" "ssgsea" "zscore" "integration"
                            diff.test = T, # Whether to stratify patients into high and low-risk groups based on the median of risk scores, and obtain the differential expression of signatures between the two groups.
                            correlation.test = T, # Whether to perform correlation analysis on the specified signature scores.
                            cor.method = "pearson", # "pearson" "kendall" "spearman"
                            cor.siganatue = c("CD_8_T_effector", "DDR") # any signatures in the signature collection, such as CD_8_T_effector and DDR in IOBR::signature_collection
) {
  #### loading the packages ########
  if (T) {
    library(IOBR)
    library(dplyr)
    library(magrittr)
    library(data.table)
    library(readr)
    library(stats)

    Sys.setenv(LANGUAGE = "en") # 显示英文报错信息
    options(stringsAsFactors = FALSE) # 禁止chr转成factor
  }

  two_subtype_logFC_Cpr_wilcox <- function(
      expr = NULL,
      treat_list = NULL,
      ctrl_list = NULL) {
    # 自定义显示进程函数
    display.progress <- function(index, totalN, breakN = 20) {
      if (index %% ceiling(totalN / breakN) == 0) {
        cat(paste(round(index * 100 / totalN), "% ", sep = ""))
      }
    }

    meanA <- meanB <- p <- fc <- lgfc <- c() # 初始化向量
    for (k in 1:nrow(expr)) {
      # display.progress(index = k,totalN = nrow(expr)) #显示进程
      a <- as.numeric(expr[k, treat_list$ID])
      b <- as.numeric(expr[k, ctrl_list$ID])
      p <- c(p, wilcox.test(a, b, na.rm = T)$p.value) # 检验表达差异. wilcoxon
      # p <- c(p,t.test(a,b,na.rm=T,alternative = "two.sided" )$p.value) # 检验表达差异, t.student
      meanA <- c(meanA, mean(a)) # treat组基因k表达均值
      meanB <- c(meanB, mean(b)) # control组基因k表达均值
      fc <- c(fc, mean(a) / mean(b)) # 计算foldchange
      lgfc <- c(lgfc, log2(mean(a) / mean(b))) # 计算log2foldchange
    }
    fdr <- p.adjust(p, method = "fdr") # 校正p值

    # 生成差异表达结果，其中log2FoldChange, pvalue, padj模仿DESeq2结果格式
    # 这里采用简单的两样本wilcox.test检验寻找显著差异表达基因。
    degs <- data.frame(
      mean_treat = meanA,
      mean_ctrl = meanB,
      FoldChange = fc,
      log2FoldChange = lgfc,
      pvalue = p,
      padj = fdr,
      row.names = rownames(expr),
      stringsAsFactors = F
    )
    return(degs)
    # write.table(degs,paste0("two_samples_logFC","_degs.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  }

  message("--- Data preprocessing ---")

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })
    return(x)
  })
  sig_list <- list()
  risk_list <- list()
  cor.result.list <- list()
  print(paste0("For dividing patients into high and low-risk groups and estimating Signature enrichment analysis, we recommend using the z-score and ssGSEA methods."))
  for (i in 1:length(inputmatrix.list)) {
    train_data <- inputmatrix.list[[i]]
    test.matrix <- train_data %>%
      magrittr::set_rownames(.$ID) %>%
      select(-c(ID, OS.time, OS)) %>%
      t()
    # Signature score estimation
    if (estima_method %in% c("zscore", "ssGSEA")) {
      sig_tme <- IOBR::calculate_sig_score(
        pdata = NULL,
        eset = test.matrix,
        signature = signature,
        method = estima_method,
        mini_gene_count = 2
      )

      sig_tme <- t(tibble::column_to_rownames(sig_tme, var = "ID"))
      resultList <- list("sig_tme" = sig_tme)
    } else if (estima_method %in% c("pca", "integration")) {
      ### For the PCA method, it is necessary to remove genes with constant expression values across all samples.
      constant_gene <- which(apply(test.matrix, 1, function(x) length(unique(x)) == 1))
      constant_gene_names <- rownames(test.matrix)[constant_columns]
      print(paste0("Some genes have constant values. For the PCA or integration method, it is necessary to remove genes with constant expression values across all samples."))

      test.matrix <- test.matrix[-constant_gene, ]

      sig_tme <- IOBR::calculate_sig_score(
        pdata = NULL,
        eset = test.matrix,
        signature = signature_collection,
        method = estima_method,
        mini_gene_count = 2
      )

      sig_tme <- t(column_to_rownames(sig_tme, var = "ID"))
      resultList <- list("sig_tme" = sig_tme)
    }
    resultList_1 <- list(resultList)
    names(resultList_1) <- names(inputmatrix.list)[i]
    sig_list <- append(sig_list, resultList_1)

    if (diff.test == T) {
      print("When Log2foldchange > 0, the signature highly expressed in High Risk Group")
      for (j in 1:length(res.by.ML.Dev.Prog.Sig$riskscore)) {
        # Divide the patients into high and low-risk groups, and calculate the difference in signature scores by IOBR for each group
        tmp <- res.by.ML.Dev.Prog.Sig$riskscore[[j]][[names(inputmatrix.list)[i]]]
        if (max(tmp$RS) == median(tmp$RS)) {
          next
        } else {
          tmp$riskgroup <- ifelse(tmp$RS > median(tmp$RS), "High", "Low")

          tmp_highrisk <- tmp %>%
            filter(riskgroup == "High") %>%
            select(ID) %>%
            as.list()

          tmp_lowrisk <- tmp %>%
            filter(riskgroup == "Low") %>%
            select(ID) %>%
            as.list()

          ## Calculate the differential signature between the high-risk and low-risk groups
          diff_sig <- two_subtype_logFC_Cpr_wilcox(expr = sig_tme, treat_list = tmp_highrisk, ctrl_list = tmp_lowrisk)

          riskresult <- list(list("Signature_data" = sig_tme, "diffenential_siganature_wilcox" = diff_sig, "Patient_risk" = tmp))

          print(paste0("diff.test", " in ", names(inputmatrix.list)[i], ": ", names(res.by.ML.Dev.Prog.Sig$riskscore[j]), " Processed Successful"))
        }
        names(riskresult) <- names(res.by.ML.Dev.Prog.Sig$riskscore)[j]
        risk_list <- append(risk_list, riskresult)
      }
    }

    if (correlation.test == T) {
      cor.sig <- c(cor.siganatue)
      sig_tme_cor <- t(sig_tme[cor.sig, ])

      sig_tme_cor <- stats::cor(sig_tme_cor,
        method = cor.method
      )
      # library(ggcorrplot)
      # ggcorrplot::ggcorrplot(sig_tme_cor,                  # Draw ggcorrplot with p-values
      #                        p.mat = cor_test_mat,
      #                        lab = TRUE,
      #                        type = "lower",
      #                        colors = c("#6395C7", "white", "#E06EAD"))
      cor.result <- list("correlation_result" = sig_tme_cor)
      names(cor.result) <- names(inputmatrix.list)[i]
      cor.result.list <- append(cor.result.list, cor.result)
      print(paste0("correlation.test", " in ", names(inputmatrix.list)[i], ": ", names(res.by.ML.Dev.Prog.Sig$riskscore[j]), " Processed Successfu"))
    }

    print(paste0(names(inputmatrix.list)[i], " Processed Successfully"))
  }

  signature_result_all <- list(sig_list, risk_list)
  names(signature_result_all) <- c("signature_list", "differential_sig_list")

  return(signature_result_all)
}
