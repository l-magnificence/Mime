#' Screening out the core variables for the prognosis with the machine learning algorithms
#' 
#' A function can be used for screening out the core features from the given candidate genes with eight machine learning algorithms. 
#' 
#' @param InputMatrix A gene expression dataframe after log2(x+1) scaled. The first three of the column names are, in order, ID,OS.time, OS. Columns starting with the fourth are gene symbols. OS.time is a numeric variable in days. OS is a numeric variable containing 0, 1. 0: Alive, 1: Dead. 
#' @param candidate_genes The input genes, that you want to screen out from, for identifying the core features.
#' @param mode  We provide three modes including 'all', 'single', and 'all_without_SVM'. The 'all' mode means using all eight methods for selecting. The 'single' mode means using only one method for running. Since SVM takes so much time, we're singling him out. The 'all_without_SVM' mode means the other seven methods used for selecting. 
#' @param seed  The seed. You can set it as any number. For example, 5201314.
#' @param single_ml The one method from the eight methods including "RSF", "Enet", "Boruta", "Xgboost", "SVM-REF", "Lasso", "CoxBoost', 'StepCox'.
#' @param nodesize The node size parameter for the RSF method. The default is 5. You can try another positive integer. For example, 10,15,20, etc. 
#'
#' @return A data frame including the methods and the core genes screened by the corresponding algorithm.
#' @export
#'
#' @examples
ML.Corefeature.Prog.Screen <- function(InputMatrix, ### 第一列ID,第二列OS.time, (day), 第三列 OS, (0/1), 从第四列开始是基因表达矩阵，经过了log2(x+1)
                                       candidate_genes,
                                       mode = NULL, # all, single,all_without_SVM
                                       seed = NULL,
                                       single_ml = NULL, # c("RSF", "Enet", "Boruta","Xgboost","SVM-REF","Lasso","CoxBoost','StepCox')
                                       nodesize = 5) {
  ### Screen out the core features via the multiple machine leaning algorithms
  ### loading the packages ####

  if (T) {
    Biocductor_packages <- c(
      "tidyverse",
      "scales",
      "Hmisc",
      "survival",
      "randomForestSRC",
      "glmnet",
      "plsRcox",
      "CoxBoost",
      "survivalsvm",
      "dplyr",
      "tibble",
      "BART",
      "miscTools",
      "compareC",
      "tidyr",
      "mixOmics",
      "data.table",
      "pbapply",
      "e1071",
      "Boruta",
      "xgboost",
      "Ckmeans.1d.dp",
      "Matrix"
    )

    lapply(Biocductor_packages, function(x) {
      library(x,
        character.only = T
        # ,  lib.loc = "/export/bioinfo-team/home/xiongzj/R/x86_64-pc-linux-gnu-library/4.1"
      )
    })
  }

  ### laoding the function ####

  if (T) {
    svmRFE.wrap <- function(test.fold, X, ...) {
      # Wrapper to run svmRFE function while omitting a given test fold
      train.data <- X[-test.fold, ]
      test.data <- X[test.fold, ]

      # Rank the features
      features.ranked <- svmRFE(train.data, ...)

      return(list(feature.ids = features.ranked, train.data.ids = row.names(train.data), test.data.ids = row.names(test.data)))
    }

    svmRFE <- function(X, k = 1, halve.above = 5000) {
      # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
      n <- ncol(X) - 1

      # Scale data up front so it doesn't have to be redone each pass
      cat("Scaling data...")
      X[, -1] <- scale(X[, -1])
      cat("Done!\n")
      flush.console()

      pb <- txtProgressBar(1, n, 1, style = 3)

      i.surviving <- 1:n
      i.ranked <- n
      ranked.list <- vector(length = n)

      # Recurse through all the features
      while (length(i.surviving) > 0) {
        if (k > 1) {
          # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)
          folds <- rep(1:k, len = nrow(X))[sample(nrow(X))]
          folds <- lapply(1:k, function(x) which(folds == x))

          # Obtain weights for each training set
          w <- lapply(folds, getWeights, X[, c(1, 1 + i.surviving)])
          w <- do.call(rbind, w)

          # Normalize each weights vector
          w <- t(apply(w, 1, function(x) x / sqrt(sum(x^2))))

          # Compute ranking criteria
          v <- w * w
          vbar <- apply(v, 2, mean)
          vsd <- apply(v, 2, sd)
          c <- vbar / vsd
        } else {
          # Only do 1 pass (i.e. regular SVM-RFE)
          w <- getWeights(NULL, X[, c(1, 1 + i.surviving)])
          c <- w * w
        }

        # Rank the features
        ranking <- sort(c, index.return = T)$ix
        if (length(i.surviving) == 1) {
          ranking <- 1
        }

        if (length(i.surviving) > halve.above) {
          # Cut features in half until less than halve.above
          nfeat <- length(i.surviving)
          ncut <- round(nfeat / 2)
          n <- nfeat - ncut

          cat("Features halved from", nfeat, "to", n, "\n")
          flush.console()

          pb <- txtProgressBar(1, n, 1, style = 3)
        } else {
          ncut <- 1
        }

        # Update feature list
        ranked.list[i.ranked:(i.ranked - ncut + 1)] <- i.surviving[ranking[1:ncut]]
        i.ranked <- i.ranked - ncut
        i.surviving <- i.surviving[-ranking[1:ncut]]

        setTxtProgressBar(pb, n - length(i.surviving))
        flush.console()
      }

      close(pb)

      return(ranked.list)
    }

    getWeights <- function(test.fold, X) {
      # Fit a linear SVM model and obtain feature weights
      train.data <- X
      if (!is.null(test.fold)) train.data <- X[-test.fold, ]

      svmModel <- svm(train.data[, -1], train.data[, 1],
        cost = 10, cachesize = 500,
        scale = F, type = "C-classification", kernel = "linear"
      )

      t(svmModel$coefs) %*% svmModel$SV
    }

    WriteFeatures <- function(results, input, save = T, file = "features_ranked.txt") {
      # Compile feature rankings across multiple folds
      featureID <- sort(apply(sapply(results, function(x) sort(x$feature, index.return = T)$ix), 1, mean), index = T)$ix
      avg.rank <- sort(apply(sapply(results, function(x) sort(x$feature, index.return = T)$ix), 1, mean), index = T)$x
      feature.name <- colnames(input[, -1])[featureID]
      features.ranked <- data.frame(FeatureName = feature.name, FeatureID = featureID, AvgRank = avg.rank)
      if (save == T) {
        write.table(features.ranked, file = file, quote = F, row.names = F)
      } else {
        features.ranked
      }
    }

    FeatSweep.wrap <- function(i, results, input) {
      # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
      svm.list <- lapply(results, function(x) {
        e1071::tune(svm,
          train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
          train.y = input[x$train.data.ids, 1],
          validation.x = input[x$test.data.ids, 1 + x$feature.ids[1:i]],
          validation.y = input[x$test.data.ids, 1],
          # Optimize SVM hyperparamters
          ranges = e1071::tune(svm,
            train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
            train.y = input[x$train.data.ids, 1],
            ranges  = list(gamma = 2^(-12:0), cost = 2^(-6:6))
          )$best.par,
          tunecontrol = tune.control(sampling = "fix")
        )$perf
      })

      error <- mean(sapply(svm.list, function(x) x$error))
      return(list(svm.list = svm.list, error = error))
    }


    PlotErrors <- function(errors, errors2 = NULL, no.info = 0.5,
                           ylim = range(c(errors, errors2), na.rm = T),
                           xlab = "Number of Features", ylab = "5 x CV Error") {
      # Makes a plot of average generalization error vs. number of top features
      AddLine <- function(x, col = "dodgerblue") {
        lines(which(!is.na(errors)), na.omit(x), col = col, lwd = 3)
        points(which.min(x), min(x, na.rm = T), col = "firebrick3")
        text(which.min(x), min(x, na.rm = T), paste(
          which.min(x), "-",
          format(min(x, na.rm = T), dig = 3)
        ), pos = 2, col = "red", cex = 1.15)
      }

      plot(errors, type = "n", ylim = ylim, xlab = xlab, ylab = ylab)
      AddLine(errors)
      if (!is.null(errors2)) AddLine(errors2, "gray30")
      abline(h = no.info, lty = 2)
    }


    Plotaccuracy <- function(errors, errors2 = NULL, no.info = 0.5,
                             ylim = range(c(errors, errors2), na.rm = T),
                             xlab = "Number of Features", ylab = "5 x CV Accuracy") {
      # Makes a plot of average generalization error vs. number of top features
      AddLine <- function(x, col = "dodgerblue") {
        lines(which(!is.na(errors)), na.omit(x), col = col, lwd = 3)
        points(which.max(x), max(x, na.rm = T), col = "firebrick3")
        text(which.max(x), max(x, na.rm = T), paste(
          which.max(x), "-",
          format(max(x, na.rm = T), dig = 3)
        ), pos = 2, col = "red", cex = 1.15)
      }

      plot(errors, type = "n", ylim = ylim, xlab = xlab, ylab = ylab)
      AddLine(errors)
      if (!is.null(errors2)) AddLine(errors2, "gray30")
      abline(h = no.info, lty = 2)
    }


    SigKMcox <- function(gene_list,
                         inputSet,
                         KM_pcutoff # KM的筛选阈值
    ) {
      print("Starting the data preprocess")
      ############### 数据预处理#######

      print("Rejecting a null value")

      library(survival)
      # 将空值的基因变成0
      # table(is.na(inputSet))
      inputSet[is.na(inputSet)] <- 0
      # table(is.na(inputSet))
      inputSet <- inputSet %>% as.data.frame()
      inputSet$OS.time <- as.numeric(inputSet$OS.time) # 时间数值化
      inputSet <- inputSet[inputSet$OS.time > 0, ] # 剔除时间为0

      # print("Correcting gene set")
      # #将基因集矫正 策略为转化从ALIAS为ENTREZID， 然后再转化为SYMBOL,避免别名
      # gene.df <- bitr(gene_list, fromType = "ALIAS",
      #                 toType = c("ENTREZID"),
      #                 OrgDb = org.Hs.eg.db)
      # gene.df1 <- bitr(gene.df$ENTREZID, fromType = "ENTREZID",
      #                  toType = c("SYMBOL"),
      #                  OrgDb = org.Hs.eg.db)
      # gene_list <- gene.df1$SYMBOL
      # write.table(gene_list, "1.ID_transformed_genelist.txt",row.names = F, quote = F)

      # 将genelist和表达矩阵的基因名称格式统一
      gene_list <- gsub("-", ".", gene_list)
      gene_list <- gsub("_", ".", gene_list)
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("-", ".", colnames(inputSet)[4:ncol(inputSet)])
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("_", ".", colnames(inputSet)[4:ncol(inputSet)])

      print("Gets the intersection of genelist and expression profile")
      # 获取genelist和表达谱的交集
      comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)], gene_list)
      # write.table(comsa1,"2.intersection_genelist_exprSet_gene.txt", row.names = F, quote = F)

      print("Processing the  input representation matrix")
      # 对输入的表达矩阵进行处理
      inputSet <- inputSet[, c("ID", "OS.time", "OS", comsa1)]

      inputSet[, c(1:2)] <- apply(inputSet[, c(1:2)], 2, as.factor)
      inputSet[, c(2:ncol(inputSet))] <- apply(inputSet[, c(2:ncol(inputSet))], 2, as.numeric)
      inputSet <- as.data.frame(inputSet)
      # rownames(inputSet) <- inputSet$ID

      print("Data preprocessing completed")
      # 自定义显示进程函数
      display.progress <- function(index, totalN, breakN = 20) {
        if (index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }


      ############### 单变量cox#######
      print("Stating the univariable cox regression")


      # KM估计筛选基因
      kmoutput <- NULL

      for (i in 1:ncol(inputSet[, 4:ncol(inputSet)])) {
        display.progress(index = i, totalN = ncol(inputSet), breakN = 20)
        g <- colnames(inputSet[, 4:ncol(inputSet)])[i]
        tmp <- inputSet[, c("OS.time", "OS", g)]
        tmp$group <- ifelse(tmp[, 3] > median(tmp[, 3]), "High", "Low")
        fitd <- survdiff(Surv(OS.time, OS) ~ group, data = tmp, na.action = na.exclude)
        p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
        kmoutput <- rbind(kmoutput, data.frame(
          gene = g,
          pvalue = p.val,
          stringsAsFactors = F
        ))
      }

      # write.csv(kmoutput,"KM_results.csv")

      print("Finished the KM selection")


      # 进行变量筛选
      selgene <- kmoutput[which(kmoutput$pvalue < KM_pcutoff), "gene"]
      return(selgene)
    }


    SigUnicox <- function(gene_list,
                          inputSet,
                          unicox_pcutoff # 单因素回归的筛选阈值
    ) {
      print("Starting the data preprocess")
      ############### 数据预处理#######

      print("Rejecting a null value")

      library(survival)
      # 将空值的基因变成0
      # table(is.na(inputSet))
      inputSet[is.na(inputSet)] <- 0
      # table(is.na(inputSet))
      inputSet <- inputSet %>% as.data.frame()
      inputSet$OS.time <- as.numeric(inputSet$OS.time) # 时间数值化
      inputSet <- inputSet[inputSet$OS.time > 0, ] # 剔除时间为0

      # print("Correcting gene set")
      # #将基因集矫正 策略为转化从ALIAS为ENTREZID， 然后再转化为SYMBOL,避免别名
      # gene.df <- bitr(gene_list, fromType = "ALIAS",
      #                 toType = c("ENTREZID"),
      #                 OrgDb = org.Hs.eg.db)
      # gene.df1 <- bitr(gene.df$ENTREZID, fromType = "ENTREZID",
      #                  toType = c("SYMBOL"),
      #                  OrgDb = org.Hs.eg.db)
      # gene_list <- gene.df1$SYMBOL
      # write.table(gene_list, "1.ID_transformed_genelist.txt",row.names = F, quote = F)

      # 将genelist和表达矩阵的基因名称格式统一
      gene_list <- gsub("-", ".", gene_list)
      gene_list <- gsub("_", ".", gene_list)
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("-", ".", colnames(inputSet)[4:ncol(inputSet)])
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("_", ".", colnames(inputSet)[4:ncol(inputSet)])

      print("Gets the intersection of genelist and expression profile")
      # 获取genelist和表达谱的交集
      comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)], gene_list)
      # write.table(comsa1,"2.intersection_genelist_exprSet_gene.txt", row.names = F, quote = F)

      print("Processing the  input representation matrix")
      # 对输入的表达矩阵进行处理
      inputSet <- inputSet[, c("ID", "OS.time", "OS", comsa1)]

      inputSet[, c(1:2)] <- apply(inputSet[, c(1:2)], 2, as.factor)
      inputSet[, c(2:ncol(inputSet))] <- apply(inputSet[, c(2:ncol(inputSet))], 2, as.numeric)
      inputSet <- as.data.frame(inputSet)
      # rownames(inputSet) <- inputSet$ID

      print("Data preprocessing completed")
      # 自定义显示进程函数
      display.progress <- function(index, totalN, breakN = 20) {
        if (index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }


      ############### 单变量cox#######
      print("Stating the univariable cox regression")

      unicox <- data.frame()
      for (i in 1:ncol(inputSet[, 4:ncol(inputSet)])) {
        display.progress(index = i, totalN = ncol(inputSet[, 4:ncol(inputSet)]))
        gene <- colnames(inputSet[, 4:ncol(inputSet)])[i]
        tmp <- data.frame(
          expr = as.numeric(inputSet[, 4:ncol(inputSet)][, i]),
          futime = inputSet$OS.time,
          fustat = inputSet$OS,
          stringsAsFactors = F
        )
        cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
        coxSummary <- summary(cox)
        unicox <- rbind.data.frame(unicox,
          data.frame(
            gene = gene,
            HR = as.numeric(coxSummary$coefficients[, "exp(coef)"])[1],
            z = as.numeric(coxSummary$coefficients[, "z"])[1],
            pvalue = as.numeric(coxSummary$coefficients[, "Pr(>|z|)"])[1],
            lower = as.numeric(coxSummary$conf.int[, 3][1]),
            upper = as.numeric(coxSummary$conf.int[, 4][1]),
            stringsAsFactors = F
          ),
          stringsAsFactors = F
        )
      }

      # write.csv(unicox,"3.unicox_results.csv")

      print("Finished the univariable cox regression")


      # 进行变量筛选
      selgene <- unicox[which(unicox$pvalue < unicox_pcutoff), "gene"]
      return(selgene)
      # write.table(selgene,paste("4.unicox_selected_cutoff_",unicox_pcutoff,"_genes.txt",sep = ""),row.names = F, quote = F)
    }
  }

  InputMatrix_pre <- InputMatrix
  colnames(InputMatrix_pre) <- gsub("-", ".", colnames(InputMatrix_pre))
  genelist.1 <- SigUnicox(gene_list = candidate_genes, inputSet = InputMatrix_pre, unicox_pcutoff = 0.05)

  genelist.2 <- SigKMcox(gene_list = genelist.1, inputSet = InputMatrix_pre, KM_pcutoff = 0.05)

  candidate_genes <- genelist.2

  print("----- finish the preprocess of the unicox and km analysis-----")

  ##### setting the pamameters ######

  rf_nodesize <- nodesize
  seed <- seed
  iter.times <- 1000
  # Checking data feasibility
  message("--- check data feasibility ---")

  # Replace '-' in column names with '.'

  candidate_genes <- gsub("-", ".", candidate_genes)
  colnames(InputMatrix) <- gsub("-", ".", colnames(InputMatrix))


  # Matching candidate genes to genes in each cohort
  common_feature <- c("ID", "OS.time", "OS", candidate_genes)
  common_feature <- intersect(common_feature, colnames(InputMatrix))

  message(paste0("---the number of the raw candidate genes is ", length(candidate_genes), " ---"))
  message(paste0("---the number of the common feature is ", length(common_feature) - 3, " ---"))

  ######### the main of the function ##########

  if (!is.na(rf_nodesize) &
    !is.na(seed) &
    mode %in% c("all", "single", "all_without_SVM") &
    identical(c("ID", "OS.time", "OS"), colnames(InputMatrix)[1:3]) &
    length(candidate_genes) > 0 &
    identical(c("ID", "OS.time", "OS"), common_feature[1:3]) &
    length(common_feature) > 3) {
    message("--- Data preprocessing ---")
    # Data preprocessing

    # Matching candidate genes to genes in each cohort

    InputMatrix <- InputMatrix[, common_feature]
    InputMatrix[, -c(1:3)] <- apply(InputMatrix[, -c(1:3)], 2, as.numeric)
    InputMatrix[, c(1:2)] <- apply(InputMatrix[, c(1:2)], 2, as.factor)
    InputMatrix[, c(2:3)] <- apply(InputMatrix[, c(2:3)], 2, as.numeric)
    InputMatrix <- InputMatrix[!is.na(InputMatrix$OS.time) & !is.na(InputMatrix$OS), ]
    InputMatrix <- InputMatrix[InputMatrix$OS.time > 0, ]

    InputMatrix[, -c(1:3)] <- apply(InputMatrix[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })

    est_dd <- as.data.frame(InputMatrix)[, common_feature[-1]]
    pre_var <- common_feature[-c(1:3)]
    selected.feature <- data.frame()


    if (mode == "all") {
      ### 1. Repeated Lasso  #############
      message("--- 1.Repeated lasso ---")
      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("1000 time lasso penalty")
      # 1000 time lasso penalty
      lasso_fea_list <- list()
      list.of.seed <- 1:iter.times

      print("This step will probably take several hours")

      lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
        set.seed(list.of.seed[x])
        cvfit <- cv.glmnet(
          x = x1,
          y = x2,
          nfolds = 10, # 10-fold交叉验证选取最优lambda
          alpha = 1, # alpha = 1 意味着 lasso
          family = "cox", # 依赖cox模型
          maxit = 1000
        )

        # optimal lambda
        fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
        if (is.element("(Intercept)", fea)) {
          lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
        } else {
          lasso_fea <- sort(fea)
        }
        return(lasso_fea)
      })

      # 输出每次运行的基因集合
      lasso_res <- NULL
      for (i in 1:iter.times) {
        lasso_res <- rbind.data.frame(lasso_res,
          data.frame(
            iteration = i,
            n.gene = length(lasso_fea_list[[i]]),
            genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
            stringsAsFactors = F
          ),
          stringsAsFactors = F
        )
      }




      genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # 根据基因出现的频次排序
      freq.cutoff <- iter.times * 0.05
      genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因. 95%


      result <- data.frame(
        method = c(rep("Lasso", length(genes))),
        selected.fea = genes
      )

      selected.feature <- rbind(selected.feature, result)



      ##### 2.Enet ###########
      message("--- 2.Enet  ---")


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("This step will probably take several hours")

      for (alpha in seq(0.1, 0.9, 0.1)) {
        message(paste0("--- 2.Enet ", alpha, "---"))

        set.seed(seed)
        # 1000 time
        fea_list <- list()
        list.of.seed <- 1:iter.times
        fea_list <- pblapply(list.of.seed, function(x) { # 大概运行2天
          set.seed(list.of.seed[x])
          cvfit <- cv.glmnet(
            x = x1,
            y = x2,
            nfolds = 10, # 10-fold交叉验证选取最优lambda
            alpha = alpha,
            family = "cox", # 依赖cox模型
            maxit = 1000
          )
          # 取出最优lambda
          fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
          if (is.element("(Intercept)", fea)) {
            lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
          } else {
            lasso_fea <- sort(fea)
          }
          return(lasso_fea)
        })


        genes <- sort(table(unlist(fea_list)), decreasing = T) # 根据基因出现的频次排序
        freq.cutoff <- iter.times * 0.05
        genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因

        result <- data.frame(
          method = c(rep(paste0("Enet", "[α=", alpha, "]"), length(genes))),
          selected.fea = genes
        )

        selected.feature <- rbind(selected.feature, result)
      }

      ##### 3.SVM-REF ###########
      message("--- 3.SVM-REF  ---")
      print("This step will probably take several hours")


      input <- est_dd[, -1]



      # 10CV (k-fold crossValidation）
      svmRFE(input, k = 10, halve.above = 100) # 分割数据，分配随机数
      nfold <- 10
      nrows <- nrow(input)
      folds <- rep(1:nfold, len = nrows)[sample(nrows)]
      folds <- lapply(1:nfold, function(x) which(folds == x))
      results <- lapply(folds, svmRFE.wrap, input, k = 10, halve.above = 100) # 特征选择
      top.features <- WriteFeatures(results, input, save = F) # 查看主要变量
      n.features <- nrow(top.features)
      if (n.features > 300) {
        n.svm <- 300
      } else {
        n.svm <- n.features
      }

      featsweep <- base::lapply(1:n.svm, FeatSweep.wrap, results, input)

      no.info <- min(prop.table(table(input[, 1])))
      errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
      fea <- top.features[1:which.min(errors), "FeatureName"]

      result <- data.frame(
        method = c(rep("SVM-REF", length(fea))),
        selected.fea = fea
      )

      selected.feature <- rbind(selected.feature, result)

      ##### 4.Boruta ###########
      set.seed(seed)

      message("--- 4.Boruta  ---")

      boruta <- Boruta(
        x = as.matrix(est_dd[, -c(1, 2)]), y = as.factor(est_dd[, c(2)]), pValue = 0.01, mcAdj = T,
        maxRuns = 1000
      )
      boruta.imp <- function(x) {
        imp <- reshape2::melt(x$ImpHistory, na.rm = T)[, -1]
        colnames(imp) <- c("Variable", "Importance")
        imp <- imp[is.finite(imp$Importance), ]

        variableGrp <- data.frame(
          Variable = names(x$finalDecision),
          finalDecision = x$finalDecision
        )

        showGrp <- data.frame(
          Variable = c("shadowMax", "shadowMean", "shadowMin"),
          finalDecision = c("shadowMax", "shadowMean", "shadowMin")
        )

        variableGrp <- rbind(variableGrp, showGrp)

        boruta.variable.imp <- merge(imp, variableGrp, all.x = T)

        sortedVariable <- boruta.variable.imp %>%
          group_by(Variable) %>%
          summarise(median = median(Importance)) %>%
          arrange(median)
        sortedVariable <- as.vector(sortedVariable$Variable)


        boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels = sortedVariable)

        invisible(boruta.variable.imp)
      }

      boruta.variable.imp <- boruta.imp(boruta)
      # head(boruta.variable.imp)
      boruta.finalVars <- data.frame(Item = getSelectedAttributes(boruta, withTentative = T), Type = "Boruta")

      result <- data.frame(
        method = c(rep("Boruta", length(boruta.finalVars$Item))),
        selected.fea = boruta.finalVars$Item
      )

      selected.feature <- rbind(selected.feature, result)


      ##### 5.Xgboost ###########

      message("--- 5.Xgboost  ---")
      train <- apply(est_dd[, -c(1)], 2, as.numeric) %>% as.data.frame()
      train_matrix <- sparse.model.matrix(OS ~ . - 1, data = train)
      train_label <- as.numeric(train$OS)
      train_fin <- list(data = train_matrix, label = train_label)
      dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)
      # 模型训练
      xgb <- xgboost(
        data = dtrain, max_depth = 6, eta = 0.5,
        objective = "binary:logistic", nround = 25
      )
      # 重要重要性排序
      importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)
      head(importance)
      importance$rel.imp <- importance$Gain / max(importance$Gain)

      # cutoff0.05
      xgboost.variable.imp <- importance[!importance$rel.imp < 0.05, ]
      xgboost.finalVars <- xgboost.variable.imp$Feature

      result <- data.frame(
        method = c(rep("Xgboost", length(xgboost.finalVars))),
        selected.fea = xgboost.finalVars
      )

      selected.feature <- rbind(selected.feature, result)


      ##### 6.RSF ###########
      message("--- 6.RSF  ---")

      fit <- rfsrc(Surv(OS.time, OS) ~ .,
        data = est_dd,
        ntree = 1000, nodesize = rf_nodesize, # 该值建议多调整
        splitrule = "logrank",
        importance = T,
        proximity = T,
        forest = T,
        seed = seed
      )
      rid <- var.select(object = fit, conservative = "high")
      rid <- rid$topvars

      result <- data.frame(
        method = c(rep("RSF", length(rid))),
        selected.fea = rid
      )

      selected.feature <- rbind(selected.feature, result)


      ##### 7.CoxBoost ###########
      message("--- 7.CoxBoost  ---")

      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
        trace = TRUE, start.penalty = 500, parallel = T
      )

      cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
        maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty
      )
      fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
        stepno = cv.res$optimal.step, penalty = pen$penalty
      )
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]
      result <- data.frame(
        method = c(rep("CoxBoost", length(rid))),
        selected.fea = rid
      )

      selected.feature <- rbind(selected.feature, result)

      ##### 8.StepCox ###########
      message("--- 8.StepCox ---")

      for (direction in c("both", "backward", "forward")) {
        fit <- stats::step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
        rid <- names(coef(fit)) # 这里不用卡P值，迭代的结果就是可以纳入的基因

        result <- data.frame(
          method = c(rep(paste0("StepCox", "+", direction), length(rid))),
          selected.fea = rid
        )

        selected.feature <- rbind(selected.feature, result)
      }

      return(selected.feature)
    } else if (mode == "single") {
      if (single_ml == "Lasso") {
        ### 1. Repeated Lasso  #############
        message("--- 1.Repeated lasso ---")
        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
        print("1000 time lasso penalty")
        # 1000 time lasso penalty
        lasso_fea_list <- list()
        list.of.seed <- 1:iter.times

        print("This step will probably take several hours")

        lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
          set.seed(list.of.seed[x])
          cvfit <- cv.glmnet(
            x = x1,
            y = x2,
            nfolds = 10, # 10-fold交叉验证选取最优lambda
            alpha = 1, # alpha = 1 意味着 lasso
            family = "cox", # 依赖cox模型
            maxit = 1000
          )

          # optimal lambda
          fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
          if (is.element("(Intercept)", fea)) {
            lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
          } else {
            lasso_fea <- sort(fea)
          }
          return(lasso_fea)
        })

        # 输出每次运行的基因集合
        lasso_res <- NULL
        for (i in 1:iter.times) {
          lasso_res <- rbind.data.frame(lasso_res,
            data.frame(
              iteration = i,
              n.gene = length(lasso_fea_list[[i]]),
              genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
              stringsAsFactors = F
            ),
            stringsAsFactors = F
          )
        }




        genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # 根据基因出现的频次排序
        freq.cutoff <- iter.times * 0.05
        genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因. 95%


        result <- data.frame(
          method = c(rep("Lasso", length(genes))),
          selected.fea = genes
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
      } else if (single_ml == "Enet") {
        ##### 2.Enet ###########
        message("--- 2.Enet  ---")


        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
        print("This step will probably take several hours")

        for (alpha in seq(0.1, 0.9, 0.1)) {
          message(paste0("--- 2.Enet ", alpha, "---"))

          set.seed(seed)
          # 1000 time
          fea_list <- list()
          list.of.seed <- 1:iter.times
          fea_list <- pblapply(list.of.seed, function(x) { # 大概运行2天
            set.seed(list.of.seed[x])
            cvfit <- cv.glmnet(
              x = x1,
              y = x2,
              nfolds = 10, # 10-fold交叉验证选取最优lambda
              alpha = alpha,
              family = "cox", # 依赖cox模型
              maxit = 1000
            )
            # 取出最优lambda
            fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
            if (is.element("(Intercept)", fea)) {
              lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
            } else {
              lasso_fea <- sort(fea)
            }
            return(lasso_fea)
          })


          genes <- sort(table(unlist(fea_list)), decreasing = T) # 根据基因出现的频次排序
          freq.cutoff <- iter.times * 0.05
          genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因

          result <- data.frame(
            method = c(rep(paste0("Enet", "[α=", alpha, "]"), length(genes))),
            selected.fea = genes
          )

          selected.feature <- rbind(selected.feature, result)
        }

        return(selected.feature)
      } else if (single_ml == "SVM-REF") {
        ##### 3.SVM-REF ###########
        message("--- 3.SVM-REF  ---")
        print("This step will probably take several hours")


        input <- est_dd[, -1]
        # 10CV (k-fold crossValidation）
        svmRFE(input, k = 10, halve.above = 100) # 分割数据，分配随机数
        nfold <- 10
        nrows <- nrow(input)
        folds <- rep(1:nfold, len = nrows)[sample(nrows)]
        folds <- lapply(1:nfold, function(x) which(folds == x))
        results <- lapply(folds, svmRFE.wrap, input, k = 10, halve.above = 100) # 特征选择
        top.features <- WriteFeatures(results, input, save = F) # 查看主要变量
        n.features <- nrow(top.features)
        if (n.features > 300) {
          n.svm <- 300
        } else {
          n.svm <- n.features
        }

        featsweep <- base::lapply(1:n.svm, FeatSweep.wrap, results, input)

        no.info <- min(prop.table(table(input[, 1])))
        errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
        fea <- top.features[1:which.min(errors), "FeatureName"]

        result <- data.frame(
          method = c(rep("SVM-REF", length(fea))),
          selected.fea = fea
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
      } else if (single_ml == "Boruta") {
        ##### 4.Boruta ###########
        set.seed(seed)

        message("--- 4.Boruta  ---")

        boruta <- Boruta(
          x = as.matrix(est_dd[, -c(1, 2)]), y = as.factor(est_dd[, c(2)]), pValue = 0.01, mcAdj = T,
          maxRuns = 1000
        )
        boruta.imp <- function(x) {
          imp <- reshape2::melt(x$ImpHistory, na.rm = T)[, -1]
          colnames(imp) <- c("Variable", "Importance")
          imp <- imp[is.finite(imp$Importance), ]

          variableGrp <- data.frame(
            Variable = names(x$finalDecision),
            finalDecision = x$finalDecision
          )

          showGrp <- data.frame(
            Variable = c("shadowMax", "shadowMean", "shadowMin"),
            finalDecision = c("shadowMax", "shadowMean", "shadowMin")
          )

          variableGrp <- rbind(variableGrp, showGrp)

          boruta.variable.imp <- merge(imp, variableGrp, all.x = T)

          sortedVariable <- boruta.variable.imp %>%
            group_by(Variable) %>%
            summarise(median = median(Importance)) %>%
            arrange(median)
          sortedVariable <- as.vector(sortedVariable$Variable)


          boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels = sortedVariable)

          invisible(boruta.variable.imp)
        }

        boruta.variable.imp <- boruta.imp(boruta)
        # head(boruta.variable.imp)
        boruta.finalVars <- data.frame(Item = getSelectedAttributes(boruta, withTentative = T), Type = "Boruta")

        result <- data.frame(
          method = c(rep("Boruta", length(boruta.finalVars$Item))),
          selected.fea = boruta.finalVars$Item
        )

        selected.feature <- rbind(selected.feature, result)
        return(selected.feature)
      } else if (single_ml == "Xgboost") {
        ##### 5.Xgboost ###########

        message("--- 5.Xgboost  ---")
        train <- apply(est_dd[, -c(1)], 2, as.numeric) %>% as.data.frame()
        train_matrix <- sparse.model.matrix(OS ~ . - 1, data = train)
        train_label <- as.numeric(train$OS)
        train_fin <- list(data = train_matrix, label = train_label)
        dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)
        # 模型训练
        xgb <- xgboost(
          data = dtrain, max_depth = 6, eta = 0.5,
          objective = "binary:logistic", nround = 25
        )
        # 重要重要性排序
        importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)
        head(importance)
        importance$rel.imp <- importance$Gain / max(importance$Gain)

        # cutoff0.05
        xgboost.variable.imp <- importance[!importance$rel.imp < 0.05, ]
        xgboost.finalVars <- xgboost.variable.imp$Feature

        result <- data.frame(
          method = c(rep("Xgboost", length(xgboost.finalVars))),
          selected.fea = xgboost.finalVars
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
      } else if (single_ml == "RSF") {
        ##### 6.RSF ###########
        message("--- 6.RSF  ---")

        fit <- rfsrc(Surv(OS.time, OS) ~ .,
          data = est_dd,
          ntree = 1000, nodesize = rf_nodesize,
          splitrule = "logrank",
          importance = T,
          proximity = T,
          forest = T,
          seed = seed
        )
        rid <- var.select(object = fit, conservative = "high")
        rid <- rid$topvars

        result <- data.frame(
          method = c(rep("RSF", length(rid))),
          selected.fea = rid
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
      } else if (single_ml == "CoxBoost") {
        ##### 7.CoxBoost ###########
        message("--- 7.CoxBoost  ---")

        set.seed(seed)
        pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
          trace = TRUE, start.penalty = 500, parallel = T
        )

        cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
          maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty
        )
        fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
          stepno = cv.res$optimal.step, penalty = pen$penalty
        )
        rid <- as.data.frame(coef(fit))
        rid$id <- rownames(rid)
        rid <- rid[which(rid$`coef(fit)` != 0), "id"]
        result <- data.frame(
          method = c(rep("CoxBoost", length(rid))),
          selected.fea = rid
        )

        selected.feature <- rbind(selected.feature, result)
        return(selected.feature)
      } else if (single_ml == "StepCox") {
        ##### 8.StepCox ###########
        message("--- 8.StepCox ---")

        for (direction in c("both", "backward", "forward")) {
          fit <- stats::step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
          rid <- names(coef(fit)) #

          result <- data.frame(
            method = c(rep(paste0("StepCox", "+", direction), length(rid))),
            selected.fea = rid
          )

          selected.feature <- rbind(selected.feature, result)
        }

        return(selected.feature)
      } else {
        warning("The parameter of the single ML is out of the bound")
      }
    } else if (mode == "all_without_SVM") {
      ### 1. Repeated Lasso  #############
      message("--- 1.Repeated lasso ---")
      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("1000 time lasso penalty")
      # 1000 time lasso penalty
      lasso_fea_list <- list()
      list.of.seed <- 1:iter.times

      print("This step will probably take several hours")

      lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
        set.seed(list.of.seed[x])
        cvfit <- cv.glmnet(
          x = x1,
          y = x2,
          nfolds = 10, # 10-fold交叉验证选取最优lambda
          alpha = 1, # alpha = 1 意味着 lasso
          family = "cox", # 依赖cox模型
          maxit = 1000
        )

        # optimal lambda
        fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
        if (is.element("(Intercept)", fea)) {
          lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
        } else {
          lasso_fea <- sort(fea)
        }
        return(lasso_fea)
      })

      # 输出每次运行的基因集合
      lasso_res <- NULL
      for (i in 1:iter.times) {
        lasso_res <- rbind.data.frame(lasso_res,
          data.frame(
            iteration = i,
            n.gene = length(lasso_fea_list[[i]]),
            genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
            stringsAsFactors = F
          ),
          stringsAsFactors = F
        )
      }




      genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # 根据基因出现的频次排序
      freq.cutoff <- iter.times * 0.05
      genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因. 95%


      result <- data.frame(
        method = c(rep("Lasso", length(genes))),
        selected.fea = genes
      )

      selected.feature <- rbind(selected.feature, result)



      ##### 2.Enet ###########
      message("--- 2.Enet  ---")


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("This step will probably take several hours")

      for (alpha in seq(0.1, 0.9, 0.1)) {
        message(paste0("--- 2.Enet ", alpha, "---"))

        set.seed(seed)
        # 1000 time
        fea_list <- list()
        list.of.seed <- 1:iter.times
        fea_list <- pblapply(list.of.seed, function(x) { # 大概运行2天
          set.seed(list.of.seed[x])
          cvfit <- cv.glmnet(
            x = x1,
            y = x2,
            nfolds = 10, # 10-fold交叉验证选取最优lambda
            alpha = alpha,
            family = "cox", # 依赖cox模型
            maxit = 1000
          )
          # 取出最优lambda
          fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
          if (is.element("(Intercept)", fea)) {
            lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
          } else {
            lasso_fea <- sort(fea)
          }
          return(lasso_fea)
        })


        genes <- sort(table(unlist(fea_list)), decreasing = T) # 根据基因出现的频次排序
        freq.cutoff <- iter.times * 0.05
        genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因

        result <- data.frame(
          method = c(rep(paste0("Enet", "[α=", alpha, "]"), length(genes))),
          selected.fea = genes
        )

        selected.feature <- rbind(selected.feature, result)
      }

      ##### 3.SVM-REF ###########
      # message('--- 3.SVM-REF  ---')
      # print('This step will probably take several hours')
      #
      #
      # input <- est_dd[,-1]
      #
      #
      #
      # # 10CV (k-fold crossValidation）
      # svmRFE(input, k = 10, halve.above = 100) #分割数据，分配随机数
      # nfold = 10
      # nrows = nrow(input)
      # folds = rep(1:nfold, len=nrows)[sample(nrows)]
      # folds = lapply(1:nfold, function(x) which(folds == x))
      # results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) #特征选择
      # top.features = WriteFeatures(results, input, save=F) #查看主要变量
      # n.features = nrow(top.features)
      # if(n.features > 300){
      #   n.svm = 300
      # }else {
      #   n.svm =n.features
      # }
      #
      # featsweep = base::lapply(1:n.svm, FeatSweep.wrap, results, input)
      #
      # no.info = min(prop.table(table(input[,1])))
      # errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
      # fea <- top.features[1:which.min(errors), "FeatureName"]
      #
      # result <-  data.frame(method = c(rep("SVM-REF", length(fea))),
      #                       selected.fea= fea)
      #
      # selected.feature <- rbind(selected.feature,result)

      ##### 3.Boruta ###########
      set.seed(seed)

      message("--- 3.Boruta  ---")

      boruta <- Boruta(
        x = as.matrix(est_dd[, -c(1, 2)]), y = as.factor(est_dd[, c(2)]), pValue = 0.01, mcAdj = T,
        maxRuns = 1000
      )
      boruta.imp <- function(x) {
        imp <- reshape2::melt(x$ImpHistory, na.rm = T)[, -1]
        colnames(imp) <- c("Variable", "Importance")
        imp <- imp[is.finite(imp$Importance), ]

        variableGrp <- data.frame(
          Variable = names(x$finalDecision),
          finalDecision = x$finalDecision
        )

        showGrp <- data.frame(
          Variable = c("shadowMax", "shadowMean", "shadowMin"),
          finalDecision = c("shadowMax", "shadowMean", "shadowMin")
        )

        variableGrp <- rbind(variableGrp, showGrp)

        boruta.variable.imp <- merge(imp, variableGrp, all.x = T)

        sortedVariable <- boruta.variable.imp %>%
          group_by(Variable) %>%
          summarise(median = median(Importance)) %>%
          arrange(median)
        sortedVariable <- as.vector(sortedVariable$Variable)


        boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels = sortedVariable)

        invisible(boruta.variable.imp)
      }

      boruta.variable.imp <- boruta.imp(boruta)
      # head(boruta.variable.imp)
      boruta.finalVars <- data.frame(Item = getSelectedAttributes(boruta, withTentative = T), Type = "Boruta")

      result <- data.frame(
        method = c(rep("Boruta", length(boruta.finalVars$Item))),
        selected.fea = boruta.finalVars$Item
      )

      selected.feature <- rbind(selected.feature, result)


      ##### 4.Xgboost ###########

      message("--- 4.Xgboost  ---")
      train <- apply(est_dd[, -c(1)], 2, as.numeric) %>% as.data.frame()
      train_matrix <- sparse.model.matrix(OS ~ . - 1, data = train)
      train_label <- as.numeric(train$OS)
      train_fin <- list(data = train_matrix, label = train_label)
      dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)
      # 模型训练
      xgb <- xgboost(
        data = dtrain, max_depth = 6, eta = 0.5,
        objective = "binary:logistic", nround = 25
      )
      # 重要重要性排序
      importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)
      head(importance)
      importance$rel.imp <- importance$Gain / max(importance$Gain)

      # cutoff0.05
      xgboost.variable.imp <- importance[!importance$rel.imp < 0.05, ]
      xgboost.finalVars <- xgboost.variable.imp$Feature

      result <- data.frame(
        method = c(rep("Xgboost", length(xgboost.finalVars))),
        selected.fea = xgboost.finalVars
      )

      selected.feature <- rbind(selected.feature, result)


      ##### 5.RSF ###########
      message("--- 5.RSF  ---")

      fit <- rfsrc(Surv(OS.time, OS) ~ .,
        data = est_dd,
        ntree = 1000, nodesize = rf_nodesize, # 该值建议多调整
        splitrule = "logrank",
        importance = T,
        proximity = T,
        forest = T,
        seed = seed
      )
      rid <- var.select(object = fit, conservative = "high")
      rid <- rid$topvars

      result <- data.frame(
        method = c(rep("RSF", length(rid))),
        selected.fea = rid
      )

      selected.feature <- rbind(selected.feature, result)


      ##### 6.CoxBoost ###########
      message("--- 6.CoxBoost  ---")

      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
        trace = TRUE, start.penalty = 500, parallel = T
      )

      cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
        maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty
      )
      fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
        stepno = cv.res$optimal.step, penalty = pen$penalty
      )
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]
      result <- data.frame(
        method = c(rep("CoxBoost", length(rid))),
        selected.fea = rid
      )

      selected.feature <- rbind(selected.feature, result)

      ##### 7.StepCox ###########
      message("--- 7.StepCox ---")

      for (direction in c("both", "backward", "forward")) {
        fit <- stats::step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
        rid <- names(coef(fit)) # 这里不用卡P值，迭代的结果就是可以纳入的基因

        result <- data.frame(
          method = c(rep(paste0("StepCox", "+", direction), length(rid))),
          selected.fea = rid
        )

        selected.feature <- rbind(selected.feature, result)
      }

      return(selected.feature)
    }
  } else {
    print("Please set the full parameters")
  }
}
