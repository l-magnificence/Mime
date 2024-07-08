#' Developing the optimal predictive model for the dichotomous variables with machine learning algorithms
#' 
#' A function can be used to develop the predictive model for dichotomous variables with seven machine learning algorithms.
#' 
#' @param train_data The training data with the 'ID' and 'Var' as the first two columns. Starting in the third column are the variables used to construct the model. 'Var' is the target predictor variable for constructing the model. 'Var' contains only Y or N.
#' @param list_train_vali_Data A list containing the training data and the other validation data. All the validation data have the same data form as the training data.
#' @param candidate_genes The candidate variables used for constructing the predictive model.
#' @param methods There are seven algorithms for developing the predictive model including 'nb', 'svmRadialWeights', 'rf', 'kknn', 'adaboost', 'LogitBoost', 'cancerclass'. 'nb':Naive Bayes algorithm. 'svmRadialWeights': Support Vector Machine (SVM). 'rf': Random Forest. 'kknn': K-nearest Neighbors.'adaboost': AdaBoost Classification Trees. 'LogitBoost':Boosted Logistic Regressions. 'cancerclass': Cancerclass. 
#' @param seed The seed you can set as any positive number, for example, 5201314.
#' @param cores_for_parallel The cores you can choose for parallel operation. The default is 12.The bigger the better if the configuration allows it.
#'
#' @return A list containing the predictive model, the AUC, the ROC, and the candidate variables, all of which are developed by each single algorithm.
#' @export
#'
#' @examples
ML.Dev.Pred.Category.Sig <- function(train_data, # cohort data used for training, the colnames of which inlcuding ID, Var, and the other candidate genes。
                                     # Var 是用于构建预测模型的目标变量，Y/N，
                                     list_train_vali_Data, # cohort data used for training, the colnames of which inlcuding ID, Var, and the other candidate genes。
                                     # Var 是用于构建预测模型的目标变量，Y/N
                                     candidate_genes = NULL,
                                     methods = NULL, # c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
                                     seed = 5201314, # 5201314
                                     cores_for_parallel = 12 #
) {
  ### === loading packages ===###
  message("---loading the packages ---")


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


  ###### loading the function #######



  model.Dev <- function(training, method, sig) {
    training <- training[, colnames(training) %in% c("Var", sig)]
    # 7 models adpoted in this study as followings:
    #' nb': navie bayes
    #' svmRadialWeights': Support Vector Machines with Class Weights
    #' rf': random forest
    #' kknn': k-Nearest Neighbors
    #' adaboost':AdaBoost Classification Trees
    #' LogitBoost':Boosted Logistic Regressions
    #' cancerclass': cancerclass

    # Grid search for parameter tuning
    Grid <- list(
      nb = expand.grid(fL = c(0, 0.5, 1, 1.5, 2.0), usekernel = TRUE, adjust = c(0.5, 0.75, 1, 1.25, 1.5)),
      svmRadialWeights = expand.grid(sigma = c(0.0005, 0.001, 0.005, 0.01, 0.05), C = c(1, 3, 5, 10, 20), Weight = c(0.1, 0.5, 1, 2, 3, 5, 10)),
      rf = expand.grid(mtry = c(2, 42, 83, 124, 165, 205, 246, 287, 328, 369)),
      kknn = expand.grid(kmax = c(5, 7, 9, 11, 13), distance = 2, kernel = "optimal"),
      adaboost = expand.grid(nIter = c(50, 100, 150, 200, 250), method = c("Adaboost.M1", "Real adaboost")),
      LogitBoost = expand.grid(nIter = c(11, 21, 31, 41, 51, 61, 71, 81, 91, 101))
    )
    TuneLength <- list(
      nb = nrow(Grid[["nb"]]),
      svmRadialWeights = nrow(Grid[["svmRadialWeights"]]),
      rf = nrow(Grid[["rf"]]),
      kknn = nrow(Grid[["kknn"]]),
      adaboost = nrow(Grid[["adaboost"]]),
      LogitBoost = nrow(Grid[["LogitBoost"]])
    )
    ls_model <- lapply(method, function(m) {
      if (m == "cancerclass") { # cancerclass is not avaliable in caret
        pData <- data.frame(class = training$Var, sample = rownames(training), row.names = rownames(training))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(training[, -1])
        Sig.Exp.train <- ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)
        predictor <- fit(Sig.Exp.train, method = "welch.test")
        model.tune <- predictor
      } else {
        f <- 5 # f folds resampling
        r <- 10 # r repeats
        n <- f * r

        # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
        seeds <- vector(mode = "list", length = n + 1)
        # the number of tuning parameter
        for (i in 1:n) seeds[[i]] <- sample.int(n = 1000, TuneLength[[m]])

        # for the last model
        seeds[[n + 1]] <- sample.int(1000, 1)


        ctrl <- trainControl(
          method = "repeatedcv",
          number = f, ## 5-folds cv
          summaryFunction = twoClassSummary, # Use AUC to pick the best model
          classProbs = TRUE,
          repeats = r, ## 10-repeats cv,
          seeds = seeds
        )



        model.tune <- train(Var ~ .,
          data = training,
          method = m,
          metric = "ROC",
          trControl = ctrl,
          tuneGrid = Grid[[m]]
        )
      }
      print(m)
      return(model.tune)
    })


    names(ls_model) <- method

    return(ls_model)
  }





  # CompareModel <- function(training, validation, method,sig){
  #
  #   training <- training[,colnames(training) %in% c('Var', sig)]
  #   validation  <- validation[,colnames(validation) %in% c('Var',sig)]
  #
  #
  #
  #
  #
  #
  #
  #
  #   auc <- lapply(ls_model,function(model.tune){
  #     if(class(model.tune) == 'predictor'){
  #       pData <- data.frame(class = validation$Var, sample = rownames(validation),row.names = rownames(validation))
  #       phenoData <- new("AnnotatedDataFrame",data=pData)
  #       Sig.Exp <- t(validation[,-1])
  #       Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
  #       prediction <- predict(model.tune, Sig.Exp.test,"N", ngenes=nrow(Sig.Exp), dist = "cor")
  #       roc <- roc(response  = prediction@prediction[,'class_membership'],
  #                  predictor = as.numeric(prediction@prediction[,'z'])
  #       )
  #       roc_result <- coords(roc, "best")
  #       auc <- data.frame(ROC=as.numeric(roc$auc), Sens = roc_result$sensitivity[1], Spec = roc_result$specificity[1])
  #     }else {
  #       prob <- predict(model.tune,validation[,-1],type = "prob")
  #       pre <- predict(model.tune,validation[,-1])
  #       test_set <- data.frame(obs = validation$Var, N = prob[,'N'], Y = prob[,'Y'], pred=pre)
  #       auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
  #     }
  #
  #     return(auc)
  #   }) %>% base::do.call(rbind,.)
  #
  #   rownames(auc) <- method
  #
  #   res <- list()
  #
  #   names(ls_model) = method
  #
  #   res[['model']] <- ls_model
  #   res[['auc']] <- auc
  #
  #
  #   return(res)
  #
  # }


  cal.model.auc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    library(dplyr)
    
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", "Y"))
    
    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    auc <- lapply(1:length(models), function(i) {
      if (models[i] == "cancerclass") {
        model.tune <- ls_model[[i]]
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
        model.tune <- ls_model[[i]]
        prob <- predict(model.tune, validation[, -1], type = "prob")
        pre <- predict(model.tune, validation[, -1])
        test_set <- data.frame(obs = validation$Var, N = prob[, "N"], Y = prob[, "Y"], pred = pre)
        auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
      }
      
      return(auc)
    }) %>% base::do.call(rbind, .)
    
    rownames(auc) <- names(ls_model)
    
    return(auc)
  }


  cal.model.roc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    library(dplyr)

    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", "Y"))

    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    roc <- lapply(1:length(models), function(i) {
      if (!models[i] == "cancerclass") {
        prob <- predict(ls_model[[models[i]]], validation[, -1], type = "prob") #
        pre <- predict(ls_model[[models[i]]], validation[, -1]) #
        test_set <- data.frame(obs = validation$Var, N = prob[, "N"], Y = prob[, "Y"], pred = pre)
        roc <- ROCit::rocit(
          score = test_set$N,
          class = test_set$obs,
          negref = "Y"
        )
      } else {
        pData <- data.frame(class = validation$Var, sample = rownames(validation), row.names = rownames(validation))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(validation[, -1])
        Sig.Exp.test <- ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)

        prediction <- predict(ls_model[[models[i]]], Sig.Exp.test, "N", ngenes = nrow(Sig.Exp), dist = "cor")
        roc <- roc(
          response = prediction@prediction[, "class_membership"],
          predictor = as.numeric(prediction@prediction[, "z"])
        )
      }
    })

    names(roc) <- models


    return(roc)
  }



  message("---loading the function---")


  common_feature <- c("ID", "Var", candidate_genes)
  list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
    colnames(x) = gsub('-','.',colnames(x))
    return(x)})
  colnames(train_data) <- gsub("-", ".", colnames(train_data))
  candidate_genes <- gsub("-", ".", candidate_genes)


  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  ##### parameters check #####


  if (
    #identical(colnames(train_data)[1:2], c("ID", "Var")) &
    #identical(common_feature[1:2], c("ID", "Var")) &
    #unique(train_data$Var) %in% c("Y", "N") &
    #length(candidate_genes) > 1 &
    all(is.element(methods, c("nb", "svmRadialWeights", "rf", "kknn", "adaboost", "LogitBoost", "cancerclass")))


  ) {
    ####### data preparation ######





    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, common_feature]
    })

    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, as.numeric)
      rownames(x) <- x$ID
      return(x)
    })

    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
      return(x)
    })


    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x <- x[!is.na(x$Var) & !is.na(x$Var), ]
      return(x)
    })

    # use the mean replace the NA
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = T)
        return(x)
      })


      return(x)
    })

    train_data <- train_data[, common_feature]
    train_data[, -c(1:2)] <- apply(train_data[, -c(1:2)], 2, as.numeric)
    train_data[, c(1:2)] <- apply(train_data[, c(1:2)], 2, as.factor)
    rownames(train_data) <- train_data$ID

    est_dd <- as.data.frame(train_data)[, common_feature[-1]]
    pre_var <- common_feature[-c(1:2)]



    print(paste0("There existing ", length(candidate_genes), " genes in candidate genes"))


    print("Intersetion of the candidate genes and the colnames of the provided data")

    print(paste0("There existing ", length(pre_var), " genes in candidate genes, colnames of training data, colnames of validation data"))




    # parallel processing


    cl <- makePSOCKcluster(cores_for_parallel)
    registerDoParallel(cl)

    res.model <- model.Dev(
      training = train_data,
      method = c("nb", "svmRadialWeights", "kknn", "rf", "adaboost", "LogitBoost", "cancerclass"),
      sig = pre_var
    )
    stopCluster(cl)



    ml.auc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.auc(res.by.model.Dev = res.model, cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.auc) <- names(list_train_vali_Data)

    ml.roc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.roc(res.by.model.Dev = res.model, cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.roc) <- names(list_train_vali_Data)

    res <- list()
    res[["model"]] <- res.model
    res[["auc"]] <- ml.auc
    res[["roc"]] <- ml.roc
    res[["sig.gene"]] <- pre_var

    return(res)
  } else {
    print("Please provide the correct parameters")
  }
}
