#' Constructing the optimal predictive model for the prognosis with machine learning algorithms
#' 
#' A function, which can be used for developing a robust predictive model with elegant performance to predict the prognosis of the patients.
#' 
#' @param train_data The training data with the ID, OS.time, and OS as the first three column names. Starting in the fourth column are the variables used to construct the model. The expression is scaled with log2(x+1). OS.time means the survival time(Day). OS means the survival status only containing 0 and 1 (1:Dead, 0:Alive). 
#' @param list_train_vali_Data A list containing the training data and some validation data. The validation data has the same format as the training data.
#' @param candidate_genes  The character vector containing the variables you just want to input for developing the predictive model. These variables should be included in the colnames of the training data.
#' @param unicox.filter.for.candi # T or F. wheather you use univariable cox regression to screen out the prognostic variables. The default is T. 
#' @param unicox_p_cutoff # The p value for the threshold of the Uni-Cox regression analysis, The default is 0.05
#' @param mode Here we provide three modes including 'all', 'single', and 'double'. 'all' means using all ten algorithms and the combinations. 'single' means using only one of the ten algorithms. 'double' means using the combination with two algorithms.
#' @param single_ml One of the ten algorithms including "RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", "survivalsvm", "Ridge", "Lasso".
#' @param alpha_for_Enet The parameter for the "Enet".  One of the values from 0.1 to 0.9. c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9). There are some conditions you could not set this parameter. [1] The mode is 'all'. [2] The mode is 'single' or 'double', but the 'Enet' is not included in the algorithms you choose. 
#' @param direction_for_stepcox The parameter for the StepCox. One  from "both", "backward", "forward". There are some conditions you could not set this parameter. [1] The mode is 'all'. [2] The mode is 'single' or 'double', but the 'StepCox' is not included in the algorithms you choose. 
#' @param double_ml1 The first algorithm of the algorithm combination. If you set the mode as the 'double', you should set this from the c('RSF', "StepCox","CoxBoost","Lasso"). Here we only provide four options.
#' @param double_ml2  The second algorithm of the algorithm combination. If you set the mode as the 'double', you should set this from the c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso"). Here we only provide ten options. The specific combinations of the algorithms are provided in the documentation.
#' @param nodesize The parameter for 'RSF'. The default is 5. You can try the positive number from 5 to 10.
#' @param seed The seed you can set as any positive integer, for example, 5201314.
#'
#' @return A list containing the developed predictive model, the risk score of the model in each data, the C index of the model in each data, and the variables used to construct the model.
#' @export
#'
#' @examples
ML.Dev.Prog.Sig = function(train_data, # cohort data used for training, the colnames of which inlcuding ID, OS.time, OS, and the other candidate genes
                           list_train_vali_Data, # a list of the validation data and the training data. The cohort data is the same to the cohort data used for training
                           # 要求队列的测序深度不能太低，将genelist和队列列名取交集之后尽量保证丢失的信息不超过20%
                           candidate_genes = NULL,
                          unicox.filter.for.candi = NULL, # 是否使用unicox 对基因进行筛选
                           unicox_p_cutoff = NULL, # 默认为0.05 unicox 筛选阈值
                           mode = NULL, # all, single, double
                           single_ml = NULL,# c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           alpha_for_Enet = NULL , # 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
                           direction_for_stepcox = NULL, #  c("both", "backward", "forward")
                           double_ml1 = NULL,#c('RSF', "StepCox","CoxBoost","Lasso")
                           double_ml2 = NULL,# c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           nodesize = NULL, # reference 5-10
                           seed = NULL, # 5201314
                           cores_for_parallel = NULL  #cores for gbm
                           ){
  
  if (is.null(alpha_for_Enet) == T){
    alpha_for_Enet<-0.1 ## default 0.35
  } else {
    alpha_for_Enet<-alpha_for_Enet
  }

  if (is.null(cores_for_parallel) == T){
    cores_for_parallel<-6 ## default 6
  } else {
    cores_for_parallel<-cores_for_parallel
  }

  
  if (is.null(direction_for_stepcox) == T){
    direction_for_stepcox<-'both' ## default 0.35
  } else {
    direction_for_stepcox<-direction_for_stepcox
  }
  
  
  
  #loading the packages
  if(T) {
    library(Matrix)
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



    Sys.setenv(LANGUAGE = "en") #显示英文报错信息
    options(stringsAsFactors = FALSE) #禁止chr转成factor
  }

    #loading the functions
  if(T){
    
    SigUnicox <- function(gene_list,
                          inputSet,
                          unicox_pcutoff #单因素回归的筛选阈值
    ){
      print("Starting the data preprocess")
      ###############数据预处理#######
      
      print("Rejecting a null value")
      
      library(survival)
      #将空值的基因变成0
      # table(is.na(inputSet))
      inputSet[is.na(inputSet)] = 0
      # table(is.na(inputSet))
      inputSet <- inputSet %>% as.data.frame()
      inputSet$OS.time <- as.numeric(inputSet$OS.time)#时间数值化
      inputSet <- inputSet[inputSet$OS.time>0,]#剔除时间为0
      
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
      
      #将genelist和表达矩阵的基因名称格式统一
      gene_list <- gsub("-",".",gene_list)
      gene_list <- gsub("_",".",gene_list)
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("-",".",colnames(inputSet)[4:ncol(inputSet)])
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("_",".",colnames(inputSet)[4:ncol(inputSet)])
      
      print("Gets the intersection of genelist and expression profile")
      #获取genelist和表达谱的交集
      comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)],gene_list)
      # write.table(comsa1,"2.intersection_genelist_exprSet_gene.txt", row.names = F, quote = F)
      
      print("Processing the  input representation matrix")
      #对输入的表达矩阵进行处理
      inputSet <- inputSet[,c("ID","OS.time","OS",comsa1)]
      
      inputSet[,c(1:2)] <- apply(inputSet[,c(1:2)],2,as.factor)
      inputSet[,c(2:ncol(inputSet))] <- apply(inputSet[,c(2:ncol(inputSet))],2,as.numeric)
      inputSet <- as.data.frame(inputSet)
      # rownames(inputSet) <- inputSet$ID
      
      print("Data preprocessing completed")
      #自定义显示进程函数
      display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  ==0  ) {
          cat(paste(round(index*100/totalN), "% ", sep=""))
        }
      }
      
      
      ###############单变量cox#######
      print("Stating the univariable cox regression")
      
      unicox <- data.frame()
      for(i in 1:ncol(inputSet[,4:ncol(inputSet)])){
        
        display.progress(index = i, totalN = ncol(inputSet[,4:ncol(inputSet)]))
        gene <- colnames(inputSet[,4:ncol(inputSet)])[i]
        tmp <- data.frame(expr = as.numeric(inputSet[,4:ncol(inputSet)][,i]),
                          futime = inputSet$OS.time,
                          fustat = inputSet$OS,
                          stringsAsFactors = F)
        cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
        coxSummary <- summary(cox)
        unicox <- rbind.data.frame(unicox,
                                   data.frame(gene = gene,
                                              HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                              z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                              pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                              lower = as.numeric(coxSummary$conf.int[,3][1]),
                                              upper = as.numeric(coxSummary$conf.int[,4][1]),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
      }
      
      # write.csv(unicox,"3.unicox_results.csv")
      
      print("Finished the univariable cox regression")
      
      
      #进行变量筛选
      selgene <- unicox[which(unicox$pvalue < unicox_pcutoff), "gene"]
      return(selgene)
      # write.table(selgene,paste("4.unicox_selected_cutoff_",unicox_pcutoff,"_genes.txt",sep = ""),row.names = F, quote = F)
    }
    
  }

  rf_nodesize <- nodesize
  seed <- seed

  # Checking data feasibility

  # Replace '-' in column names with '.'
  list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
    colnames(x) = gsub('-','.',colnames(x))
    return(x)})

  candidate_genes = gsub('-','.',candidate_genes)
  colnames(train_data) = gsub('-','.',colnames(train_data))

  # Matching candidate genes to genes in each cohort
  common_feature = c('ID', 'OS.time', 'OS',candidate_genes)

  for (i in names(list_train_vali_Data)) {
    common_feature = intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  message(paste0('---the number of the raw candidate genes is ', length(candidate_genes),' ---'))
  message(paste0('---the number of the common feature is ', length(common_feature)-3,' ---'))
  
  returnIDtoRS = function(rs.table.list, rawtableID){
    
    for (i in names(rs.table.list)) {
      rs.table.list[[i]] $ID = rawtableID[[i]]$ID
      rs.table.list[[i]] = rs.table.list[[i]] %>% dplyr::select('ID', everything())
    }
    
    return(rs.table.list)
  }
  
  

  if(!is.na(rf_nodesize)&
     !is.na(seed)&
     mode %in% c('all','single','double') &
     identical(c('ID','OS.time','OS'),colnames(train_data)[1:3])&
     length(candidate_genes)>0 &
     identical(c('ID','OS.time','OS'),common_feature[1:3]) &
     length(common_feature)>3 ){

    message('--- Data preprocessing ---')

    # Data preprocessing

    # Matching candidate genes to genes in each cohort
    common_feature = c('ID', 'OS.time', 'OS',candidate_genes)

    for (i in names(list_train_vali_Data)) {
      common_feature = intersect(common_feature, colnames(list_train_vali_Data[[i]]))
    }


    list_train_vali_Data <- lapply(list_train_vali_Data, function(x){x[, common_feature]})

    list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
      x[,-c(1:3)] <- apply(x[,-c(1:3)],2,as.numeric)
      return(x)})

    list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
      x[,c(1:2)] <- apply(x[,c(1:2)],2,as.factor)
      return(x)})

    list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
      x[,c(2:3)] <- apply(x[,c(2:3)],2,as.numeric)
      return(x)})

    list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
      x <- x[!is.na(x$OS.time)&!is.na(x$OS), ]
      return(x)})

    list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
      x <- x[x$OS.time > 0, ]
      return(x)})
    # use the mean replace the NA
    list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
      x[,-c(1:3)] <-apply(x[,-c(1:3)],2,function(x){
        x[is.na(x)]=mean(x,na.rm = T)
        return(x)
      })


      return(x)})

    train_data = train_data[,common_feature]
    train_data[,-c(1:3)] <- apply(train_data[,-c(1:3)],2,as.numeric)
    train_data[,c(1:2)] <- apply(train_data[,c(1:2)],2,as.factor)
    train_data[,c(2:3)] <- apply(train_data[,c(2:3)],2,as.numeric)

    if(is.null(unicox_p_cutoff)){
      unicox_p_cutoff =0.05
    } else {
      unicox_p_cutoff = unicox_p_cutoff
    }
    
    
    if(is.null(unicox.filter.for.candi)){
      unicox.filter.for.candi =T
    } else {
      unicox.filter.for.candi = unicox.filter.for.candi
    }
    
    if(unicox.filter.for.candi){
      
      cd.gene = common_feature[-c(1:3)]
      
      cd.gene1 = SigUnicox(gene_list = cd.gene,inputSet =train_data,unicox_pcutoff = unicox_p_cutoff )
      
      message(paste0('---the number of the final unicox filtered candidate genes is ', length(cd.gene1),' ---'))
      print(cd.gene1)
      cd.gene2 = c(common_feature[1:3],cd.gene1)
      common_feature =cd.gene2
      train_data =  as.data.frame(train_data)[,common_feature]
      
    } else{
      message(paste0('---the number of the final not unicox filtered candidate genes is ', length(common_feature)-3,' ---'))
      
      
    }


    est_dd <- as.data.frame(train_data)[, common_feature[-1]]
    val_dd_list <- lapply(list_train_vali_Data, function(x){x[, common_feature[-1]]})
    pre_var = common_feature[-c(1:3)]

    if(mode == 'all') {

      # return a dataframe  containing the C-index of models with different
      # machine learning algorithms in different queues

      result <- data.frame()
      ml.res = list()
      riskscore = list()

      # 1-1.RSF --------------------------------------------------------------
      message("---1-1 RSF ---")
      set.seed(seed)
      fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize,  #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)



      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- 'RSF'
      result <- rbind(result, cc)
      ml.res[[ 'RSF']] =fit
      riskscore[[ 'RSF']] = rs


      # 1-2.RSF + CoxBoost #########
      message("---1-2.RSF + CoxBoost ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize,  #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")

      rid <- rid$topvars

      if(length(rid)>1) {

        est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
        val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
        set.seed(seed)
        pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                                    trace=TRUE, start.penalty = 500, parallel = T)


        cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                              maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
        fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                        stepno = cv.res$optimal.step, penalty = pen$penalty)
        rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('RSF + ','CoxBoost')
        result <- rbind(result, cc)
        ml.res[[ paste0('RSF + ','CoxBoost')]] =fit

        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[ paste0('RSF + ','CoxBoost')]] = rs

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }


      # 1-3.RSF + Enet####

      message("---1-3.RSF + Enet ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars


      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      x1 <- as.matrix(est_dd2[, rid])
      x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
      for (alpha in seq(0.1, 0.9, 0.1)) {
        set.seed(seed)
        fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
        rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
        result <- rbind(result, cc)
        ml.res[[paste0('RSF + ', 'Enet', '[α=', alpha, ']')]] =fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[paste0('RSF + ', 'Enet', '[α=', alpha, ']')]] = rs

      }

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }

      # 1-4.RSF + GBM #####
      message("---1-4.RSF + GBM ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize,  #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars

      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                 n.trees = 10000,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10, n.cores = cores_for_parallel)

      # find index for number trees with minimum CV error
      best <- which.min(fit$cv.error)
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                 n.trees = best,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10,n.cores = cores_for_parallel)
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
      cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('RSF + ', 'GBM')
      result <- rbind(result, cc)
      ml.res[[ paste0('RSF + ', 'GBM')]] =  list( 'fit' = fit, 'best'= best)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('RSF + ', 'GBM')]] = rs

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }

      # 1-5.RSF + Lasso ####
      set.seed(seed)

      message("---1-5.RSF + Lasso ---")

      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars

      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      x1 <- as.matrix(est_dd2[, rid])
      x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})

      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('RSF + ', 'Lasso')
      result <- rbind(result, cc)

      ml.res[[ paste0('RSF + ', 'Lasso')]] =fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('RSF + ', 'Lasso')]] = rs

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }

      ## 1-6.RSF + plsRcox

      message("---1-6.RSF + plsRcox ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")



      rid <- rid$topvars

      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      set.seed(seed)
      cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)


      fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))

      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('RSF + ', 'plsRcox')
      result <- rbind(result, cc)

      ml.res[[ paste0('RSF + ', 'plsRcox')]] =fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('RSF + ', 'plsRcox')]] = rs

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }



      ## 1-7.RSF + Ridge

      message("---1-7.RSF + Ridge ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")

      rid <- rid$topvars

      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      x1 <- as.matrix(est_dd2[, rid])
      x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
      set.seed(seed)

      fit = cv.glmnet(x1, x2,
                      nfold=10, #10-fold cross-validation
                      family = "cox", alpha = 0,
                      type.measure = "class")
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('RSF + ', 'Ridge')
      result <- rbind(result, cc)

      ml.res[[ paste0('RSF + ', 'Ridge')]] =fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('RSF + ', 'Ridge')]] = rs


      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }


      ## 1-8.RSF + StepCox
      message("---1-8.RSF + StepCox ---")


      set.seed(seed)
      fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars


      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      for (direction in c("both", "backward", "forward")) {
        fit <- step(coxph(Surv(OS.time, OS)~., est_dd2), direction = direction)
        rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
        result <- rbind(result, cc)
        ml.res[[ paste0('RSF + ', 'StepCox', '[', direction, ']')]] =fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[  paste0('RSF + ', 'StepCox', '[', direction, ']')]] = rs

      }

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }



      ## 1-9.RSF + SuperPC

      message("---1-9.RSF + SuperPC ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, ##
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars

      if(length(rid)>1) {


      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
                   censoring.status = est_dd2$OS,
                   featurenames = colnames(est_dd2)[-c(1, 2)])
      set.seed(seed)
      fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
      repeat {
        tryCatch({
          cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                               n.fold = 10,
                               n.components = 3,
                               min.features = 2,
                               max.features = nrow(data$x),
                               compute.fullcv = TRUE,
                               compute.preval =TRUE)
          break
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          cat("Retrying...\n")
          Sys.sleep(1)
        })
      }
      rs <- lapply(val_dd_list2, function(w){
        test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status=w$OS, featurenames = colnames(w)[-c(1, 2)])
        ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
        rr <- as.numeric(ff$v.pred)
        rr2 <- cbind(w[, 1:2], RS = rr)
        return(rr2)
      })

      cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
        rownames_to_column('ID')
      cc$Model <- paste0('RSF + ', 'SuperPC')
      result <- rbind(result, cc)
      ml.res[[ paste0('RSF + ', 'SuperPC')]] = list(fit, cv.fit)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('RSF + ', 'SuperPC')]] = rs



      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }

      ## 1-10.RSF + survival-SVM

      message("---1-10.RSF + survival-SVM ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars
      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      fit = survivalsvm(Surv(OS.time, OS)~., data= est_dd2, gamma.mu = 1)
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('RSF + ', 'survival-SVM')
      result <- rbind(result,cc)

      ml.res[[ paste0('RSF + ', 'survival-SVM')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('RSF + ', 'survival-SVM')]] = rs


      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }
      # 2.Enet ---------------------------------------------------------------------
      message("---2.Enet ---")

      # 2.Enet
      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      for (alpha in seq(0.1, 0.9, 0.1)) {
        set.seed(seed)
        message(paste0('--- 2.Enet', '[α=', alpha, '] ---'))

        fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
        rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
        cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('Enet', '[α=', alpha, ']')
        result <- rbind(result, cc)
        ml.res[[ paste0('Enet', '[α=', alpha, ']')]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[  paste0('Enet', '[α=', alpha, ']')]] = rs


      }


      # 3.StepCox -------------------------------------------------------------
      message("---3.StepCox ---")


      # 3.StepCox
      for (direction in c("both", "backward", "forward")) {

        message( paste0('---3.StepCox', '[', direction, ']---'))

        fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction)



        rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('StepCox', '[', direction, ']')
        result <- rbind(result, cc)
        ml.res[[ paste0('StepCox', '[', direction, ']')]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[   paste0('StepCox', '[', direction, ']')]] = rs


      }


      #### direction = both ####
      direction = 'both'

      if(T){
        fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
        rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因

        if(length(rid)>1) {

          est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)
          message( paste0('---3.StepCox', '[', direction, ']', ' + CoxBoost ---'))

          pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                      trace=TRUE, start.penalty = 500, parallel = T)
          cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
          cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
          result <- rbind(result, cc)

          ml.res[[paste0('StepCox', '[', direction, ']', ' + CoxBoost')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[   paste0('StepCox', '[', direction, ']', ' + CoxBoost')]] = rs




          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          for (alpha in seq(0.1, 0.9, 0.1)) {
            set.seed(seed)
            message(  paste0('--- 3.StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, '] ---'))


            fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
            result <- rbind(result, cc)

            ml.res[[paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

            riskscore[[  paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')]] = rs


          }


          set.seed(seed)
          message( paste0('--- 3.StepCox', '[', direction, ']', ' + GBM ---'))


          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = 10000,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          # find index for number trees with minimum CV error
          best <- which.min(fit$cv.error)

          set.seed(seed)

          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = best,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
          cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction, ']', ' + GBM')]] =  list( 'fit' = fit, 'best'= best)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + GBM')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + Lasso ---'))
          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold=10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction, ']', ' + Lasso')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + Lasso')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + plsRcox ---'))
          set.seed(seed)
          cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
          fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                         event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + plsRcox')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + plsRcox')]] = rs


          message( paste0('--- 3.StepCox', '[', direction, ']', ' + Ridge ---'))

          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #-fold cross-validation
                          family = "cox", alpha = 0,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + Ridge')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + Ridge')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + RSF ---'))
          set.seed(seed)
          fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
                       ntree = 1000, nodesize = rf_nodesize, #
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)


          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
          result <- rbind(result, cc)
          ml.res[[  paste0('StepCox', '[', direction, ']', ' + RSF')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + RSF')]] = rs


          message( paste0('--- 3.StepCox', '[', direction, ']', ' + SuperPC ---'))

          data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
                       censoring.status = est_dd2$OS,
                       featurenames = colnames(est_dd2)[-c(1,2)])
          set.seed(seed)
          fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
          repeat {
            tryCatch({
              cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                   n.fold = 10,
                                   n.components = 3,
                                   min.features = 2,
                                   max.features = nrow(data$x),
                                   compute.fullcv = TRUE,
                                   compute.preval =TRUE)
              break
            }, error = function(e) {
              cat("Error:", conditionMessage(e), "\n")
              cat("Retrying...\n")
              Sys.sleep(1)
            })
          }
          rs <- lapply(val_dd_list2, function(w){
            test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
            ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
            rr <- as.numeric(ff$v.pred)
            rr2 <- cbind(w[,1:2], RS = rr)
            return(rr2)
          })
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + SuperPC')]] = list('fit'=fit,'cv.fit'=cv.fit)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[   paste0('StepCox', '[', direction, ']', ' + SuperPC')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + survival-SVM ---'))

          fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + survival-SVM')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + survival-SVM')]] = rs

        } else {
          warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
        }

      }

      #### direction = backward ####
      direction = 'backward'

      if(T){
        fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
        rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因

        if(length(rid)>1) {

          est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)
          message( paste0('---3.StepCox', '[', direction, ']', ' + CoxBoost ---'))

          pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                      trace=TRUE, start.penalty = 500, parallel = T)
          cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
          cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
          result <- rbind(result, cc)

          ml.res[[paste0('StepCox', '[', direction, ']', ' + CoxBoost')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[   paste0('StepCox', '[', direction, ']', ' + CoxBoost')]] = rs




          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          for (alpha in seq(0.1, 0.9, 0.1)) {
            set.seed(seed)
            message(  paste0('--- 3.StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, '] ---'))


            fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
            result <- rbind(result, cc)

            ml.res[[paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

            riskscore[[  paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')]] = rs


          }


          set.seed(seed)
          message( paste0('--- 3.StepCox', '[', direction, ']', ' + GBM ---'))


          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = 10000,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          # find index for number trees with minimum CV error
          best <- which.min(fit$cv.error)

          set.seed(seed)

          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = best,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
          cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction, ']', ' + GBM')]] =  list( 'fit' = fit, 'best'= best)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + GBM')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + Lasso ---'))
          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold=10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction, ']', ' + Lasso')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + Lasso')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + plsRcox ---'))
          set.seed(seed)
          cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
          fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                         event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + plsRcox')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + plsRcox')]] = rs


          message( paste0('--- 3.StepCox', '[', direction, ']', ' + Ridge ---'))

          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #-fold cross-validation
                          family = "cox", alpha = 0,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + Ridge')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + Ridge')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + RSF ---'))
          set.seed(seed)
          fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
                       ntree = 1000, nodesize = rf_nodesize, #
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
          result <- rbind(result, cc)
          ml.res[[  paste0('StepCox', '[', direction, ']', ' + RSF')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + RSF')]] = rs


          message( paste0('--- 3.StepCox', '[', direction, ']', ' + SuperPC ---'))

          data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
                       censoring.status = est_dd2$OS,
                       featurenames = colnames(est_dd2)[-c(1,2)])
          set.seed(seed)
          fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
          repeat {
            tryCatch({
              cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                   n.fold = 10,
                                   n.components = 3,
                                   min.features = 2,
                                   max.features = nrow(data$x),
                                   compute.fullcv = TRUE,
                                   compute.preval =TRUE)
              break
            }, error = function(e) {
              cat("Error:", conditionMessage(e), "\n")
              cat("Retrying...\n")
              Sys.sleep(1)
            })
          }
          rs <- lapply(val_dd_list2, function(w){
            test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
            ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
            rr <- as.numeric(ff$v.pred)
            rr2 <- cbind(w[,1:2], RS = rr)
            return(rr2)
          })
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + SuperPC')]] = list('fit'=fit,'cv.fit'=cv.fit)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[   paste0('StepCox', '[', direction, ']', ' + SuperPC')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + survival-SVM ---'))

          fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + survival-SVM')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + survival-SVM')]] = rs

        } else {
          warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
        }

      }


      #### direction = forward ####
      direction = 'forward'

      if(T){
        fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
        rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因

        if(length(rid)>1) {

          est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)
          message( paste0('---3.StepCox', '[', direction, ']', ' + CoxBoost ---'))

          pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                      trace=TRUE, start.penalty = 500, parallel = T)
          cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
          cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
          result <- rbind(result, cc)

          ml.res[[paste0('StepCox', '[', direction, ']', ' + CoxBoost')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[   paste0('StepCox', '[', direction, ']', ' + CoxBoost')]] = rs




          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          for (alpha in seq(0.1, 0.9, 0.1)) {
            set.seed(seed)
            message(  paste0('--- 3.StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, '] ---'))


            fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
            result <- rbind(result, cc)

            ml.res[[paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

            riskscore[[  paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')]] = rs


          }


          set.seed(seed)
          message( paste0('--- 3.StepCox', '[', direction, ']', ' + GBM ---'))


          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = 10000,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          # find index for number trees with minimum CV error
          best <- which.min(fit$cv.error)

          set.seed(seed)

          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = best,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
          cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction, ']', ' + GBM')]] =  list( 'fit' = fit, 'best'= best)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + GBM')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + Lasso ---'))
          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold=10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction, ']', ' + Lasso')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + Lasso')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + plsRcox ---'))
          set.seed(seed)
          cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
          fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                         event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + plsRcox')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + plsRcox')]] = rs


          message( paste0('--- 3.StepCox', '[', direction, ']', ' + Ridge ---'))

          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #-fold cross-validation
                          family = "cox", alpha = 0,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + Ridge')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + Ridge')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + RSF ---'))
          set.seed(seed)
          fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
                       ntree = 1000, nodesize = rf_nodesize, #
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
          result <- rbind(result, cc)
          ml.res[[  paste0('StepCox', '[', direction, ']', ' + RSF')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + RSF')]] = rs


          message( paste0('--- 3.StepCox', '[', direction, ']', ' + SuperPC ---'))

          data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
                       censoring.status = est_dd2$OS,
                       featurenames = colnames(est_dd2)[-c(1,2)])
          set.seed(seed)
          fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
          repeat {
            tryCatch({
              cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                   n.fold = 10,
                                   n.components = 3,
                                   min.features = 2,
                                   max.features = nrow(data$x),
                                   compute.fullcv = TRUE,
                                   compute.preval =TRUE)
              break
            }, error = function(e) {
              cat("Error:", conditionMessage(e), "\n")
              cat("Retrying...\n")
              Sys.sleep(1)
            })
          }
          rs <- lapply(val_dd_list2, function(w){
            test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
            ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
            rr <- as.numeric(ff$v.pred)
            rr2 <- cbind(w[,1:2], RS = rr)
            return(rr2)
          })
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + SuperPC')]] = list('fit'=fit,'cv.fit'=cv.fit)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[   paste0('StepCox', '[', direction, ']', ' + SuperPC')]] = rs



          message( paste0('--- 3.StepCox', '[', direction, ']', ' + survival-SVM ---'))

          fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction, ']', ' + survival-SVM')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[  paste0('StepCox', '[', direction, ']', ' + survival-SVM')]] = rs

        } else {
          warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
        }

      }

      # # 4.CoxBoost ------------------------------------------------------------
      message("---4.CoxBoost ---")

      ## 4-1.CoxBoost
      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)


      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
      cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost')
      result <- rbind(result, cc)

      ml.res[[ paste0('CoxBoost') ]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('CoxBoost')]] = rs

      ## 4-2.CoxBoost + Enet


      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)


      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)`!=0), "id"]

      if(length(rid)>1) {


      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      x1 <- as.matrix(est_dd2[, rid])
      x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
      for (alpha in seq(0.1, 0.9, 0.1)) {
        set.seed(seed)
        message( paste0('--- 4.CoxBoost', ' + Enet', '[α=', alpha, '] ---'))
        fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
        rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')
        result <- rbind(result, cc)
        ml.res[[ paste0('CoxBoost', ' + Enet', '[α=', alpha, ']') ]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[  paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')]] = rs

      }

      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }

      ## 4-3.CoxBoost + GBM
      set.seed(seed)
      message( paste0('--- 4.CoxBoost + ', 'GBM ---'))

      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)


      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)`!=0), "id"]

      if(length(rid)>1) {



      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                 n.trees = 10000,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10, n.cores = cores_for_parallel)
      # find index for number trees with minimum CV error
      best <- which.min(fit$cv.error)
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                 n.trees = best,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10,n.cores = cores_for_parallel)
      rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
      cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost + ', 'GBM')
      result <- rbind(result, cc)
      ml.res[[  paste0('CoxBoost + ', 'GBM')]] = list('fit'=fit, 'best'=best)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('CoxBoost + ', 'GBM')]] = rs


      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }

      ## 4-4.CoxBoost + Lasso
      set.seed(seed)

      message( paste0('--- 4.CoxBoost + ', 'Lasso ---'))

      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)



      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno=cv.res$optimal.step, penalty=pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]

      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      x1 <- as.matrix(est_dd2[, rid])
      x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost + ', 'Lasso')
      result <- rbind(result, cc)
      ml.res[[  paste0('CoxBoost + ', 'Lasso')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('CoxBoost + ', 'Lasso')]] = rs

      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }

      ## 4-5.CoxBoost + plsRcox
      set.seed(seed)
      message( paste0('--- 4.CoxBoost + ', 'plsRcox ---'))

      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)

      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]


      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      set.seed(seed)
      cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)

      fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))

      rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost + ', 'plsRcox')
      result <- rbind(result, cc)
      ml.res[[   paste0('CoxBoost + ', 'plsRcox')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('CoxBoost + ', 'plsRcox')]] = rs

      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }


      ## 4-6.CoxBoost + Ridge
      set.seed(seed)
      message( paste0('--- 4.CoxBoost + ', 'Ridge ---'))

      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)


      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]

      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      x1 <- as.matrix(est_dd2[, rid])
      x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold=10, #10-fold cross-validation
                      family = "cox", alpha = 0,
                      type.measure = "class")
      rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
      cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost + ', 'Ridge')
      result <- rbind(result, cc)

      ml.res[[   paste0('CoxBoost + ', 'Ridge')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('CoxBoost + ', 'Ridge')]] = rs


      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }


      ## 4-7.CoxBoost + StepCox
      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)


      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]

      if(length(rid)>1) {


      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      for (direction in c("both", "backward", "forward")) {

        message( paste0('--- 4.CoxBoost + ', 'StepCox', '[', direction, '] ---'))

        fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
        rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
        result <- rbind(result, cc)
        ml.res[[   paste0('CoxBoost + ', 'StepCox', '[', direction, ']')]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[ paste0('CoxBoost + ', 'StepCox', '[', direction, ']')]] = rs

      }
      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }


      ## 4-8.CoxBoost + SuperPC
      set.seed(seed)
    message(paste0('--- 4.CoxBoost + ', 'SuperPC ---'))

      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)



      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                            maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]

      if(length(rid)>1) {


      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
                   featurenames = colnames(est_dd2)[-c(1,2)])
      set.seed(seed)
      fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
      repeat {
        tryCatch({
          cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                               n.fold = 10,
                               n.components = 3,
                               min.features = 2,
                               max.features = nrow(data$x),
                               compute.fullcv = TRUE,
                               compute.preval =TRUE)
          break
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          cat("Retrying...\n")
          Sys.sleep(1)
        })
      }
      rs <- lapply(val_dd_list2, function(w){
        test <- list(x=t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
        ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
        rr <- as.numeric(ff$v.pred)
        rr2 <- cbind(w[,1:2], RS = rr)
        return(rr2)
      })
      cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost + ', 'SuperPC')
      result <- rbind(result, cc)

      ml.res[[   paste0('CoxBoost + ', 'SuperPC')]] = list('fit'=fit, 'cv.fit'=cv.fit)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('CoxBoost + ', 'SuperPC')]] = rs

      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }
      ## 4-9.CoxBoost + survival-SVM

      message(paste0('--- 4.CoxBoost + ', 'survival-SVM ---'))

      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)


      cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rid <- as.data.frame(coef(fit))
      rid$id <- rownames(rid)
      rid <- rid[which(rid$`coef(fit)` != 0), "id"]

      if(length(rid)>1) {


      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      fit = survivalsvm(Surv(OS.time, OS)~., data = est_dd2, gamma.mu = 1)
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
      result <- rbind(result, cc)
      ml.res[[  paste0('CoxBoost + ', 'survival-SVM')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[ paste0('CoxBoost + ', 'survival-SVM')]] = rs


      } else {
        warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
      }

      # 5.plsRcox ---------------------------------------------------------------
      message("---5.plsRcox ---")

      set.seed(seed)
      cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$OS.time, status = est_dd$OS), nt = 10, verbose = FALSE)

      fit <- plsRcox(est_dd[,pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))

      rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('plsRcox')
      result <- rbind(result, cc)
      ml.res[[   paste0('plsRcox')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('plsRcox')]] = rs

      # 6.superpc ---------------------------------------------------------------
      message("---6.superpc ---")

      data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 2)])
      set.seed(seed)
      fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
      repeat {
        tryCatch({
          cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                               n.fold = 10,
                               n.components = 3,
                               min.features = 2,
                               max.features = nrow(data$x),
                               compute.fullcv = TRUE,
                               compute.preval =TRUE)
          break
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          cat("Retrying...\n")
          Sys.sleep(1)
        })
      }
      rs <- lapply(val_dd_list, function(w){
        test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
        ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
        rr <- as.numeric(ff$v.pred)
        rr2 <- cbind(w[,1:2], RS = rr)
        return(rr2)
      })
      cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('SuperPC')
      result <- rbind(result, cc)
      ml.res[[   paste0('SuperPC')]] = list('fit'=fit, 'cv.fit'=cv.fit)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('SuperPC')]] = rs


      # 7.GBM -------------------------------------------------------------------
      message("---7.GBM ---")

      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd, distribution = 'coxph',
                 n.trees = 10000,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10, n.cores = cores_for_parallel)
      # find index for number trees with minimum CV error
      best <- which.min(fit$cv.error)
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd, distribution = 'coxph',
                 n.trees = best,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10, n.cores = cores_for_parallel)
      rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('GBM')
      result <- rbind(result, cc)
      ml.res[[   paste0('GBM')]] = list('fit' =fit , 'best' = best)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('GBM')]] = rs



      # 8.survivalsvm -----------------------------------------------------------
      message("---8.survivalsvm ---")

      fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd, gamma.mu = 1)
      rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('survival - SVM')
      result <- rbind(result, cc)

      ml.res[[ paste0('survival - SVM')]] = fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('survival - SVM')]] = rs




      # 9.Ridge -----------------------------------------------------------------
      message("---9.Ridge ---")


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
      cv.fit = cv.glmnet(x1, x2,
                        nfold = 10, #10-fold cross-validation
                        family = "cox",
                        type.measure = "class"
      )

      rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cv.fit$lambda.min)))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Ridge')
      result <- rbind(result, cc)

      ml.res[[ paste0('Ridge')]] =  list('fit'=fit, 'cv.fit'=cv.fit)

      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

      riskscore[[  paste0('Ridge')]] = rs

      # 10.Lasso ----------------------------------------------------------------
      message("---10.Lasso ---")

      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso')
      result <- rbind(result, cc)
      ml.res[[ paste0('Lasso')]] =  fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[  paste0('Lasso')]] = rs

      ## 10.1.Lasso + CoxBoost

      message( paste0('--- 10.Lasso + ','CoxBoost ---'))

      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min

      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      # rid <- rid[-1]

      if(length(rid)>1) {

      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      set.seed(seed)
      pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                  trace = TRUE, start.penalty = 500, parallel = T)

      cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                            maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
      fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                      stepno = cv.res$optimal.step, penalty = pen$penalty)
      rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
      cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso + ','CoxBoost')
      result <- rbind(result, cc)

      ml.res[[  paste0('Lasso + ','CoxBoost')]] =  fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[  paste0('Lasso + ','CoxBoost')]] = rs

      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }

      ## 10.2.Lasso + GBM

      message( paste0('--- 10.Lasso + ','GBM ---'))

      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #例文描述：10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min

      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      # rid <- rid[-1]
      if(length(rid)>1) {


      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                 n.trees = 10000,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10, n.cores = cores_for_parallel)
      # find index for number trees with minimum CV error
      best <- which.min(fit$cv.error)
      set.seed(seed)
      fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                 n.trees = best,
                 interaction.depth = 3,
                 n.minobsinnode = 10,
                 shrinkage = 0.001,
                 cv.folds = 10, n.cores = cores_for_parallel)
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
      cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso + ', 'GBM')
      result <- rbind(result, cc)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[  paste0('Lasso + ','GBM')]] = rs

      ml.res[[  paste0('Lasso + ', 'GBM')]] =  list('fit' =fit , 'best' = best)

      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }
      ## 10.3.Lasso + plsRcox

      message( paste0('--- 10.Lasso + ','plsRcox ---'))

      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min


      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      # rid <- rid[-1]
      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      set.seed(seed)
      cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)

      fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))


      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso + ', 'plsRcox')
      result <- rbind(result, cc)
      ml.res[[  paste0('Lasso + ', 'plsRcox')]] =  fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[  paste0('Lasso + ','plsRcox')]] = rs

      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }
      ## 10.4.Lasso + RSF

      message( paste0('--- 10.Lasso + ','RSF ---'))

      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min

      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      rid<-rid[-1]

      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      set.seed(seed)
      fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso', ' + RSF')
      result <- rbind(result, cc)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[  paste0('Lasso + ','RSF')]] = rs

      ml.res[[  paste0('Lasso', ' + RSF') ]] =  fit
      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }
      ## 10.5.Lasso + stepcox


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min



      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      # rid <- rid[-1]
      if(length(rid)>1) {


      est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      for (direction in c("both", "backward", "forward")) {

        message( paste0('--- 10.Lasso + ', 'StepCox', '[', direction, '] ---'))

        fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
        rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
        result <- rbind(result, cc)
        ml.res[[  paste0('Lasso + ', 'StepCox', '[', direction, ']') ]] =  fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[  paste0('Lasso + ', 'StepCox', '[', direction, ']')]] = rs

      }
      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }

      ## 10.6.Lasso + superPC
      message( paste0('--- 10.Lasso + ','superPC ---'))


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min

      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      # rid <- rid[-1]

      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
      data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
                   featurenames = colnames(est_dd2)[-c(1,2)])
      set.seed(seed)
      fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
      repeat {
        tryCatch({
          cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                               n.fold = 10,
                               n.components = 3,
                               min.features = 2,
                               max.features = nrow(data$x),
                               compute.fullcv = TRUE,
                               compute.preval =TRUE)
          break
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          cat("Retrying...\n")
          Sys.sleep(1)
        })
      }
      rs <- lapply(val_dd_list2, function(w){
        test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
        ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
        rr <- as.numeric(ff$v.pred)
        rr2 <- cbind(w[, 1:2], RS = rr)
        return(rr2)
      })
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso + ', 'SuperPC')
      result <- rbind(result, cc)

      ml.res[[   paste0('Lasso + ', 'SuperPC') ]] =  list('fit'=fit, 'cv.fit'=cv.fit)
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[   paste0('Lasso + ', 'SuperPC')]] = rs

      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }

      ## 10.7.Lasso + survival-SVM

      message( paste0('--- 10.Lasso + ','survival-SVM ---'))

      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      set.seed(seed)
      fit = cv.glmnet(x1, x2,
                      nfold = 10, #10-fold cross-validation
                      family = "cox", alpha = 1,
                      type.measure = "class")
      fit$lambda.min

      myCoefs <- coef(fit, s = "lambda.min");
      rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
      # rid <- rid[-1]

      if(length(rid)>1) {

      est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
      val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
      fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
      rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
      cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
        rownames_to_column('ID')
      cc$Model <- paste0('Lasso + ', 'survival-SVM')
      result <- rbind(result, cc)

      ml.res[[  paste0('Lasso + ', 'survival-SVM')]] =  fit
      rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
      
      riskscore[[   paste0('Lasso + ', 'survival-SVM')]] = rs


      } else {
        warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
      }

      # Return the result --------------------------------------------------------------

      result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

      return(result.cindex.fit)


    } else if (mode == 'single'

               ){ # return the trained model

      result <- data.frame()
      ml.res = list()
      riskscore = list()

      if(single_ml %in% 'RSF') {


        # -----RSF -------------------------------------------------------------------


        message('--- Start RSF ---')

        set.seed(seed)
        fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
                     ntree = 1000, nodesize = rf_nodesize,  #
                     splitrule = 'logrank',
                     importance = T,
                     proximity = T,
                     forest = T,
                     seed = seed)
        rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- 'RSF'
        result <- rbind(result, cc)
        ml.res[[ 'RSF']] =fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[ 'RSF']] = rs


        result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)

      } else if (single_ml %in%  'Enet'){

        # -----Enet -------------------------------------------------------------------

        message('--- Start Enet ---')

        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha_for_Enet, nfolds = 10) ### nfold is set as 10

          rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('Enet', '[α=', alpha_for_Enet, ']')
          result <- rbind(result, cc)
          ml.res[[ paste0('Enet', '[α=', alpha_for_Enet, ']')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('Enet', '[α=', alpha_for_Enet, ']')]] = rs

              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)



      } else if (single_ml %in% 'StepCox' & direction_for_stepcox %in% c("both", "backward", "forward")){
        # -----StepCox -------------------------------------------------------------------

        message('--- Start StepCox ---')
        set.seed(seed)

        fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction_for_stepcox)

        rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']')
        result <- rbind(result, cc)
        ml.res[[ paste0('StepCox', '[', direction_for_stepcox, ']')]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[   paste0('StepCox', '[', direction_for_stepcox, ']')]] = rs


            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)

      } else if (single_ml %in%  'CoxBoost') {
        # -----CoxBoost -------------------------------------------------------------------

        message('--- Start CoxBoost ---')
        set.seed(seed)
        pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                    trace = TRUE, start.penalty = 500, parallel = T)


        cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
        fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                        stepno = cv.res$optimal.step, penalty = pen$penalty)
        rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
        cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('CoxBoost')
        result <- rbind(result, cc)

        ml.res[[ paste0('CoxBoost') ]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[ paste0('CoxBoost')]] = rs

            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)

      } else if (single_ml %in%  'plsRcox') {

        # -----plsRcox -------------------------------------------------------------------

        message('--- Start plsRcox ---')

        set.seed(seed)
        cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$OS.time, status = est_dd$OS), nt = 10, verbose = T)

        fit <- plsRcox(est_dd[,pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))

        rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('plsRcox')
        result <- rbind(result, cc)
        ml.res[[   paste0('plsRcox')]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[  paste0('plsRcox')]] = rs

            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)


      } else if (single_ml %in% 'superpc'){
        # -----superpc -------------------------------------------------------------------

        message('--- Start superpc ---')
        set.seed(seed)

        data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 2)])
        set.seed(seed)
        fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
        repeat {
          tryCatch({
            cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                 n.fold = 10,
                                 n.components = 3,
                                 min.features = 2,
                                 max.features = nrow(data$x),
                                 compute.fullcv = TRUE,
                                 compute.preval =TRUE)
            break
          }, error = function(e) {
            cat("Error:", conditionMessage(e), "\n")
            cat("Retrying...\n")
            Sys.sleep(1)
          })
        }

        rs <- lapply(val_dd_list, function(w){
          test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
          ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
          rr <- as.numeric(ff$v.pred)
          rr2 <- cbind(w[,1:2], RS = rr)
          return(rr2)
        })
        cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('SuperPC')
        result <- rbind(result, cc)
        ml.res[[   paste0('SuperPC')]] = list('fit'=fit, 'cv.fit'=cv.fit)
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[  paste0('SuperPC')]] = rs


            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)



        } else if (single_ml %in%  'GBM'){
        # -----GBM -------------------------------------------------------------------

        message('--- Start GBM ---')

        set.seed(seed)
        fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd, distribution = 'coxph',
                   n.trees = 10000,
                   interaction.depth = 3,
                   n.minobsinnode = 10,
                   shrinkage = 0.001,
                   cv.folds = 10, n.cores = cores_for_parallel)
        # find index for number trees with minimum CV error
        best <- which.min(fit$cv.error)
        set.seed(seed)
        fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd, distribution = 'coxph',
                   n.trees = best,
                   interaction.depth = 3,
                   n.minobsinnode = 10,
                   shrinkage = 0.001,
                   cv.folds = 10, n.cores = cores_for_parallel)

        rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('GBM')
        result <- rbind(result, cc)
        ml.res[[   paste0('GBM')]] = list('fit' =fit , 'best' = best)
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[  paste0('GBM')]] = rs


            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)



      } else if (single_ml %in% 'survivalsvm'){
        # -----survivalsvm -------------------------------------------------------------------

        message('--- Start survivalsvm ---')

        set.seed(seed)
        fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd, gamma.mu = 1)
        rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('survival - SVM')
        result <- rbind(result, cc)

        ml.res[[ paste0('survival - SVM')]] = fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[  paste0('survival - SVM')]] = rs
            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)



      } else if (single_ml %in% 'Ridge'){
        # -----Ridge -------------------------------------------------------------------

        message('--- Start Ridge ---')


        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
        set.seed(seed)
        fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
        cv.fit = cv.glmnet(x1, x2,
                          nfold = 10, #10-fold cross-validation
                          family = "cox",
                          type.measure = "class"
        )


        rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cv.fit$lambda.min)))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('Ridge')
        result <- rbind(result, cc)

        ml.res[[ paste0('Ridge')]] =  list('fit'=fit, 'cv.fit'=cv.fit)
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        

        riskscore[[  paste0('Ridge')]] = rs

            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)


      } else if (single_ml %in%  'Lasso'){
        # -----Lasso -------------------------------------------------------------------

        message('--- Start Lasso ---')


        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
        set.seed(seed)
        fit = cv.glmnet(x1, x2,
                        nfold = 10, #10-fold cross-validation
                        family = "cox", alpha = 1,
                        type.measure = "class")
        rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('Lasso')
        result <- rbind(result, cc)
        ml.res[[ paste0('Lasso')]] =  fit
        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
        
        riskscore[[  paste0('Lasso')]] = rs

            result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

        return(result.cindex.fit)

      }
      # -----Finish and return the fit -------------------------------------------------------------------

      message('--- Finish and return the fit ---')



    } else if (mode == 'double' &
               double_ml1 %in% c('RSF', "StepCox","CoxBoost","Lasso")&
               double_ml2 %in% c('RSF',"Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")){

      # double_ml1 == RSF

      result <- data.frame()
      ml.res = list()
      riskscore = list()


      if(double_ml1 =='RSF'&
         double_ml2 %in% c("Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")) {

        if(double_ml2== 'Enet'& alpha_for_Enet %in% seq(0.1, 0.9, 0.1)){

          message(paste0('--- ',double_ml1,' + ',double_ml2, ' ',alpha_for_Enet,' ---' ))
          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")
          rid <- rid$topvars

          if(length(rid)>1) {


            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            x1 <- as.matrix(est_dd2[, rid])
            x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            for (alpha in seq(0.1, 0.9, 0.1)) {
              set.seed(seed)
              fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
              rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
              cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
                rownames_to_column('ID')
              cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha_for_Enet, ']')
              result <- rbind(result, cc)
              ml.res[[paste0('RSF + ', 'Enet', '[α=', alpha_for_Enet, ']')]] =fit
              rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
              
              riskscore[[paste0('RSF + ', 'Enet', '[α=', alpha_for_Enet, ']')]] = rs

                  result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

              return(result.cindex.fit)

            }


          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }



        } else if (double_ml2 =='StepCox' & direction_for_stepcox %in% c("both", "backward", "forward")){


          message(paste0('--- ',double_ml1,' + ',double_ml2, ' ',direction_for_stepcox,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")


          rid <- rid$topvars

          if(length(rid)>1) {


          est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            fit <- step(coxph(Surv(OS.time, OS)~., est_dd2), direction = direction_for_stepcox)


            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'StepCox', '[', direction_for_stepcox, ']')
            result <- rbind(result, cc)
            ml.res[[ paste0('RSF + ', 'StepCox', '[', direction_for_stepcox, ']')]] =fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('RSF + ', 'StepCox', '[', direction_for_stepcox, ']')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 == 'CoxBoost'){

          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")

          rid <- rid$topvars

          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            set.seed(seed)
            pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                                        trace=TRUE, start.penalty = 500, parallel = T)


            cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                                  maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
            fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                            stepno = cv.res$optimal.step, penalty = pen$penalty)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ','CoxBoost')
            result <- rbind(result, cc)
            ml.res[[ paste0('RSF + ','CoxBoost')]] =fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('RSF + ','CoxBoost')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)

          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        }else if (double_ml2 == 'plsRcox'){

          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")



          rid <- rid$topvars

          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            set.seed(seed)
            cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)


            fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))

            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'plsRcox')
            result <- rbind(result, cc)

            ml.res[[ paste0('RSF + ', 'plsRcox')]] =fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('RSF + ', 'plsRcox')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='superpc'){
          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")


          rid <- rid$topvars

          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
                         censoring.status = est_dd2$OS,
                         featurenames = colnames(est_dd2)[-c(1, 2)])
            set.seed(seed)
            fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
            repeat {
              tryCatch({
                cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                     n.fold = 10,
                                     n.components = 3,
                                     min.features = 2,
                                     max.features = nrow(data$x),
                                     compute.fullcv = TRUE,
                                     compute.preval =TRUE)
                break
              }, error = function(e) {
                cat("Error:", conditionMessage(e), "\n")
                cat("Retrying...\n")
                Sys.sleep(1)
              })
            }
            rs <- lapply(val_dd_list2, function(w){
              test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status=w$OS, featurenames = colnames(w)[-c(1, 2)])
              ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
              rr <- as.numeric(ff$v.pred)
              rr2 <- cbind(w[, 1:2], RS = rr)
              return(rr2)
            })

            cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'SuperPC')
            result <- rbind(result, cc)
            ml.res[[ paste0('RSF + ', 'SuperPC')]] = list(fit, cv.fit)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('RSF + ', 'SuperPC')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)
          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2=='GBM'){

          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")


          rid <- rid$topvars


          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                       n.trees = 10000,
                       interaction.depth = 3,
                       n.minobsinnode = 10,
                       shrinkage = 0.001,
                       cv.folds = 10, n.cores = cores_for_parallel)

            # find index for number trees with minimum CV error
            best <- which.min(fit$cv.error)
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                       n.trees = best,
                       interaction.depth = 3,
                       n.minobsinnode = 10,
                       shrinkage = 0.001,
                       cv.folds = 10,n.cores = cores_for_parallel)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
            cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'GBM')
            result <- rbind(result, cc)
            ml.res[[ paste0('RSF + ', 'GBM')]] =list( 'fit' = fit, 'best'= best)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('RSF + ', 'GBM')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)

          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='survivalsvm'){

          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")


          rid <- rid$topvars

          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            fit = survivalsvm(Surv(OS.time, OS)~., data= est_dd2, gamma.mu = 1)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'survival-SVM')
            result <- rbind(result,cc)

            ml.res[[ paste0('RSF + ', 'survival-SVM')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('RSF + ', 'survival-SVM')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 == 'Ridge'){
          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))


          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")

          rid <- rid$topvars
          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            x1 <- as.matrix(est_dd2[, rid])
            x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            set.seed(seed)

            fit = cv.glmnet(x1, x2,
                            nfold=10, #10-fold cross-validation
                            family = "cox", alpha = 0,
                            type.measure = "class")
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'Ridge')
            result <- rbind(result, cc)

            ml.res[[ paste0('RSF + ', 'Ridge')]] =fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('RSF + ', 'Ridge')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='Lasso'){

          message(paste0('--- ',double_ml1,' + ',double_ml2,' ---' ))

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                       ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rid <- var.select(object = fit, conservative = "high")


          rid <- rid$topvars
          if(length(rid)>1) {

            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            x1 <- as.matrix(est_dd2[, rid])
            x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2,
                            nfold = 10, #10-fold cross-validation
                            family = "cox", alpha = 1,
                            type.measure = "class")
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})

            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('RSF + ', 'Lasso')
            result <- rbind(result, cc)

            ml.res[[ paste0('RSF + ', 'Lasso')]] =fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('RSF + ', 'Lasso')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
          }

        }



      }else if(double_ml1 =='StepCox'&
               double_ml2 %in% c("CoxBoost",'Enet', 'RSF',"plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")&
               direction_for_stepcox %in% c("both", "backward", "forward")){

        if(double_ml2 =='CoxBoost'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))

          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

            set.seed(seed)
            message( paste0('---3.StepCox', '[', direction_for_stepcox, ']', ' + CoxBoost ---'))

            pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                        trace=TRUE, start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                  maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                            stepno = cv.res$optimal.step, penalty = pen$penalty)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + CoxBoost')
            result <- rbind(result, cc)

            ml.res[[paste0('StepCox', '[', direction_for_stepcox, ']', ' + CoxBoost')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[   paste0('StepCox', '[', direction_for_stepcox, ']', ' + CoxBoost')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }



        } else if(double_ml2 =='Enet'){

          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))

          if(length(rid)>1) {


          est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})


          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha_for_Enet, nfolds = 10)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + Enet', '[α=', alpha_for_Enet, ']')
            result <- rbind(result, cc)

            ml.res[[paste0('StepCox', '[', direction_for_stepcox, ']', ' + Enet', '[α=', alpha_for_Enet, ']')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + Enet', '[α=', alpha_for_Enet, ']')]] = rs
                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }


        }else if(double_ml2 =='RSF'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)
          fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
                       ntree = 1000, nodesize = rf_nodesize,
                       splitrule = 'logrank',
                       importance = T,
                       proximity = T,
                       forest = T,
                       seed = seed)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + RSF')
          result <- rbind(result, cc)
          ml.res[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + RSF')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + RSF')]] = rs

              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }

        }else if(double_ml2 =='plsRcox'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)


          cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
          fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                         event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + plsRcox')
          result <- rbind(result, cc)




          ml.res[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + plsRcox')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + plsRcox')]] = rs


          result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }


        }else if(double_ml2 =='superpc'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))
          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)
          data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
                       censoring.status = est_dd2$OS,
                       featurenames = colnames(est_dd2)[-c(1,2)])
          set.seed(seed)
          fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
          repeat {
            tryCatch({
              cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                   n.fold = 10,
                                   n.components = 3,
                                   min.features = 2,
                                   max.features = nrow(data$x),
                                   compute.fullcv = TRUE,
                                   compute.preval =TRUE)
              break
            }, error = function(e) {
              cat("Error:", conditionMessage(e), "\n")
              cat("Retrying...\n")
              Sys.sleep(1)
            })
          }
          rs <- lapply(val_dd_list2, function(w){
            test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
            ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
            rr <- as.numeric(ff$v.pred)
            rr2 <- cbind(w[,1:2], RS = rr)
            return(rr2)
          })
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + SuperPC')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + SuperPC')]] = list('fit'=fit,'cv.fit'=cv.fit)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[   paste0('StepCox', '[', direction_for_stepcox, ']', ' + SuperPC')]] = rs

              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }


        }else if(double_ml2 =='GBM'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          set.seed(seed)
          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = 10000,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)
          # find index for number trees with minimum CV error
          best <- which.min(fit$cv.error)
          set.seed(seed)
          fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
                     n.trees = best,
                     interaction.depth = 3,
                     n.minobsinnode = 10,
                     shrinkage = 0.001,
                     cv.folds = 10,n.cores = cores_for_parallel)

          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
          cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + GBM')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction_for_stepcox, ']', ' + GBM')]] = list('fit' =fit , 'best' = best)
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + GBM')]] = rs

              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }


        }else if(double_ml2 =='survivalsvm'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))
          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)

          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + survival-SVM')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + survival-SVM')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + survival-SVM')]] = rs


              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }

        }else if(double_ml2 =='Ridge'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))
          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})



          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #10-fold cross-validation
                          family = "cox", alpha = 0,
                          type.measure = "class")
          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + Ridge')
          result <- rbind(result, cc)

          ml.res[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + Ridge')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + Ridge')]] = rs


              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }

        }else if(double_ml2 =='Lasso'){
          set.seed(seed)
          message(paste0('--- ',double_ml1,' ', direction_for_stepcox ,' + ',double_ml2,' ---' ))

          fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction_for_stepcox)
          rid <- names(coef(fit))
          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})

          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold=10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")

          rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('StepCox', '[', direction_for_stepcox, ']', ' + Lasso')
          result <- rbind(result, cc)

          ml.res[[ paste0('StepCox', '[', direction_for_stepcox, ']', ' + Lasso')]] = fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
          
          riskscore[[  paste0('StepCox', '[', direction_for_stepcox, ']', ' + Lasso')]] = rs

              result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

          return(result.cindex.fit)

          } else {
            warning('The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2')
          }

        }




      } else if (double_ml1 =='CoxBoost'&
                 double_ml2 %in% c("Enet", "StepCox","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")){

        if( double_ml2 =='Enet' & alpha_for_Enet %in% seq(0.1, 0.9, 0.1)){
          message(paste0('--- ',double_ml1, ' + ',double_ml2, ' ',alpha_for_Enet ,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)


          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)`!=0), "id"]

          if(length(rid)>1) {


          est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
          x1 <- as.matrix(est_dd2[, rid])
          x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha_for_Enet, nfolds = 10)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha_for_Enet, ']')
            result <- rbind(result, cc)
            ml.res[[ paste0('CoxBoost', ' + Enet', '[α=', alpha_for_Enet, ']') ]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('CoxBoost', ' + Enet', '[α=', alpha_for_Enet, ']')]] = rs
                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)

          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='GBM'){
          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)


          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)`!=0), "id"]
          if(length(rid)>1) {


            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                       n.trees = 10000,
                       interaction.depth = 3,
                       n.minobsinnode = 10,
                       shrinkage = 0.001,
                       cv.folds = 10, n.cores = cores_for_parallel)
            # find index for number trees with minimum CV error
            best <- which.min(fit$cv.error)
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                       n.trees = best,
                       interaction.depth = 3,
                       n.minobsinnode = 10,
                       shrinkage = 0.001,
                       cv.folds = 10,n.cores = cores_for_parallel)
            rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
            cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'GBM')
            result <- rbind(result, cc)
            ml.res[[  paste0('CoxBoost + ', 'GBM')]] = list('fit' =fit , 'best' = best)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('CoxBoost + ', 'GBM')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='Lasso'){
          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)



          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                          stepno=cv.res$optimal.step, penalty=pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)` != 0), "id"]
          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            x1 <- as.matrix(est_dd2[, rid])
            x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2,
                            nfold = 10, #10-fold cross-validation
                            family = "cox", alpha = 1,
                            type.measure = "class")
            rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'Lasso')
            result <- rbind(result, cc)
            ml.res[[  paste0('CoxBoost + ', 'Lasso')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('CoxBoost + ', 'Lasso')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }


        } else if (double_ml2 =='plsRcox'){
          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)

          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)` != 0), "id"]
          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            set.seed(seed)
            cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)

            fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))

            rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'plsRcox')
            result <- rbind(result, cc)
            ml.res[[   paste0('CoxBoost + ', 'plsRcox')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('CoxBoost + ', 'plsRcox')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='Ridge'){


          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)


          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)` != 0), "id"]

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            x1 <- as.matrix(est_dd2[, rid])
            x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2,
                            nfold=10, #10-fold cross-validation
                            family = "cox", alpha = 0,
                            type.measure = "class")
            rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'Ridge')
            result <- rbind(result, cc)

            ml.res[[   paste0('CoxBoost + ', 'Ridge')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('CoxBoost + ', 'Ridge')]] = rs
                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='StepCox' & direction_for_stepcox %in% c("both", "backward", "forward")){


          message(paste0('--- ',double_ml1, ' + ',double_ml2, ' ', direction_for_stepcox ,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)


          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)` != 0), "id"]
          if(length(rid)>1) {


          est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction_for_stepcox)
            rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction_for_stepcox, ']')
            result <- rbind(result, cc)
            ml.res[[   paste0('CoxBoost + ', 'StepCox', '[', direction_for_stepcox, ']')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('CoxBoost + ', 'StepCox', '[', direction_for_stepcox, ']')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }


        } else if (double_ml2 =='superpc'){
          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          set.seed(seed)
          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)



          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                                maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)` != 0), "id"]
          if(length(rid)>1) {


            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
                         featurenames = colnames(est_dd2)[-c(1,2)])
            set.seed(seed)
            fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
            repeat {
              tryCatch({
                cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                     n.fold = 10,
                                     n.components = 3,
                                     min.features = 2,
                                     max.features = nrow(data$x),
                                     compute.fullcv = TRUE,
                                     compute.preval =TRUE)
                break
              }, error = function(e) {
                cat("Error:", conditionMessage(e), "\n")
                cat("Retrying...\n")
                Sys.sleep(1)
              })
            }
            rs <- lapply(val_dd_list2, function(w){
              test <- list(x=t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
              ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
              rr <- as.numeric(ff$v.pred)
              rr2 <- cbind(w[,1:2], RS = rr)
              return(rr2)
            })
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'SuperPC')
            result <- rbind(result, cc)

            ml.res[[   paste0('CoxBoost + ', 'SuperPC')]] = list('fit'=fit, 'cv.fit'=cv.fit)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('CoxBoost + ', 'SuperPC')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }



        } else if (double_ml2 =='survivalsvm'){
          set.seed(seed)

          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                      trace = TRUE, start.penalty = 500, parallel = T)


          cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                                maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
          fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                          stepno = cv.res$optimal.step, penalty = pen$penalty)
          rid <- as.data.frame(coef(fit))
          rid$id <- rownames(rid)
          rid <- rid[which(rid$`coef(fit)` != 0), "id"]
          if(length(rid)>1) {


            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            fit = survivalsvm(Surv(OS.time, OS)~., data = est_dd2, gamma.mu = 1)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
            result <- rbind(result, cc)
            ml.res[[  paste0('CoxBoost + ', 'survival-SVM')]] = fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[ paste0('CoxBoost + ', 'survival-SVM')]] = rs


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)

          } else {
            warning('The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2')
          }

        }



      } else if (double_ml1 =='Lasso'&
                 double_ml2 %in% c( "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","RSF")){

        if(double_ml2 =='CoxBoost'){

          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min

          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          # rid <- rid[-1]
          if(length(rid)>1) {


            est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            set.seed(seed)
            pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                        trace = TRUE, start.penalty = 500, parallel = T)

            cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                                  maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                            stepno = cv.res$optimal.step, penalty = pen$penalty)
            rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
            cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso + ','CoxBoost')
            result <- rbind(result, cc)

            ml.res[[  paste0('Lasso + ','CoxBoost')]] =  fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('Lasso + ','CoxBoost')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='GBM'){

          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min

          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          # rid <- rid[-1]

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                       n.trees = 10000,
                       interaction.depth = 3,
                       n.minobsinnode = 10,
                       shrinkage = 0.001,
                       cv.folds = 10, n.cores = cores_for_parallel)
            # find index for number trees with minimum CV error
            best <- which.min(fit$cv.error)
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
                       n.trees = best,
                       interaction.depth = 3,
                       n.minobsinnode = 10,
                       shrinkage = 0.001,
                       cv.folds = 10, n.cores = cores_for_parallel)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
            cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso + ', 'GBM')
            result <- rbind(result, cc)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('Lasso + ','GBM')]] = rs

            ml.res[[  paste0('Lasso + ', 'GBM')]] =  list('fit' =fit , 'best' = best)

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='plsRcox'){

          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min


          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          # rid <- rid[-1]
          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            set.seed(seed)
            cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)

            fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))


            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso + ', 'plsRcox')
            result <- rbind(result, cc)
            ml.res[[  paste0('Lasso + ', 'plsRcox')]] =  fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('Lasso + ','plsRcox')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='RSF'){


          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #例文描述：10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min

          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          rid<-rid[-1]

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
                         ntree = 1000, nodesize = rf_nodesize, #
                         splitrule = 'logrank',
                         importance = T,
                         proximity = T,
                         forest = T,
                         seed = seed)
            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso', ' + RSF')
            result <- rbind(result, cc)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('Lasso + ','RSF')]] = rs

            ml.res[[  paste0('Lasso', ' + RSF') ]] =  fit


                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='StepCox' & direction_for_stepcox %in% c("both", "backward", "forward")){


          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ', direction_for_stepcox,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min



          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          # rid <- rid[-1]

          if(length(rid)>1) {

          est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
          val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction_for_stepcox)

            rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction_for_stepcox, ']')
            result <- rbind(result, cc)
            ml.res[[  paste0('Lasso + ', 'StepCox', '[', direction_for_stepcox, ']') ]] =  fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[  paste0('Lasso + ', 'StepCox', '[', direction_for_stepcox, ']')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }

        } else if (double_ml2 =='superpc'){
          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, # 10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min

          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          # rid <- rid[-1]

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
            data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
                         featurenames = colnames(est_dd2)[-c(1,2)])
            set.seed(seed)
            fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
            repeat {
              tryCatch({
                cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                                     n.fold = 10,
                                     n.components = 3,
                                     min.features = 2,
                                     max.features = nrow(data$x),
                                     compute.fullcv = TRUE,
                                     compute.preval =TRUE)
                break
              }, error = function(e) {
                cat("Error:", conditionMessage(e), "\n")
                cat("Retrying...\n")
                Sys.sleep(1)
              })
            }
            rs <- lapply(val_dd_list2, function(w){
              test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
              ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
              rr <- as.numeric(ff$v.pred)
              rr2 <- cbind(w[, 1:2], RS = rr)
              return(rr2)
            })
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso + ', 'SuperPC')
            result <- rbind(result, cc)

            ml.res[[   paste0('Lasso + ', 'SuperPC') ]] =  list('fit'=fit, 'cv.fit'=cv.fit)
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[   paste0('Lasso + ', 'SuperPC')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)


          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }


        } else if (double_ml2 =='survivalsvm'){


          message(paste0('--- ',double_ml1, ' + ',double_ml2,' ---' ))

          x1 <- as.matrix(est_dd[, pre_var])
          x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
          set.seed(seed)
          fit = cv.glmnet(x1, x2,
                          nfold = 10, #例文描述：10-fold cross-validation
                          family = "cox", alpha = 1,
                          type.measure = "class")
          fit$lambda.min

          myCoefs <- coef(fit, s = "lambda.min");
          rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
          # rid <- rid[-1]

          if(length(rid)>1) {

            est_dd2 <- train_data[,c('OS.time', 'OS', rid)]
            val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[,c('OS.time', 'OS', rid)]})
            fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
            rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
            cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
              rownames_to_column('ID')
            cc$Model <- paste0('Lasso + ', 'survival-SVM')
            result <- rbind(result, cc)

            ml.res[[  paste0('Lasso + ', 'survival-SVM')]] =  fit
            rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
            
            riskscore[[   paste0('Lasso + ', 'survival-SVM')]] = rs

                result.cindex.fit = list('Cindex.res'=result, 'ml.res'=ml.res, 'riskscore' = riskscore, 'Sig.genes' = pre_var)

            return(result.cindex.fit)



          } else {
            warning('The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2')
          }




        }






      }




}






 message("--- The analysis has been completed ---")


  } else {

    print('Please provide the full parameters, verify that the column names of the supplied queues match the conditions
          (ID, OS.time, OS for the first through third columns, respectively) ,  verify there exist common genes in
          all cohorts')

  }




}






















