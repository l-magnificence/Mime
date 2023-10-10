


# 1.0 利用ml建立预后模型 ----------------------------------------------------------

setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res")
# dir.create('1.Prog.Model')
setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model")
rm(list = ls())

##

 
####################### 建立预后模型（张炜，已经完成） ###################################################
# source('/export3/zhangw/Project_Cross/Project_Mime/Proj/code/Prognostic.model.con.R')

####################### Cindex 可视化（刘宏伟 已经完成） ###################################################

load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")
source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

source('/export3/zhangw/Project_Cross/Project_Mime/Proj/data/Glioma.cohort.R')
list_train_vali_Data = list(TCGA= sur.matrix.TCGA.Glioma,
                            CGGA.325 = sur.matrix.CGGA_RNAseq_325_FPKM_Glioma,
                            CGGA.693 = sur.matrix.CGGA_RNAseq_693_FPKM_Glioma,
                            CGGA.1018 = sur.matrix.CGGA_RNAseq_1018_FPKM_Glioma,
                            CGGA.array = sur.matrix.CGGA_array_Glioma,
                            GLASS_TP = sur.matrix.GLASS_TP_Glioma,
                            GLASS_R1 = sur.matrix.GLASS_R1_Glioma,
                            GSE108474 = sur.matrix.GSE108474_Glioma,
                            GSE16011 = sur.matrix.GSE16011_Glioma,
                            GSE43289 =sur.matrix.GSE43289_Glioma,
                            GSE7696 = sur.matrix.GSE7696_Glioma
)
save(list_train_vali_Data,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cindex_dis_all.pdf',width = 10,height = 15,onefile = F)
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.2)
dev.off()

### RSF + survival−SVM 为最佳 的model
# 单独展示RSF + survival-SVM 在所有队列中的cindex

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cindex_specific_model.pdf',width = 5,height = 5.5,onefile = F)
cindex_dis_select(res,
                  model="RSF + survival-SVM",
                  order= names(list_train_vali_Data))
dev.off()

rm(list = ls())

############# compared with other clinical and molecular variables in predicting prognosis  #### 
## Not applicable so far
## PMID: 35145098 参考文献

####################### KM曲线 生存分析 （刘宏伟 已经完成）  ###################################################

### KM曲线 生存分析 
# 将 RSF + survival−SVM 生存分析的km曲线在所有11个队列中，直接提取risk score table就行

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")

survplot <- vector("list",11) 
for (i in c(1:11)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "RSF + survival-SVM",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}

library(patchwork)

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/sur_km.pdf',width = 15,height = 20,onefile = F)
survplot[[1]]+survplot[[2]]+survplot[[3]]+survplot[[4]]+survplot[[5]]+survplot[[6]]+
  survplot[[7]]+survplot[[8]]+survplot[[9]]+survplot[[10]]+survplot[[11]]+
  plot_layout(ncol = 3)
dev.off()

####################### 计算auc （张炜 刘宏伟 已经完成） ###################################################

source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_AUC_ml_res.R')

all.auc.1y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1)
save(all.auc.1y,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.1y.Rdata")

all.auc.3y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3)
save(all.auc.3y,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.3y.Rdata")

all.auc.5y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5)
save(all.auc.5y,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.5y.Rdata")

rm(list = ls())
#######################  auc 可视化（刘宏伟 已经完成） ###################################################

#
##
### auc 可视化
# 类似于C-index在所有model在所有队列中的表达，

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Monamdel/all.auc.1y.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.3y.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.5y.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc1y_dis_all.pdf',width = 10,height = 15,onefile = F)
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.2,
            year=1)
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc3y_dis_all.pdf',width = 10,height = 15,onefile = F)
auc_dis_all(all.auc.3y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.2,
            year=3)
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc5y_dis_all.pdf',width = 10,height = 15,onefile = F)
auc_dis_all(all.auc.5y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.2,
            year=5)
dev.off()

# 或者展示RSF + survival-SVM 在所有队列中的auc， 这里可能没有3年或者5年auc，队列中生存时间不够
## PMID: 35145098 参考文献

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc_specific_model.pdf',width = 10,height = 5,onefile = F)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="RSF + survival-SVM",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
dev.off()

####################### roc 可视化 （刘宏伟 已经完成）  ###################################################

#
##
### roc 可视化 
# 展示RSF + survival-SVM 在所有队列中？

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc1y_roc.pdf',width = 7.5,height = 7.0,onefile = F)
roc_vis(all.auc.1y,
        model_name = "RSF + survival-SVM",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc3y_roc.pdf',width = 7.5,height = 7.0,onefile = F)
roc_vis(all.auc.3y,
        model_name = "RSF + survival-SVM",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=3)
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc5y_roc.pdf',width = 7.5,height = 7.0,onefile = F)
roc_vis(all.auc.5y,
        model_name = "RSF + survival-SVM",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=5)
dev.off()

####################### 将RSF+SuperPC 的结果进行meta 分析 （张炜）   ################################

## 将RSF+SuperPC 的结果进行meta 分析

# step1 计算RS 的单因素回归 结果

mm = res$riskscore$`RSF + SuperPC`


unicox.cal = lapply(mm,function(x){

  tmp <- x
  cox <- coxph(Surv(OS.time, OS) ~ RS, data = tmp)
  coxSummary <- summary(cox)
  
  unicox = c( as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
              as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
              as.numeric(coxSummary$conf.int[,3][1]),
              as.numeric(coxSummary$conf.int[,4][1])
  )
  
  return(unicox)})

unicox.rs.res = unicox.cal  %>% do.call(rbind,.)
colnames(unicox.rs.res) = c('HR','pvalue',"LCI",'HCI')




#
## 
### 单因素回归结果 可视化美化，如果可以的话？ （刘宏伟）



# step2 meta analysis
library(meta)
input <- unicox.rs.res %>% as.data.frame() %>%
  mutate(`Hazard Ratio(95%CI)` = paste(HR,'(',LCI,'-',HCI,')',sep=""))
input$Group = rownames(input)
# 对不服从正态分布的HR对数转换
lnhr <- log(input$HR)
lnlci <- log(input$LCI)
lnhci <- log(input$HCI)
selnhr <- (lnhci-lnlci)/(2*1.96)
metamodel = metagen(TE = lnhr, seTE = selnhr,
                    sm="HR",            # 待合并的效应量
                    data=input,  # 输入数据
                    studlab=Group)     # 研究标签


pdf("meta_random_rs.rsf.superpc.pdf", width = 12, height = 5)
forest(metamodel,
       comb.fixed = F, comb.random = T,
       layout = 'revman5') # 套用RevMan 5风格
dev.off()

####################### 将RSF + survival-SVM 的结果进行多因素回归  （刘宏伟） ################################

#
##
### 多因素回归

# 在TCGA, CGGA325,CGGA693, CGGA1018等中进行多因素回归， 变量的选择？



####################### 将RSF+SuperPC 的结果和已经发表的文章的signature进行比较 （张炜，刘宏伟，张益浩）  ################################
### 这个需要等张益浩收集完所有的signature， 主要是张益浩，其余人辅助收集
### 张炜负责将signature 计算成riskscore
### 刘宏伟 负责可视化


## PMID: 35145098 参考文献

# 热图展示所有 signature的预后关系

# 比较 Cindex, auc(1,3,5?), 单因素回归的HR和p值，某几个队列中多因素回归的结果



####################### 就预后模型这一块还有什么要补充的吗 ########################################################
### 欢迎补充





# 2.0 利用1.0中建立的RSF+SuperPC预后模型进行免疫浸润的分析（张益浩） -----------------------------------------------------

### 这里应该是同一套代码，在TCGA, CGGA325, CGGA693, CGGA1018 四个队列中都跑一遍 



 
# 3.0 利用1.0中建立的RSF+SuperPC预后模型进行富集的分析（刘宏伟） -----------------------------------------------------



# 4.0 分类变量的预测模型的建立--这里是预测免疫治疗反应与否(张炜)  -----------------------------------------------------



# 5.0 核心feature的选择（和预后相关）(张炜已完成)  -----------------------------------------------------

setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res")
# dir.create('5.CoreFeature')
setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/5.CoreFeature")




































































