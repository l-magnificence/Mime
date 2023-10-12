


# 1.0 利用ml建立预后模型 ----------------------------------------------------------

setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res")
# dir.create('1.Prog.Model')
setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model")
rm(list = ls())

##

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
####################### auc 可视化（刘宏伟 已经完成） ###################################################

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

####################### 将optimal.model 的结果进行meta 分析 （张炜 已经完成）   ################################

## 将optimal.model 的结果进行meta 分析

# step1 计算RS 的单因素回归 结果

source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_unicox_ml_res.R')

optimal.model = 'RSF + survival-SVM'

mm = res[["riskscore"]][[optimal.model]]


unicox.rs.res = cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = optimal.model)



#
## 
### 单因素回归结果 可视化美化，如果可以的话？ （刘宏伟）



# step2 meta analysis
source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_unicox_meta_ml_res.R')

metamodel = cal_unicox_meta_ml_res(input = unicox.rs.res)


pdf("meta_random_rs.rsf.superpc.pdf", width = 12, height = 5)
forest(metamodel,
       comb.fixed = F, comb.random = T,
       layout = 'revman5') # 套用RevMan 5风格
dev.off()

####################### 将optimal.model的结果进行多因素回归  （刘宏伟） ################################

#
##
### 多因素回归

# 在TCGA, CGGA325,CGGA693, CGGA1018等中进行多因素回归， 变量的选择？



####################### 将optimal.model的结果和已经发表的文章的signature进行比较 （张炜，刘宏伟，张益浩）  ################################
####################### 收集完所有glioma 相关的signature（张益浩已完成）#############################
# load('/export3/zhangw/Project_Cross/Project_Mime/data/sig/glioma.sig.Rdata')
# load('/export3/zhangw/Project_Cross/Project_Mime/data/sig/lgg.sig.Rdata')
# load('/export3/zhangw/Project_Cross/Project_Mime/data/sig/gbm.sig.Rdata')

#######################  将signature 计算成riskscore （张炜已完成）####################################
# rm(list = ls())
# gc()
# source('/export3/zhangw/Project_Cross/Project_Mime/Proj/data/Glioma.cohort.R')
# 
# list_train_vali_Data = list(TCGA= sur.matrix.TCGA.Glioma,
#                             CGGA.325 = sur.matrix.CGGA_RNAseq_325_FPKM_Glioma,
#                             CGGA.693 = sur.matrix.CGGA_RNAseq_693_FPKM_Glioma,
#                             CGGA.1018 = sur.matrix.CGGA_RNAseq_1018_FPKM_Glioma,
#                             CGGA.array = sur.matrix.CGGA_array_Glioma,
#                             GLASS_R1 = sur.matrix.GLASS_R1_Glioma,
#                             GLASS_TP = sur.matrix.GLASS_TP_Glioma,
#                             GSE108474 = sur.matrix.GSE108474_Glioma,
#                             GSE16011 = sur.matrix.GSE16011_Glioma,
#                             GSE43289 =sur.matrix.GSE43289_Glioma,
#                             GSE7696 = sur.matrix.GSE7696_Glioma
# )
load("~/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

source("/export3/zhangw/Project_Cross/Project_Mime/Function/cal_RS_pre.prog.sig.R")
# rs.glioma = cal_RS_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data)
# rs.gbm = cal_RS_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data)
# rs.lgg = cal_RS_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data)
# rs.lgg.GBM = cal_RS_pre.prog.sig(type.sig = c('LGG','GBM'),list_input_data = list_train_vali_Data)
rs.glioma.lgg.gbm = cal_RS_pre.prog.sig(type.sig = c('LGG','GBM','Glioma'),list_input_data = list_train_vali_Data)
save(rs.glioma.lgg.gbm,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/rs.glioma.lgg.gbm.Rdata")

source("/export3/zhangw/Project_Cross/Project_Mime/Function/cal_cindex_pre.prog.sig.R")
# cc.gbm = cal_cindex_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data)
# cc.lgg = cal_cindex_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data)
# cc.glioma= cal_cindex_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data)
cc.glioma.lgg.gbm = cal_cindex_pre.prog.sig(type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data)
save(cc.glioma.lgg.gbm,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cc.glioma.lgg.gbm.Rdata")

source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_auc_pre.prog.sig.R')

auc.glioma.lgg.gbm.1 = cal_auc_pre.prog.sig(type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data,AUC_time = 1)
save(auc.glioma.lgg.gbm.1,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.1.Rdata")

auc.glioma.lgg.gbm.3 = cal_auc_pre.prog.sig(type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data,AUC_time = 3)
save(auc.glioma.lgg.gbm.3,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.3.Rdata")

auc.glioma.lgg.gbm.5 = cal_auc_pre.prog.sig(type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data,AUC_time = 5)
save(auc.glioma.lgg.gbm.5,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.5.Rdata")

rm(list = ls())
# auc.glioma.1 = cal_auc_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data,AUC_time = 1)
# auc.glioma.3 = cal_auc_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data,AUC_time = 3)
# auc.glioma.5 = cal_auc_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data,AUC_time = 5)
# 
# 
# auc.gbm.1 = cal_auc_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data,AUC_time = 1)
# auc.gbm.3 = cal_auc_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data,AUC_time = 3)
# auc.gbm.5 = cal_auc_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data,AUC_time = 5)
# 
# 
# auc.lgg.1 = cal_auc_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data,AUC_time = 1)
# auc.lgg.3 = cal_auc_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data,AUC_time = 3)
# auc.lgg.5 = cal_auc_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data,AUC_time = 5)

########################## 可视化（刘宏伟 已经完成）########################################

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

load("~/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

## PMID: 35145098 参考文献
# 热图展示所有 signature的预后关系, 单因素回归的HR和p值
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/rs.glioma.lgg.gbm.Rdata")
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/hr_comp.pdf',width = 36,height = 6,onefile = F)
HR_com(rs.glioma.lgg.gbm,
       res,
       model_name="RSF + survival-SVM",
       dataset=names(list_train_vali_Data))
dev.off()

# 比较 Cindex
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/cc.glioma.lgg.gbm.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cindex_comp.pdf',width = 32,height = 13,onefile = F)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="RSF + survival-SVM",
            dataset=names(list_train_vali_Data))
dev.off()

# 比较 auc(1,3,5)
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.1.Rdata")
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.3.Rdata")
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.5.Rdata")
Load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.1y.Rdata")
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.3y.Rdata")
load("~/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.5y.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc_comp_1y.pdf',width = 32,height = 13,onefile = F)
auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name="RSF + survival-SVM",
         dataset=names(list_train_vali_Data))
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc_comp_3y.pdf',width = 32,height = 13,onefile = F)
auc_comp(auc.glioma.lgg.gbm.3,
         all.auc.3y,
         model_name="RSF + survival-SVM",
         dataset=names(list_train_vali_Data))
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc_comp_5y.pdf',width = 32,height = 13,onefile = F)
auc_comp(auc.glioma.lgg.gbm.5,
         all.auc.5y,
         model_name="RSF + survival-SVM",
         dataset=names(list_train_vali_Data))
dev.off()

####################### 就预后模型这一块还有什么要补充的吗 ########################################################
### 欢迎补充





# 2.0 利用1.0中建立的optimal.model预后模型进行免疫浸润的分析（张益浩） -----------------------------------------------------

### 这里应该是同一套代码，在TCGA, CGGA325, CGGA693, CGGA1018 四个队列中都跑一遍 



 
# 3.0 利用1.0中建立的optimal.model预后模型进行富集的分析（张益浩） -----------------------------------------------------



# 4.0 分类变量的预测模型的建立--这里是预测免疫治疗反应与否(张炜计算,已经完成, 可视化刘宏伟)  -----------------------------------------------------

setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/4.ICI_response")

#######################利用ml建立预测免疫治疗反应的model (张炜 已经完成)############################################
# source('/export3/zhangw/Project_Cross/Project_Mime/Proj/code/ICI.response.model.R')

load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/4.ICI_response/ICIresponse.model.all.p1.Rdata")

#######################可视化建立的预测ICI反应的model (刘宏伟)############################################



# 5.0 核心feature的选择（和预后相关）(计算张炜，已完成， 可视化刘宏伟)  -----------------------------------------------------

#######################核心feature的选择（和预后相关) (张炜 已完成)############################################


setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res")
# dir.create('5.CoreFeature')
setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/5.CoreFeature")

# 计算
# source('/export3/zhangw/Project_Cross/Project_Mime/Proj/code/ML.CoreFeature.R')

load('/export3/zhangw/Project_Cross/Project_Mime/Proj/res/5.CoreFeature/feature.all.res.Rdata')

#######################可视化 (刘宏伟)################################################################

# 6.0 在核心feature中选择某一个关键基因进行单基因分析的一整套(预后,差异,回归,龙哥的paper那一套)（刘宏伟） ----------------------------------------------------



































































