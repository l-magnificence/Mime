


# 1.0 利用ml建立预后模型 ----------------------------------------------------------

# setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res")
# # dir.create('1.Prog.Model')
# setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model")
# rm(list = ls())

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

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/sur_km.pdf',width = 10,height = 28,onefile = F)
survplot[[1]]+survplot[[2]]+survplot[[3]]+survplot[[4]]+survplot[[5]]+survplot[[6]]+
  survplot[[7]]+survplot[[8]]+survplot[[9]]+survplot[[10]]+survplot[[11]]+
  plot_layout(ncol = 2)
dev.off()

####################### 计算auc （张炜 刘宏伟 已经完成） ###################################################

source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_AUC_ml_res.R')

all.auc.1y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                            auc_cal_method="KM")
save(all.auc.1y,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.1y_km.Rdata")

all.auc.3y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                            auc_cal_method="KM")
save(all.auc.3y,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.3y_km.Rdata")

all.auc.5y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                            auc_cal_method="KM")
save(all.auc.5y,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.5y_km.Rdata")

rm(list = ls())
####################### auc 可视化（刘宏伟 已经完成） ###################################################

### auc 可视化
# 类似于C-index在所有model在所有队列中的表达，

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.1y_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.3y_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.5y_km.Rdata")
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
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")
optimal.model = 'RSF + survival-SVM'

unicox.rs.res = cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = optimal.model,type ='categorical')
# step2 meta analysis
source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_unicox_meta_ml_res.R')

metamodel = cal_unicox_meta_ml_res(input = unicox.rs.res)
save(metamodel,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/metamodel.Rdata")

####################### 单因素回归结果 可视化（刘宏伟 已经完成） #######################

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/metamodel.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

tiff("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/meta_rs_rsf_survivalSVM.tiff", units ="in",
     width = 9.0, height = 4.5,res = 600,compression ='zip')
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
dev.off()

rm(list = ls())

####################### 将optimal.model的结果进行多因素回归  （刘宏伟 已经完成） ################################

library(ezcox)
library(forestploter)
library(grid)
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")
# load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")
optimal.model = 'RSF + survival-SVM'
multic_cox<-data.frame()
multic_cox_model<-list()

### 多因素回归
# TCGA
datasets<-"TCGA"
meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")
rs <- res[["riskscore"]][[optimal.model]][[datasets]]
rownames(rs)<-rs$ID
rownames(rs)<-gsub("\\.","-",rownames(rs))
meta<-meta[rownames(rs),]
meta<-cbind(meta,rs)

data<-dplyr::select(meta,c("Age (years at diagnosis)","Gender","Grade","IDH status","1p/19q codeletion","MGMT promoter status",
                           "OS","OS.time","RS"))
colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                  "status","time","RS")
data$Age<-as.numeric(data$Age)
data<-data[data$IDH !="NA",]
data<-data[data$`Chr1p/19q` !="NA",]
data<-data[data$MGMT !="NA",]

data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
data$Age<-ifelse(data$Age>50,">50","<50")
data$Age<-factor(data$Age,levels = c("<50",">50"))
data[data$Gender=="female","Gender"]<-"Female"
data[data$Gender=="male",'Gender']<-"Male"
data$Gender<-factor(data$Gende,levels=c("Female","Male"))
data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
data[data$IDH=="WT","IDH"]<-"Wildtype"
data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
data[data$`Chr1p/19q`=="non-codel","Chr1p/19q"]<-"Non-codel"
data[data$`Chr1p/19q`=="codel","Chr1p/19q"]<-"Codel"
data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))

ezcox_res <- ezcox(data,
                   covariates = c("Risk score"),
                   controls = c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"),
                   return_models=TRUE)

tmp<-as.data.frame(ezcox_res[["res"]])
tmp$Variable<-c("Risk score","Gender","Age","Grade","","IDH","Chr1p/19q","MGMT")
tmp$Cohorts<-c(datasets,"","","","","","","")
multic_cox<-rbind(multic_cox,tmp)
multic_cox_model[[datasets]]<-get_models(ezcox_res)
show_models(multic_cox_model[[datasets]])

#CGGA
for (i in c("CGGA.325","CGGA.693","CGGA.1018")) {
  datasets<-i
  meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/CGGA_data//processed_data/cgga_clinic.rds")
  rs <- res[["riskscore"]][[optimal.model]][[datasets]]
  rownames(rs)<-rs$ID
  rownames(meta)<-meta$CGGA_ID
  meta<-meta[rownames(rs),]
  meta<-cbind(meta[,-7],rs)
  
  data<-dplyr::select(meta,c("Age","Gender","Grade","IDH_mutation_status","1p19q_codeletion_status","MGMTp_methylation_status",
                             "OS","OS.time","RS"))
  colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                    "status","time","RS")
  data$Age<-as.numeric(data$Age)
  data<-data[! is.na(data$Grade),]
  data<-data[! is.na(data$IDH),]
  data<-data[! is.na(data$MGMT),]
  data<-data[! is.na(data$`Chr1p/19q`),]
  
  data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
  data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
  data$Age<-ifelse(data$Age>50,">50","<50")
  data$Age<-factor(data$Age,levels = c("<50",">50"))
  data$Gender<-factor(data$Gende,levels=c("Female","Male"))
  data[data$Grade=="WHO II","Grade"]<-"G2"
  data[data$Grade=="WHO III","Grade"]<-"G3"
  data[data$Grade=="WHO IV","Grade"]<-"G4"
  data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
  data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
  data[data$MGMT=="un-methylated","MGMT"]<-"Unmethylated"
  data[data$MGMT=="methylated","MGMT"]<-"Methylated"
  data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
  data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))
  
  ezcox_res <- ezcox(data,
                     covariates = c("Risk score"),
                     controls = c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"),
                     return_models=TRUE)
  
  tmp<-as.data.frame(ezcox_res[["res"]])
  tmp$Variable<-c("Risk score","Gender","Age","Grade","","IDH","Chr1p/19q","MGMT")
  tmp$Cohorts<-c(datasets,"","","","","","","")
  multic_cox<-rbind(multic_cox,tmp)
  multic_cox_model[[datasets]]<-get_models(ezcox_res)
  show_models(multic_cox_model[[datasets]])
  
}

# gse16011
# datasets<-"GSE16011"
# meta <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/GSE16011_meta.rds")
# rs <- res[["riskscore"]][[optimal.model]][[datasets]]
# rownames(rs)<-rs$ID
# meta<-meta[rownames(rs),]
# meta<-cbind(meta,rs)
# 
# data<-dplyr::select(meta,c("Ageat diagnosis","Gender","grade","IDH1 (R132) mutation","1p","19q",
#                            "OS","OS.time","RS"))
# colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p","Chr19q",
#                   "status","time","RS")
# data$Age<-as.numeric(data$Age)
# data<-data[! is.na(data$Grade),]
# data<-data[! is.na(data$IDH),]
# data<-data[! is.na(data$Chr1p),]
# data<-data[! is.na(data$Chr19q),]
# 
# data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
# data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
# data$Age<-ifelse(data$Age>50,">50","<50")
# data$Age<-factor(data$Age,levels = c("<50",">50"))
# data$Gender<-factor(data$Gende,levels=c("Female","Male"))
# data[data$Grade=="gradeII","Grade"]<-"G2"
# data[data$Grade=="gradeIII","Grade"]<-"G3"
# data[data$Grade=="gradeIV","Grade"]<-"G4"
# data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
# data[data$IDH=="no mutation","IDH"]<-"Wildtype"
# data[data$IDH=="mutation","IDH"]<-"Mutant"
# data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
# data$`Chr1p/19q`<-ifelse(data$Chr1p %in% c("LOH","partial LOH") & data$Chr19q %in% c("LOH","partial LOH"),"Codel","Non-codel")
# data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))
# 
# ezcox_res <- ezcox(data,
#                    covariates = c("Risk score"),
#                    controls = c("Gender","Age","Grade","IDH","Chr1p/19q"),
#                    return_models=TRUE)
# show_models(get_models(ezcox_res))

save(multic_cox,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/multic_cox.Rdata")
save(multic_cox_model,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/multic_cox_model.Rdata")

## plot
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/multic_cox.Rdata")
dt<-multic_cox
dt$se <- (log(dt$upper_95) - log(dt$HR))/1.96
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$`HR (95% CI)` <- ifelse(is.na(dt$n_contrast), "",
                           sprintf("%.2f (%.2f - %.2f)",
                                   dt$HR, dt$lower_95, dt$upper_95))
dt$P <- ifelse(dt$p.value<0.001, "P<0.001",sprintf("%.3f",dt$p.value))
colnames(dt)[c(3:6,17)]<-c("Contrast","Reference","Number of contrast","Number of reference","P value")

tm <- forest_theme(core=list(bg_params=list(fill =c(rep("#3182BDFF",8),rep("#E6550DFF",8),rep("#31A354FF",8),rep("#756BB1FF",8)),
                                            alpha =rep(c(0.7,rep(0.5,7)),4))),
                   base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 16,
                   ci_col = "#762a83",
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 1,
                   footnote_fontface = "italic",
                   footnote_col = "red")

p <-forestploter::forest(dt[,c(13,1,4,3,6,5,15:17)],
                     est = dt$HR,
                     lower = dt$lower_95, 
                     upper = dt$upper_95,
                     sizes = dt$se,
                     #is_summary = c(rep(FALSE, nrow(dt)-2), TRUE,TRUE),
                     ci_column = 7,
                     ref_line = 1,
                     arrow_lab = c("Better", "Worse"),
                     #xlim = c(0, 1.5),
                     #ticks_at = c(0.5, 1, 2, 5,7.5),
                     x_trans="log2",
                     footnote = " Multivariate Cox Regression",
                     theme = tm)
p <- add_text(p, text = "Multivariate Cox regression in different cohorts",
              part = "header",
              row = 0,
              col = 4:7,
              just = c("center"),
              gp = gpar(fontface = "bold"))

p <- add_border(p, 
                part = "header", 
                row = c(0,1),
                gp = gpar(lwd = 1))


tiff("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/multi_rs_rsf_survivalSVM.tiff", units ="in",
     width = 10.5, height = 8.5,res = 600,compression ='zip')
p                
dev.off()

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
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

source("/export3/zhangw/Project_Cross/Project_Mime/Function/cal_RS_pre.prog.sig.R")
# rs.glioma = cal_RS_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data)
# rs.gbm = cal_RS_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data)
# rs.lgg = cal_RS_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data)
# rs.lgg.GBM = cal_RS_pre.prog.sig(type.sig = c('LGG','GBM'),list_input_data = list_train_vali_Data)
rs.glioma.lgg.gbm = cal_RS_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('LGG','GBM','Glioma'),list_input_data = list_train_vali_Data)
save(rs.glioma.lgg.gbm,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/rs.glioma.lgg.gbm.Rdata")

source("/export3/zhangw/Project_Cross/Project_Mime/Function/cal_cindex_pre.prog.sig.R")
# cc.gbm = cal_cindex_pre.prog.sig(type.sig = 'GBM',list_input_data = list_train_vali_Data)
# cc.lgg = cal_cindex_pre.prog.sig(type.sig = 'LGG',list_input_data = list_train_vali_Data)
# cc.glioma= cal_cindex_pre.prog.sig(type.sig = 'Glioma',list_input_data = list_train_vali_Data)
cc.glioma.lgg.gbm = cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data)
save(cc.glioma.lgg.gbm,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cc.glioma.lgg.gbm.Rdata")

source('/export3/zhangw/Project_Cross/Project_Mime/Function/cal_auc_pre.prog.sig.R')

auc.glioma.lgg.gbm.1 = cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data,AUC_time = 1,
                                            auc_cal_method = 'KM')
save(auc.glioma.lgg.gbm.1,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.1_km.Rdata")

auc.glioma.lgg.gbm.3 = cal_auc_pre.prog.sig(type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data,AUC_time = 3,
                                            auc_cal_method = 'KM')
save(auc.glioma.lgg.gbm.3,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.3_km.Rdata")

auc.glioma.lgg.gbm.5 = cal_auc_pre.prog.sig(type.sig = c('Glioma','LGG','GBM'),list_input_data = list_train_vali_Data,AUC_time = 5,
                                            auc_cal_method = 'KM')
save(auc.glioma.lgg.gbm.5,file="/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.5_km.Rdata")

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

load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

## PMID: 35145098 参考文献
# 热图展示所有 signature的预后关系, 单因素回归的HR和p值
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/rs.glioma.lgg.gbm.Rdata")
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/hr_comp.pdf',width = 36,height = 6,onefile = F)
HR_com(rs.glioma.lgg.gbm,
       res,
       model_name="RSF + survival-SVM",
       dataset=names(list_train_vali_Data),
       type = "categorical")
dev.off()

# 比较 Cindex
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cc.glioma.lgg.gbm.Rdata")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/cindex_comp.pdf',width = 32,height = 13,onefile = F)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="RSF + survival-SVM",
            dataset=names(list_train_vali_Data))
dev.off()

# 比较 auc(1,3,5)
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.1_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.3_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/auc.glioma.lgg.gbm.5_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.1y_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.3y_km.Rdata")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/all.auc.5y_km.Rdata")

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

####################### show model gene expression in p1 knockdown（刘宏伟 已经完成） ########################################################
### 欢迎补充
library(ComplexHeatmap)
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/1.Prog.Model/101ml.res.Rdata")
gene<- res[["ml.res"]][["RSF + survival-SVM"]][["var.names"]]

primaryGBMcell_expr <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/primaryGBMcell_expr.rds")
primaryGBMcell_meta <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/primaryGBMcell_meta.rds")

p1_primaryGBMcell_meta<-primaryGBMcell_meta[primaryGBMcell_meta$Treatment %in% c("shRNA1_PIEZO1_72h",
                                                                                 "shRNA2_PIEZO1_72h",
                                                                                 "shRNA_control_72h"),]
p1_primaryGBMcell_meta[p1_primaryGBMcell_meta=="shRNA_control_72h"]<-"Control"
p1_primaryGBMcell_meta[p1_primaryGBMcell_meta=="shRNA2_PIEZO1_72h"]<-"shRNA2"
p1_primaryGBMcell_meta[p1_primaryGBMcell_meta=="shRNA1_PIEZO1_72h"]<-"shRNA1"

p1_primaryGBMcell_expr <- primaryGBMcell_expr[,rownames(p1_primaryGBMcell_meta)]
p1_primaryGBMcell_expr <- p1_primaryGBMcell_expr[gene,]

annotation_colors = list(Treatment=c("Control"="#374E55","shRNA1"="#DF8F44","shRNA2"="#00A1D5"),
                         Cell_line=c("GBM001"="#B24745","GBM005"="#79AF97"))

meta<-p1_primaryGBMcell_meta[p1_primaryGBMcell_meta$Cell_line=="GBM001",]
data<-p1_primaryGBMcell_expr[,rownames(meta)]
p1<-pheatmap(data,scale = "row", cluster_cols = F,cluster_rows  = T,
         col=colorRampPalette(c(paletteer::paletteer_c("grDevices::Spectral", 30,direction = -1)))(30),name="Expression",
         show_colnames = F,border=F,border_color ="white",
         annotation_colors = annotation_colors,gaps_col = 9,
         annotation_names_col = F,
         annotation_col = dplyr::select(meta,c("Cell_line","Treatment")))

meta<-p1_primaryGBMcell_meta[p1_primaryGBMcell_meta$Cell_line=="GBM005",]
data<-p1_primaryGBMcell_expr[,rownames(meta)]
p2<-pheatmap(data,scale = "row", cluster_cols = F,cluster_rows  = T,
             col=colorRampPalette(c(paletteer::paletteer_c("grDevices::Spectral", 30,direction = -1)))(30),name="Expression",
             show_colnames = F,border=F,border_color ="white",
             annotation_colors = annotation_colors,gaps_col = 9,
             annotation_col = dplyr::select(meta,c("Cell_line","Treatment")))

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/1.Prog.Model/model_gene_p1_heatmap.pdf',
          width = 5.5,height = 3.5,onefile = F)
p1+p2
dev.off()


# 2.0 利用1.0中建立的optimal.model预后模型进行免疫浸润的分析（张益浩） -----------------------------------------------------

### 这里应该是同一套代码，在TCGA, CGGA325, CGGA693, CGGA1018 四个队列中都跑一遍 



 
# 3.0 利用1.0中建立的optimal.model预后模型进行富集的分析（张益浩） -----------------------------------------------------



# 4.0 分类变量的预测模型的建立--这里是预测免疫治疗反应与否(张炜计算,已经完成, 可视化刘宏伟)  -----------------------------------------------------

setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/4.ICI_response")

#######################利用ml建立预测免疫治疗反应的model (张炜 已经完成)############################################
# dir.create('4.ICI_response')
# source('/export3/zhangw/Project_Cross/Project_Mime/Proj/code/ICI.response.model.R')

#######################可视化建立的预测ICI反应的model (刘宏伟 已经完成)############################################
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/4.ICI_response/ICIresponse.model.all.p1.Rdata")
source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/4.ICI_response/auc_all.pdf',width = 6.5,height = 3.8,onefile = F)
auc_vis_category_all(res.p1.ici,dataset = c("training","validation"),
                     order= c("training","validation"))
dev.off()

methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
plot_list<-list()
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.p1.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/4.ICI_response/roc_all.pdf',width = 12,height = 12,onefile = F)
aplot::plot_list(gglist=plot_list,ncol=3)
dev.off()

# 5.0 核心feature的选择（和预后相关）(计算张炜，已完成， 可视化刘宏伟)  -----------------------------------------------------

#######################核心feature的选择（和预后相关) (张炜 已完成)############################################

setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res")
# dir.create('5.CoreFeature')
setwd("/export3/zhangw/Project_Cross/Project_Mime/Proj/res/5.CoreFeature")

# 计算
# source('/export3/zhangw/Project_Cross/Project_Mime/Proj/code/ML.CoreFeature.R')

#######################可视化 (刘宏伟 已经完成)################################################################

source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")
load('/export3/zhangw/Project_Cross/Project_Mime/Proj/res/5.CoreFeature/feature.all.res.Rdata')

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/5.CoreFeature/core_feature_intersect_tcga.pdf',
          width = 15,height = 7.5,onefile = F)
core_feature_select(res.feature.all)
dev.off()

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/5.CoreFeature/core_feature_intersect_tcga_rank.pdf',
          width = 4.8,height = 7,onefile = F)
core_feature_rank(res.feature.all, top=50)
dev.off()

# 6.0 在核心feature中选择某一个关键基因进行单基因分析的一整套(预后,差异,回归,龙哥的paper那一套)（刘宏伟） ----------------------------------------------------
gene<-"MCAM"
## expression iin piezo1 knockdown cell line
primaryGBMcell_expr <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/primaryGBMcell_expr.rds")
primaryGBMcell_meta <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/primaryGBMcell_meta.rds")

p1_primaryGBMcell_meta<-primaryGBMcell_meta[primaryGBMcell_meta$Treatment %in% c("shRNA1_PIEZO1_72h",
                                                                                 "shRNA2_PIEZO1_72h",
                                                                                 "shRNA_control_72h"),]
p1_primaryGBMcell_meta[p1_primaryGBMcell_meta=="shRNA_control_72h"]<-"Control"
p1_primaryGBMcell_meta[p1_primaryGBMcell_meta=="shRNA2_PIEZO1_72h"]<-"shRNA2"
p1_primaryGBMcell_meta[p1_primaryGBMcell_meta=="shRNA1_PIEZO1_72h"]<-"shRNA1"

p1_primaryGBMcell_expr <- primaryGBMcell_expr[,rownames(p1_primaryGBMcell_meta)]
p1_primaryGBMcell_meta$gene<-t(p1_primaryGBMcell_expr[gene,])[,1]

library(ggplot2)
library(ggpubr)
comparisons<-list(c("Control","shRNA1"),c("Control","shRNA2"))

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/5.CoreFeature/mcam_primarycell.pdf',
          width = 3.5,height = 5,onefile = F)
ggplot(p1_primaryGBMcell_meta, aes(Treatment, gene,fill =Treatment ))+
  scale_fill_manual(values =c("#374E55","#DF8F44","#00A1D5"))+
  geom_violin(alpha=0.25, size=1,color="NA",trim=FALSE)+
  geom_boxplot( outlier.size = -1, color="black",lwd=0.2, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 0.6) +
  facet_wrap(~Cell_line)+
  theme_bw()+
  labs(title = gene)+
  ylab("Expression") +
  xlab("") +
  stat_compare_means(comparisons = comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
dev.off()

## survival in cohorts
source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")
load("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/data/Glioma.cohort.Rdata")

survplot <- vector("list",11) 
for (i in c(1:11)) {
  print(survplot[[i]]<-core_feature_sur(gene, 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/5.CoreFeature/gene_km.pdf',width = 10,height = 28,onefile = F)
aplot::plot_list(gglist=survplot,ncol=2)
dev.off()

## correlation in cohorts
source("/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/code/plot_function.R")

dataset_col<-c("#3182BDFF","#E6550DFF","#31A354FF","#756BB1FF","#636363FF","#6BAED6FF","#FD8D3CFF","#74C476FF",
                          "#9E9AC8FF","#969696FF","#9ECAE1FF","#FDAE6BFF","#A1D99BFF","#BCBDDCFF","#BDBDBDFF","#C6DBEFFF",
                          "#FDD0A2FF","#C7E9C0FF","#DADAEBFF","#D9D9D9FF")

corplot <- list()
for (i in c(1:4,6:11)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                                       color = dataset_col[i],
                                       feature1="PIEZO1",
                                       feature2="MCAM",
                                       method="pearson"))
}
corplot<-corplot[-5]

cairo_pdf('/export/bioinfo-team/home/liuhw/bioinfo_mill/Mime_proj/res/5.CoreFeature/gene_cor.pdf',width = 10,height = 25,onefile = F)
aplot::plot_list(gglist=corplot,ncol=2)
dev.off()
