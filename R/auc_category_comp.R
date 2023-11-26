#' Compare the AUC of our model with previously published models for predicting ICI response
#'
#' Creates a distribution plot of AUC among different mdoels 
#'
#' @param object Output of function ML.Dev.Pred.Category.Sig
#' @param object2 Output of function cal_auc_previous_sig
#' @param model_name Method name in object used to compare
#' @param dataset_col If NULL, color values are set to the default colors. Otherwise, you can specify consistent number of color value for datasets
#' @param dataset A vector of names for all datasets 
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' auc_category_comp(res.ici,
#'             auc.other.pre,
#'             model_name="svmRadialWeights",
#'             dataset=names(list_train_vali_Data))
#' 
#' 
auc_category_comp<-function(object,# output of ML.Dev.Pred.Category.Sig
                      object2, # output of cal_auc_previous_sig
                      model_name, ## c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
                      dataset_col=NULL, # color value for cohort
                      dataset # input datasets name
                      
){
  
  library(ggplot2)
  library(aplot)
  
  if (is.null(dataset_col) == T){
    dataset_col<-c("#3182BDFF","#E6550DFF","#31A354FF","#756BB1FF","#636363FF","#6BAED6FF","#FD8D3CFF","#74C476FF",
                              "#9E9AC8FF","#969696FF","#9ECAE1FF","#FDAE6BFF","#A1D99BFF","#BCBDDCFF","#BDBDBDFF","#C6DBEFFF",
                              "#FDD0A2FF","#C7E9C0FF","#DADAEBFF","#D9D9D9FF") ## default 20 color values
  } else {
    dataset_col<-dataset_col
  }
  
  method <- model_name
  
  tmp<-data.frame()
  for (i in 1:length(dataset)) {
    auc<-object[["auc"]][[dataset[i]]]
    auc$method<-rownames(auc)
    auc$dataset<-dataset[i]
    rownames(auc)<-NULL
    tmp<-rbind(tmp,auc)
  }
  tmp<-tmp[tmp$method==method,]
  tmp<-tmp[,c(1,5,4)]
  colnames(tmp)<-c("AUC","Dataset","Model")
  tmp$text<-"red"
    
  tmp2<-data.frame()
  for (i in 1:length(object2)) {
    auc<-as.data.frame(object2[[i]])
    auc$Dataset<-rownames(auc)
    auc$Model<-names(object2)[i]
    rownames(auc)<-NULL
    tmp2<-rbind(tmp2,auc)
  }
  tmp2$text<-"black"
    
  auc_comp_other<-rbind(tmp,tmp2)
  
  plot_list<-list()
  for (t in 1:length(dataset)) {
    auc_comp_other_select<-auc_comp_other[auc_comp_other$Dataset==dataset[t],]
    auc_comp_other_select<-auc_comp_other_select[order(auc_comp_other_select$AUC,decreasing = F),]
    auc_comp_other_select$Model<-factor(auc_comp_other_select$Model,levels = unique(auc_comp_other_select$Model))
    
    plot_list[[dataset[t]]]<- 
      ggplot(auc_comp_other_select, aes(x=Model, y=AUC)) +
      geom_segment( aes(x=Model, xend=Model, y=0, yend=AUC,color=Dataset)) +
      geom_point(aes(color=Dataset))+
      scale_color_manual(values=dataset_col[t],name="Cohort")+
      theme(panel.grid=element_blank(),
            axis.text.y.left  = (element_text(color=auc_comp_other_select$text)),
            panel.border = element_rect(colour = "black", fill=NA,size = 0.3),
            legend.position = "",
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill="white"))+
      #scale_y_continuous(position = "right")+
      labs(y="AUC",x="",title = "")+
      coord_flip()
  }
  
  p1<-aplot::plot_list(gglist=plot_list,ncol=length(dataset))
  print(p1)

}


