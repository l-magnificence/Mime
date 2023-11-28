#' Compare immune cell infiltration between high-risk and low-risk group determined by specific model
#'
#' Creates a heatmap of immune cell infiltration among different methods 
#'
#' @param object Output of function ML.Dev.Prog.Sig 
#' @param object2 Output of function TME_deconvolution_all
#' @param model_name Model name in object used to compare
#' @param col If NULL, color values are set to the default colors. Otherwise, you can specify two color values for high and low risk group
#' @param dataset Names for specific dataset 
#'
#' @return a heatmap object
#' @export
#'
#' @examples
#' immuno_heatmap(res,
#'             devo,
#'             model_name="StepCox[backward] + plsRcox",
#'             dataset="Dataset1")
#' 
#' 
immuno_heatmap<-function(object,# output of ML.Dev.Prog.Sig
                         object2, # output of TME_deconvolution_all
                         model_name, ## input specific model name
                         col = NULL, # two color value for high risk and low risk  group
                         dataset # input datasets name
                            
){
  
  library(gplots)
  library(viridis)
  library(tidyverse)
  library(dplyr)
  library(ComplexHeatmap)
  
  if (is.null(col) == T) {
    col <- c("#868686", "#B24745") ## default color value
  } else {
    col <- col
  }
  
  risk <- object[["riskscore"]][[model_name]][[dataset]]
  hmdat <- object2[[dataset]][["tme_combine"]]
  
  risk$risk <- ifelse(risk$RS > median(risk$RS),"high","low")
  hmdat <- lapply(hmdat, function(tibble_i) {
    tibble_i <-  tibble_i |> as.data.frame() |> tibble::column_to_rownames("cell_type")
    return(tibble_i)
  })
  
  ht_opt$message = FALSE
  
  samorder <- risk[order(risk$RS),]$ID
  annCol <- data.frame(RiskScore = scale(risk$RS),
                       RiskType = risk$risk,
                       row.names = risk$ID,
                       stringsAsFactors = F)
  
  
  annColors <- list("RiskScore" = viridis::rocket(64),
                    "RiskType" = c("high" = col[2],"low" = col[1]))
  
  color <- list(greenred(64),bluered(64),turbo(64),inferno(64),viridis(64),magma(64),plasma(64),mako(64),cividis(64),rocket(64))
  
  if (length(hmdat) > 1) {
    plot_list<-list()
    plot_list[[1]]<-ComplexHeatmap::pheatmap(mat = t(as.matrix(hmdat[[1]][,samorder])),
                                             border_color = NA,
                                             color = color[[1]],
                                             cluster_rows = F,
                                             cluster_cols = F,
                                             show_rownames = F,
                                             show_colnames = T,
                                             annotation_row = annCol[samorder,,drop = F],
                                             annotation_colors = annColors,
                                             cellwidth = 10,
                                             cellheight = 0.8,
                                             gaps_row = table(annCol$RiskType)[2],
                                             name = names(hmdat)[1])
    for (i in 2:length(hmdat)) {
      plot_list[[i]]<-ComplexHeatmap::pheatmap(mat = t(as.matrix(hmdat[[i]][,samorder])),
                                               border_color = NA,
                                               color = color[[i]],
                                               cluster_rows = F,
                                               cluster_cols = F,
                                               show_rownames = F,
                                               show_colnames = T,
                                               # annotation_row = annCol[samorder,,drop = F],
                                               # annotation_colors = annColors,
                                               cellwidth = 10,
                                               cellheight = 0.8,
                                               gaps_row = table(annCol$RiskType)[2],
                                               name = names(hmdat)[i])
    }
    
    ht <- plot_list[[1]]
    for (i in 2:length(plot_list)) {
      ht <- ht + plot_list[[i]] 
    }
    
  } else {
    plot_list<-list()
    for (i in 1:length(hmdat)) {
      plot_list[[i]]<-ComplexHeatmap::pheatmap(mat = t(as.matrix(hmdat[[i]][,samorder])),
                                               border_color = NA,
                                               color = color[[i]],
                                               cluster_rows = F,
                                               cluster_cols = F,
                                               show_rownames = F,
                                               show_colnames = T,
                                               annotation_row = annCol[samorder,,drop = F],
                                               annotation_colors = annColors,
                                               cellwidth = 10,
                                               cellheight = 0.8,
                                               gaps_row = table(annCol$RiskType)[2],
                                               name = names(hmdat)[i])
    }
    ht <- plot_list[[1]]
  }
  
  p1<-draw(ht,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
  print(p1)
  
}


