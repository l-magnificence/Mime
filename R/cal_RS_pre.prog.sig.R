#' Title
#'
#' @param use_your_own_collected_sig
#' @param collected_sig_table
#' @param type.sig
#' @param list_input_data
#'
#' @return
#' @export
#'
#' @examples
cal_RS_pre.prog.sig <- function(use_your_own_collected_sig, # 是否使用您自己收集的signature， T or F
                                collected_sig_table, # 列名分别为
                                # "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol"
                                # 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'
                                type.sig, ### prognostic signature 的类型，c('Glioma','LGG','GBM')， c('Glioma'),c('Glioma','LGG')
                                list_input_data # list of the cohorts(ID,OS.time, OS····)经过了log2（x+1）转化
) {
  library(tidyverse)
  library(survival)

  if (use_your_own_collected_sig) {
    sig.input <- collected_sig_table
  } else {
    pre.prog.sig <- MIME::pre.prog.sig

    if (all(type.sig %in% names(pre.prog.sig))) {
      if (length(type.sig) == 1) {
        sig.input <- pre.prog.sig[[type.sig[1]]]
      } else {
        sig.input <- pre.prog.sig[[type.sig[1]]]

        for (i in 2:length(type.sig)) {
          sig.input <- rbind(sig.input, pre.prog.sig[[type.sig[i]]])
        }
      }
    } else {
      stop("please provide correct type.sig")
    }
  }


  sig.input$Coef <- as.numeric(sig.input$Coef)
  sig.input$symbol <- gsub("-", ".", sig.input$symbol)



  # Replace '-' in column names with '.'
  list_input_data <- lapply(list_input_data, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    return(x)
  })

  common_feature <- c("ID", "OS.time", "OS", unique(sig.input$symbol))

  # for (i in names(list_input_data)) {
  #   common_feature = intersect(common_feature, colnames(list_input_data[[i]]))
  # }

  returnIDtoRS <- function(rs.table.list, rawtableID) {
    for (i in names(rs.table.list)) {
      rs.table.list[[i]]$ID <- rawtableID[[i]]$ID
      rs.table.list[[i]] <- rs.table.list[[i]] %>% dplyr::select("ID", everything())
    }

    return(rs.table.list)
  }

  list_input_data <- lapply(list_input_data, function(x) {
    not.gene <- common_feature[which(!common_feature %in% colnames(x))]
    cons.mat <- as.data.frame(matrix(rep(0, length(not.gene) * nrow(x)), nrow = nrow(x), ncol = length(not.gene)))
    colnames(cons.mat) <- not.gene


    x <- cbind(x, cons.mat)

    x <- x[, common_feature]
    return(x)
  })




  list_input_data <- lapply(list_input_data, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, as.numeric)
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x[, c(2:3)] <- apply(x[, c(2:3)], 2, as.numeric)
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x <- x[!is.na(x$OS.time) & !is.na(x$OS), ]
    return(x)
  })

  list_input_data <- lapply(list_input_data, function(x) {
    x <- x[x$OS.time > 0, ]
    return(x)
  })
  # use the mean replace the NA
  list_input_data <- lapply(list_input_data, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })


    return(x)
  })




  model.name <- unique(sig.input$model)

  val_dd_list <- lapply(list_input_data, function(x) {
    x[, c("OS.time", "OS", common_feature)]
  })


  rs.table <- lapply(model.name, function(z) {
    coef.tab <- sig.input[sig.input$model == z, c("Coef", "symbol")]

    val_dd_list2 <- lapply(val_dd_list, function(x) {
      x[, c("OS.time", "OS", coef.tab$symbol)]
    })


    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2],
        RS = apply(as.data.frame(x[, -c(1:2)]), 1, function(x) {
          x %*% coef.tab$Coef
        })
      )
    })
    rs <- returnIDtoRS(rs.table.list = rs, rawtableID = list_input_data)
    return(rs)
  })
  names(rs.table) <- model.name


  return(rs.table)
}