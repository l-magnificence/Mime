% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_RS_ml_res.R
\name{cal_RS_ml_res}
\alias{cal_RS_ml_res}
\title{Calculate risk scores of Machine Learning Models in all data}
\usage{
cal_RS_ml_res(
  res.by.ML.Dev.Prog.Sig = NULL,
  train_data,
  inputmatrix.list,
  mode = NULL,
  single_ml = NULL,
  double_ml1 = NULL,
  double_ml2 = NULL,
  alpha_for_Enet = NULL,
  direction_for_stepcox = NULL
)
}
\arguments{
\item{res.by.ML.Dev.Prog.Sig}{The results of function ML.Dev.Prog.Sig}

\item{train_data}{The training data using in ML.Dev.Prog.Sig}

\item{inputmatrix.list}{A list contain the data frames (colnames:ID,OS.time,OS,other variables), log2(x+1)， OS.time(day), OS(0/1)}

\item{mode}{Choose MF models: 'all', 'single', 'double'}

\item{single_ml}{If the mode is set to "single", you must fill in the following models: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso").}

\item{double_ml1}{If the mode is set to "double", you need to fill in the modeling methods here: c('RSF', "StepCox","CoxBoost","Lasso").}

\item{double_ml2}{If the mode is set to "double", you need to fill in the modeling methods here: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")}

\item{alpha_for_Enet}{One of the values from 0.1 to 0.9. c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9). There are some conditions you could not set this parameter. [1] The mode is 'all'. [2] The mode is 'single' or 'double', but the 'Enet' is not included in the algorithms you choose.}

\item{direction_for_stepcox}{The parameter for the StepCox. One  from "both", "backward", "forward". There are some conditions you could not set this parameter. [1] The mode is 'all'. [2] The mode is 'single' or 'double', but the 'StepCox' is not included in the algorithms you choose.}
}
\value{
A list of the risk score calculated by the predictive model in each data in the data of the input matrix list.
}
\description{
Calculate risk scores of Machine Learning Models in all data
}
