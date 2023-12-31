% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TME_deconvolution_all.R
\name{TME_deconvolution_all}
\alias{TME_deconvolution_all}
\title{Tumor microenvironment deconvolution using immunedeconv}
\usage{
TME_deconvolution_all(
  inputmatrix.list,
  deconvolution_method = c("xcell", "epic", "abis", "estimate", "cibersort",
    "cibersort_abs"),
  microarray_names = "none",
  indications = NULL,
  tumor = TRUE,
  column = "gene_symbol",
  rmgenes = NULL,
  scale_mrna = TRUE,
  expected_cell_types = NULL,
  ...
)
}
\arguments{
\item{inputmatrix.list}{A gene expression dataframe. The first three of the column names are, in order, ID,OS.time, OS. Columns starting with the fourth are gene symbols. OS.time is a numeric variable in days. OS is a numeric variable containing 0, 1. 0: Alive, 1: Dead.}

\item{deconvolution_method}{Deconvolution Methods in immunedeconv::deconvolution_methods}

\item{microarray_names}{Please tell us which datasets are microarray, use the names of elements in inputmatrix.list. such as c("CGGA.array", "GSE108474", "GSE16011", "GSE43289", "GSE7696") if none, enter "none".}

\item{indications}{a character vector with one indication per
sample for TIMER and ConsensusTME(immunedeconv::timer_available_cancers).
Argument is ignored for all other methods.}

\item{tumor}{use a signature matrix/procedure optimized for tumor samples,
if supported by the method. Currently affects EPIC and}

\item{column}{Only relevant if `gene_expression` is an ExpressionSet. Defines in which column
of fData the HGNC symbol can be found.}

\item{rmgenes}{a character vector of gene symbols. Exclude these genes from the analysis.
Use this to exclude e.g. noisy genes.}

\item{scale_mrna}{logical. If FALSE, disable correction for mRNA content of different cell types.
This is supported by methods that compute an absolute score (EPIC)}

\item{expected_cell_types}{Limit the analysis to the cell types given in this list. If the cell
types present in the sample are known *a priori*, setting this can improve results for
xCell (see https://github.com/grst/immunedeconv/issues/1).}

\item{...}{arguments passed to the respective method}
}
\value{
`list` containing `data.frame` with `cell_type` as first column and a column with the
    calculated cell fractions for each sample.

a list containing deconvolution scores in each cohorts and ML methods
}
\description{
Tumor microenvironment deconvolution using immunedeconv
}
\examples{
test.devo <- TME_deconvolution_all(list_train_vali_Data)
}
