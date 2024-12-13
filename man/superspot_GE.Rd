% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCimplify_SpatialDLS.R
\name{superspot_GE}
\alias{superspot_GE}
\title{Simplification of spatial transcriptomic dataset}
\usage{
superspot_GE(
  MC,
  ge,
  groups,
  mode = c("average", "sum"),
  weights = NULL,
  do.median.norm = FALSE
)
}
\arguments{
\item{MC}{Metaspot object obtained from SCimplify_SpatialDLS function}

\item{ge}{gene expression matrix (or any coordinate matrix) with genes as rows and cells as cols}

\item{groups}{vector of membership (assignment of single-cell to metacells)}

\item{mode}{string indicating whether to average or sum up \code{ge} within metacells}

\item{weights}{vector of a cell weight (NULL by default), used for computing average gene expression withing cluster of metaells}

\item{do.median.norm}{whether to normalize by median value (FALSE by default)}
}
\value{
a matrix of simplified (averaged withing groups) data with ncol equal to number of groups and nrows as in the initial dataset
}
\description{
This function converts (i.e., averages or sums up) gene-expression matrix of spatial transcriptomic data into a gene expression
matrix of metacells
}
