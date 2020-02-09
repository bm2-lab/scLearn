\name{Cell_type_filter}
\alias{Cell_type_filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering out cell types and/or cell time point.
}
\description{
Filtering out cell types and/or cell time point without the minimal amount of cells.
}
\usage{
Cell_type_filter(expression_profile, sample_information_cellType,sample_information_timePoint=NULL, min_cell_number = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile. The row is gene and The column is sample.
}
  \item{sample_information_cellType}{
  A character vector showing the cell type of each sample. The column name of the vector is the sample name.
}
  \item{sample_information_timePoint}{
  A character vector showing the time point of each sample. The column name of the vector is the sample name. The default is NULL.
}
  \item{min_cell_number}{
  The minimal cell number to filter cell type. Default is 10.
}
}
\value{
\item{expression_profile}{A dataframe showing expression profile after quality control.}
\item{sample_information_cellType}{A character vector showing the cell type of each sample after quality control.}
\item{sample_information_timePoint}{A character vector showing the time point of each sample after quality control.}

}
\author{
Bin Duan
}
