\name{Feature_selection_M3Drop}
\alias{Feature_selection_M3Drop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selecting important features for the following analysis with M3Drop package.
}
\description{
M3Drop is a method to select genes with high dropout rate. Here, we use this method to select features for our following analysis.
}
\usage{
Feature_selection_M3Drop(expression_profile, log_normalized = TRUE, threshold = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile. The row is gene and The column is sample.
}
  \item{log_normalized}{
  TRUE by default. If TRUE, it means the data has been logNormalized before.
}
  \item{threshold}{
  The threshold to select significant genes. Default is 0.05.
}
}
\value{
A vecter showing the selected gene names.
}
\author{
Bin Duan
}
