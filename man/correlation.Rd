\name{correlation}
\alias{correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating the correlation of each sample of a matrix.
}
\description{
Sometimes, we need calculate the correlation of each sample of a matrix. This function can calculate the correlation among each sample of a matrix with metrics(pearson,spearman,cosin or euclidean) by parallel.
}
\usage{
correlation(matrix, method = c("pearson", "spearman", "cosin", "euclidean"), cpu_num = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{
  A matrix with column is sample, row is feature.
}
  \item{method}{
  The metric to calculate the correlation. Four metrics("pearson", "spearman", "cosin", "euclidean") are supported.
}
  \item{cpu_num}{
  The number of cpu used for parallel calculation.
}
}
\value{
A matrix showing the correlation with column and row are both sample.
}
\author{
Bin Duan
}
