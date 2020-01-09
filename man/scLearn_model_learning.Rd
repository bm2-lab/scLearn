\name{scLearn_model_learning}
\alias{scLearn_model_learning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training the learning-based model with scLearn method.
}
\description{
Training the learning-based model with scLearn method.
}
\usage{
scLearn_model_learning(high_varGene_names, expression_profile, sample_information, bootstrap_times = 10, cutoff = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{high_varGene_names}{
  A vecter showing the selected genes.
}
  \item{expression_profile}{
  A dataframe showing the reference expression profile. The row is gene and The column is sample.
}
  \item{sample_information}{
  A character vector showing the cell type of each sample. The column name of the vector is the sample name.
}
  \item{bootstrap_times}{
  The times for bootstrapping which should be at least larger than 1. Default is 10.
}
  \item{cutoff}{
  The cutoff for selecting similarity threshold for each cell type. Default is 0.01.
}
}
\value{
\item{high_varGene_names }{A vecter showing the selected genes.}
\item{simi_threshold_learned }{A list showing the similarity threshold for each bootstrapping.}
\item{feature_matrix_learned}{A list showing the learned feature for each bootstrapping}
\item{trans_matrix_learned}{A list showing the learned metric for each bootstrapping}
}
\author{
Bin Duan
}
