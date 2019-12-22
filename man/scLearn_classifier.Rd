\name{scLearn_classifier}
\alias{scLearn_classifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training the classifier with scLearn method.
}
\description{
Training the classifier with scLearn method.
}
\usage{
scLearn_classifier(high_varGene_names, expression_profile, sample_information, bootstrap_times = 10, cutoff = 0.01)
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
  A character vector showing the cell type of each sample. The vector's column name is the sample name.
}
  \item{bootstrap_times}{
  The times for bootstrapping. Default is 10.
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
