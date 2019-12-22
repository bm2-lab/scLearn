\name{sankey_plot}
\alias{sankey_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the sankey graph.
}
\description{
Plotting the sankey graph.
}
\usage{
sankey_plot(predict_result, sample_information_reference, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predict_result}{
  A data frame showing the predicted result from function "scLearn_predict".
}
  \item{sample_information_reference}{
  A character vector showing the cell type of each reference sample. The vector's column name is the sample name.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
}
\value{
A data frame showing the predicted result after a little modification.
}
\author{
Bin Duan
}
