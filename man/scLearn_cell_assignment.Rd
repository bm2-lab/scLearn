\name{scLearn_cell_assignment}
\alias{scLearn_cell_assignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To predict the cell assignment result with trained scLearn model.
}
\description{
To predict the cell assignment result with trained scLearn model.
}
\usage{
scLearn_cell_assignment(scLearn_model_learning_result, expression_profile_query, vote_rate = 0.6, diff = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scLearn_model_learning_result}{
  The result calculated by function "scLearn_model_learning".
}
  \item{expression_profile_query}{
  A dataframe showing the query expression profile. The row is gene and The column is sample.
}
  \item{vote_rate}{
  For 10 times bootstrappings, the percentage of 10 results that has the same assigned cell type. Default is 0.6. If you want a less strict assignment result for "unassigned", you can set it smaller, such as 0.5.
}
  \item{diff}{
  For similar cell types, it is hard to distinguish. If a query cell, if its similarity value to different cell types is too close, then it will be assigned as "unassigned". Default is 0.05. If you want a less strict assignment result for "unassigned", you can set it smaller,such as 0.01.
}
}
\value{
A dataframe with two columns: "predicted cell type" and "sample name".
}
\author{
Bin Duan
}
