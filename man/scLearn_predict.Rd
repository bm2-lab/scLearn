\name{scLearn_predict}
\alias{scLearn_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To predict the cell assignment result with trained classifier.
}
\description{
To predict the cell assignment result with trained classifier.
}
\usage{
scLearn_predict(scLearn_classifier_result, expression_profile_query, vote_rate = 0.6, diff = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scLearn_classifier_result}{
  The classifier trained before from function "scLearn_classifier".
}
  \item{expression_profile_query}{
  A dataframe showing the query expression profile. The row is gene and The column is sample.
}
  \item{vote_rate}{
  Default is 0.6.
}
  \item{diff}{
  Default is 0.05.
}
}
\value{
A dataframe with two columns: "predicted cell type" and "sample name".
}
\author{
Bin Duan
}
