\name{DrawCluster}
\alias{DrawCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drawing the cluster graph with "tsne" or "umap".
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DrawCluster(data, label = NULL, point_size = 1, method = c("tsne", "umap"), draw_cluster_text = TRUE, calculated = TRUE, pca = TRUE, perplexity = 100, plot = TRUE, seed = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A matrix with column is sample, row is feature.
}
  \item{label}{
  A vector showing the sample label corresponding the colnames of the matrix "data".
}
  \item{point_size}{
  The size of point. Default is 1.
}
  \item{method}{
  The method to plot which should be one of "tsne" and "umap". Default is "tsne".
}
  \item{draw_cluster_text}{
  TRUE by default. If TRUE, the cluster label will be drawed in the corresponding place in the graph.
}
  \item{calculated}{
  TRUE by default. If FALSE, "tsne" or "umap" won't be performed, the graph will be draw with data offered directly.
}
  \item{pca}{
  TRUE by default. If FALSE, pca won't be performed.
}
  \item{perplexity}{
  Default is 100.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{seed}{
  Default is 1.
}
}
\value{
\item{p }{The ggplot object for plotting.}
\item{x }{A matrix showing coordinates calculated by "tsne" or "umap".}
\item{cell_group}{A vector showing the cell label.}
%% ...
}
\author{
Bin Duan
}
