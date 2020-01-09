\name{Cell_qc}
\alias{Cell_qc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quality control for single cell transcriptome.
}
\description{
scLearn evaluates the cell quality by three factors, i.e. the number of genes detected (default 500 to 10000), the number of unique molecular identifier induced (default above 1500), and the percentage of mitochondrial genes detected (default below 10\% among all the detected genes). Only cells satified the whole three factors are retained.
}
\usage{
Cell_qc(expression_profile, sample_information = NULL, species = "Hs", gene_low = 500, gene_high = 10000, mito_high = 0.1, umi_low = 1500, umi_high = Inf, logNormalize = TRUE, plot = FALSE, plot_path = "./quality_control.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile. The row is gene and The column is sample.
}
  \item{sample_information}{
  A character vector showing the cell type of each sample. The column name of the vector is the sample name.
}
  \item{species}{
  The species the cells belong to. Currently species="Hs" for homo sapiens or species="Mm" for mus musculus are available. It is used to detect mitochondrial genes, so species doesn't matter if data won't be considered the percentage of mitochondrial genes.
}
  \item{gene_low}{
  The minimun gene number for cells. Cells with genes below this threshold are filtered(default 500).
}
  \item{gene_high}{
  The maximum gene number for cells. cells with genes above this threshold are filtered(default 10000).
}
  \item{mito_high}{
  The maximum percentage of mitochondrial genes detected(default 10\% among all the detected genes).
}
  \item{umi_low}{
  The minimum number of unique molecular identifier induced (default 1500).
}
  \item{umi_high}{
  The maximun number of unique molecular identifier induced (default Inf)
}
  \item{logNormalize}{
  TRUE by default. If FALSE, the data will not be performed log normalized.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{plot_path}{
    The path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{expression_profile}{A dataframe showing expression profile after quality control.}
\item{sample_information}{A character vector showing the cell type of each sample after quality contro.}
}
\author{
Bin Duan
}

