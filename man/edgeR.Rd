\name{edgeR}
\alias{edgeR}
\title{handler for performing edgeR on subsampled RNA-Seq data}
\usage{
edgeR(count.matrix, treatment, pair = NULL)
}
\arguments{
  \item{count.matrix}{Matrix of unnormalized counts}

  \item{treatment}{vector describing the treatment of each
  sample (column)}

  \item{pair}{Pair of treatments to compare. If there are
  only two treatments, this is optional.}
}
\value{
A data frame with one row for each gene, with the columns

\item{coefficient}{log2 fold change between the conditions}
\item{pvalue}{P-value using a negative binomial test
comparing the conditions}
}
\description{
This encapsulates the edgeR \code{exactTest} function. It
takes a count matrix, potentially subsampled, a vector
describing the treatment of each sample, and the pair of
conditions to compare, and performs a negative binomial
exact test for differential expression between each pair of
conditions.
}

