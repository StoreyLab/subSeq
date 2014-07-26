\name{DESeq2}
\alias{DESeq2}
\title{handler for performing DESeq2 on subsampled RNA-Seq data}
\usage{
DESeq2(count.matrix, treatment)
}
\arguments{
  \item{count.matrix}{Matrix of unnormalized counts}

  \item{treatment}{vector describing the treatment of each
  sample (column)}
}
\value{
A data frame with one row for each gene, with the columns

\item{coefficient}{fold change between the conditions}
\item{pvalue}{P-value using a negative binomial test
comparing the conditions}
}
\description{
This encapsulates the DESeq2 \code{nbinomWaldTest}
function. It takes a count matrix, potentially subsampled,
a vector describing the treatment of each sample, and the
pair of conditions to compare, and performs a negative
binomial Wald test for differential expression between each
pair of conditions.
}

