\name{edgeR.glm}
\alias{edgeR.glm}
\title{handler for performing an edgeR glm on subsampled RNA-Seq data
This encapsulates the edgeR \code{glmFit} and \code{glmLRT}
functions. It takes a count matrix, potentially subsampled, a vector describing
the group (usually treatment) of each sample, and a model matrix, and performs
a negative binomial regression to find each coefficient, followed by a
likelihood ratio test to test for significance of each}
\usage{
edgeR.glm(count.matrix, mm, column = 2, group = NULL)
}
\arguments{
  \item{count.matrix}{Matrix of unnormalized counts}

  \item{group}{vector describing the group of each sample}

  \item{mm}{model matrix describing one or more regressors
  to fit}

  \item{column}{on which column of the model matrix should
  the likelihood ratio test should be performed}
}
\value{
A data frame with one row for each gene, with the columns

\item{effect.size}{log fold change between the conditions}
\item{pvalue}{P-value using a likelihood ratio test
comparing the conditions}
}
\description{
handler for performing an edgeR glm on subsampled RNA-Seq
data This encapsulates the edgeR \code{glmFit} and
\code{glmLRT} functions. It takes a count matrix,
potentially subsampled, a vector describing the group
(usually treatment) of each sample, and a model matrix, and
performs a negative binomial regression to find each
coefficient, followed by a likelihood ratio test to test
for significance of each
}

