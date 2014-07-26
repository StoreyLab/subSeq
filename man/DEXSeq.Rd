\name{DEXSeq}
\alias{DEXSeq}
\title{handler for performing DEXSeq on subsampled RNA-Seq data}
\usage{
DEXSeq(count.matrix, design, geneIDs, exonIDs)
}
\arguments{
  \item{count.matrix}{Matrix of unnormalized counts. The
  row names should represent the exon names}

  \item{design}{Study design data frame for DEXSeq}

  \item{geneIDs}{vector of gene IDs- one for each row of
  the count matrix (repeated for each exon in a gene)}

  \item{exonIDs}{vector of exon IDs- one for each row of
  the count matrix}
}
\value{
A data frame with one row for each gene. This contains all
columns from the DEUresultTable function in DEXSeq.
}
\description{
This encapsulates the DEXSeq dispersion estimation and
\code{testForDEU} functions. It takes a count matrix,
potentially subsampled, and a data table representing the
study design, and tests for differential exon usage.
}
\examples{
data( "pasillaExons", package="pasilla" )
count.matrix = pasillaExons
}

