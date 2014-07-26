\name{getSeed}
\alias{getSeed}
\title{Extract the global random seed from a subsamples object}
\usage{
getSeed(ss)
}
\arguments{
  \item{ss}{A subsamples object, returned from the
  \code{subsample} function, or a summary of that object}
}
\description{
A subsamples object, or a summary.subsamples object, does
not contain the subsampled count matrix at each depth (as
it would take too much space and is rarely used). However,
as it saves the random seed used to generate the count
matrix, the count matrix at any depth can be retrieved.
This can be done for a subsamples object \code{ss} by
retrieving the seed with \code{getSeed(ss)}. If this seed
is provided to the subsample function, then the same
matrices will be generated when the proportion is the same.

This is useful for adding additional methods or subsampling
depths to an existing subsamples object (after which they
can be combined with \code{combineSubsamples}).
}
\examples{
data(hammer)

hammer.counts = hammer
}

