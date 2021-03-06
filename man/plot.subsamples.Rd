% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plot.summary.subsamples.R
\name{plot.subsamples}
\alias{plot.subsamples}
\title{plot metrics as a function of subsampled read depth}
\usage{
\method{plot}{subsamples}(x, ...)
}
\arguments{
\item{x}{a \code{subsamples} object}

\item{...}{further arguments passed to or from other methods.}
}
\value{
plot a subSeq object
}
\description{
Plot the number of genes found significant, the Spearman correlation of the
effect size estimates with the full experiment, and the empirical false
discovery rate as a function of the subsampled read depth. This determines
whether these metrics saturate, which indicates that the experiment has an
appropriate sequencing depth.
}
\details{
This is an alias for the \link{plot.summary.subsamples} function, so that
plotting can be done directly on the subsamples object. We recommend using
\code{summary(ss)} first, so that the summary operation does not have to
be performed each time the figure is plotted, and so the summary object
can be examined on its own.
}
\examples{
if (interactive()) {
# import the subsampling object (see ?subsample to see how ss is created)
data(ss)

# plot subsample object
plot(ss)
}
}

