\name{plot.summary.subsamples}
\alias{plot.summary.subsamples}
\title{plot metrics as a function of subsampled read depth}
\usage{
\method{plot}{summary.subsamples}(x, ...)
}
\arguments{
  \item{x}{a \code{summary.subsamples} object}

  \item{...}{further arguments passed to or from other
  methods.}
}
\description{
Plot the number of genes found significant, the Spearman
correlation of the effect size estimates with the full
experiment, and the empirical false discovery rate as a
function of the subsampled read depth. This determines
whether these metrics saturate, which indicates that the
experiment has an appropriate sequencing depth.
}

