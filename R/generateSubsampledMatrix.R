#' Generate the read matrix corresponding to a particular level
#' 
#' @description
#' 
#' Generate a subsampled matrix from an original count matrix. This can be used
#' to perform read subsampling analyses, (though generally the \code{subsample}
#' function is recommended).
#' 
#' It is also useful for reproducing the results of an earlier run (see Details).
#' 
#' @param counts Original matrix of read counts
#' @param proportion The specific proportion to subsample
#' @param seed A subsampling seed, which can be extracted from a subsamples
#' or summary.subsamples object. If not given, doesn't set the seed.
#' @param replication Replicate number: allows performing multiple deterministic 
#' replications at a given subsampling proportion
#' 
#' @details
#' 
#' A subsamples object, or a summary.subsamples object, does not contain the
#' subsampled count matrix at each depth (as it would take too much space and
#' is rarely used). However, as it saves the random seed used to generate the
#' count matrix, the count matrix at any depth can be retrieved. This can be
#' done for a subsamples object \code{ss} by retrieving the seed with
#' \code{getSeed(ss)}. When given along with the original counts, the
#' proportion, and the replication number (if more than one subsampling was done
#' at each proportion) this produces the same matrix as was used in the analysis.
#' 
#' The seed is calculated deterministically using an md5 hash of three combined
#' values: the global seed used for the subsampling object, the subsampling
#' proportion, and the replication # for that proportion.
#' 
#' @examples
#' 
#' data(hammer)
#' 
#' hammer.counts = Biobase::exprs(hammer)[, 1:4]
#' hammer.design = Biobase::pData(hammer)[1:4, ]
#' hammer.counts = hammer.counts[rowSums(hammer.counts) >= 5, ]
#' 
#' ss = subsample(hammer.counts, c(.01, .1, 1), treatment=hammer.design$protocol,
#'                  method=c("edgeR", "DESeq2", "voomLimma"))
#' 
#' seed = getSeed(ss)
#' 
#' # generate the matrices used at each subsample
#' subm.01 = generateSubsampledMatrix(hammer.counts, .01, seed)
#' subm.1 = generateSubsampledMatrix(hammer.counts, .1, seed)
#' 
#' @import digest
#' @return subsamples matrix at specified subsampling proportion
#' @export
generateSubsampledMatrix <- function(counts, indeces, proportions, seed, replication=1) {
  if (!missing(seed)) {
    # calculate seed using an md5 hash of the global seed and the
    # subsampling proportions
    s = readBin(digest(c(seed, proportions, replication, indeces), raw=TRUE), "integer")
    set.seed(s)
  }
  else if (is.null(seed)) {
    stop(paste("Given a NULL seed: Probably was an error retrieving",
               "the seed from the desired object"))
  }
  # apply random binomial sampling to each cell
  # keep row names
  rns <- rownames(counts)[indeces]
  n = nrow(counts)
  ret <- apply( t( as.array( as.numeric(inds))), 2, function(ind){
    rbinom(n, counts[ , ind], proportions[ ind])
  })

  rownames(ret) <- rns
  ret
}

#' Extract the global random seed from a subsamples object
#' 
#' @description
#' 
#' A subsamples object, or a summary.subsamples object, does not contain the
#' subsampled count matrix at each depth (as it would take too much space and
#' is rarely used). However, as it saves the random seed used to generate the
#' count matrix, the count matrix at any depth can be retrieved. This can be
#' done for a subsamples object \code{ss} by retrieving the seed with
#' \code{getSeed(ss)}. If this seed is provided to the subsample function, then
#' the same matrices will be generated when the proportion is the same.
#' 
#' This is useful for adding additional methods or subsampling depths to an
#' existing subsamples object (after which they can be combined with
#' \code{combineSubsamples}).
#' 
#' @param ss A subsamples object, returned from the \code{subsample} function,
#' or a summary of that object
#' 
#' @examples
#' 
#' data(hammer)
#' 
#' hammer.counts = Biobase::exprs(hammer)[, 1:4]
#' hammer.design = Biobase::pData(hammer)[1:4, ]
#' hammer.counts = hammer.counts[rowSums(hammer.counts) >= 5, ]
#' 
#' ss = subsample(hammer.counts, c(.01, .1, 1), treatment=hammer.design$protocol,
#'                  method=c("edgeR", "DESeq2", "voomLimma"))
#' 
#' seed = getSeed(ss)
#' @return get seed of subSeq object
#' @export
getSeed <- function(ss) {
  attr(ss, "seed")
}
