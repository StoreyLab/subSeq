#' calculate summary statistics for each subsampled depth in a subsamples object
#' 
#' Given a subsamples object, calculate a metric for each depth that summarizes
#' the power, the specificity, and the accuracy of the effect size estimates at
#' that depth.
#' 
#' To perform these calculations, one must compare each depth to an "oracle" depth,
#' which, if not given explicitly, is assumed to be the highest subsampling depth.
#' This thus summarizes how closely each agrees with the full experiment: if very
#' low-depth subsamples still agree, it means that the depth is high enough that
#' the depth does not make a strong qualitative difference.
#' 
#' @param object a subsamples object
#' @param oracle a subsamples object of one depth showing what each depth should 
#' be compared to; if NULL, each will be compared to the highest depth
#' @param FDR.level A false discovery rate used to calculate the number of genes
#' found significant at each level
#' @param average If TRUE, averages over replications at each method+depth
#' combination before returning
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A summary object, which is a \code{data.table}
#' with one row for each subsampling depth, containing the metrics
#' 
#' \item{significant}{number of genes found significant at the given FDR}
#' \item{pearson}{Pearson correlation of the coefficient estimates with the oracle}
#' \item{spearman}{Spearman correlation of the coefficient estimates with the oracle}
#' \item{MSE}{mean squared error between the coefficient estimates at the oracle}
#' \item{estFDP}{estimated FDP: the estimated false discovery proportion, as calculated from the
#' average oracle local FDR within genes found significant at this depth}
#' \item{rFDP}{relative FDP: the proportion of genes found significant at this depth that were not found
#' significant in the oracle}
#' \item{percent}{the percentage of genes found significant in the oracle that
#' were found significant at this depth}
#' 
#' @details Note that selecting average=TRUE averages the depths of the replicates
#' (as two subsamplings with identical proportions may have different depths by
#' chance). This may lead to depths that are not integers.
#' 
#' @examples
#' 
#' data(hammer)
#' 
#' hammer.counts = hammer@@assayData$exprs[, 1:4]
#' hammer.design = hammer@@phenoData@@data[1:4, ]
#' hammer.counts = hammer.counts[rowSums(hammer.counts) >= 5, ]
#' 
#' ss = subsample(hammer.counts, c(.01, .1, 1), treatment=hammer.design$protocol,
#'                  method=c("edgeR", "DESeq2", "voomLimma"))
#' 
#' ss.summary = summary(ss)
#' 
#' @importFrom qvalue lfdr
#' 
#' @export
summary.subsamples <-
function(object, oracle=NULL, FDR.level=.05, average=FALSE, ...) {
    # find the oracle for each method
    object$method = as.character(object$method)
    if (is.null(oracle)) {
        oracles = lapply(unique(object$method), function(m) {
            object[method == m, ][depth == max(depth), ]
        })
    }
    else {
        # oracle is the same for all methods
        oracles = lapply(unique(object$method), function(m) oracle)
    }
    names(oracles) = unique(object$method)

    # calculate lfdr for each oracle
    for (n in names(oracles)) {
        # calculate for all p-values that aren't equal to 1 separately
        p = oracles[[n]]$pvalue
        non1.lfdr = lfdr(p[p != 1])
        oracles[[n]]$lfdr = max(non1.lfdr)
        oracles[[n]]$lfdr[p != 1] = non1.lfdr
    }

    get.stats = function(i, q, coef, m) {
        o = oracles[[as.character(m[1])]][match(i, ID), ]

        # statistics
        valid = !is.na(coef) & !is.na(o$coefficient)
        osig = o$qvalue < FDR.level
        list(significant=sum(q < FDR.level),
              pearson=cor(coef, o$coefficient, use="complete.obs"),
              spearman=cor(coef, o$coefficient, use="complete.obs", method="spearman"),
              MSE=mean((coef[valid] - o$coefficient[valid])^2),
              estFDP=mean(o$lfdr[q < FDR.level]),
              rFDP=mean((o$qvalue > FDR.level)[q < FDR.level]),
              percent=mean(q[osig] < FDR.level))
    }

    ret = object[, get.stats(ID, qvalue, coefficient, method),
                by=c("depth", "proportion", "method", "replication")]

    # any case where none are significant, the estFDP/rFDP should be 0 (not NaN)
    # since technically there were no false discoveries
    ret$estFDP[ret$significant == 0] = 0
    ret$rFDP[ret$significant == 0] = 0

    if (average) {
        ret = ret[, list(depth=mean(depth), significant=mean(significant), pearson=mean(pearson),
                         spearman=mean(spearman), MSE=mean(MSE), rFDP=mean(rFDP),
                         estFDP=mean(estFDP), percent=mean(percent)), by=c("proportion", "method")]
    }
    
    class(ret) = c("summary.subsamples", "data.table", "data.frame")
    attr(ret, "seed") = attr(object, "seed")
    ret
}
