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
#' @import dplyr
#' @import tidyr
#' 
#' @export
summary.subsamples <-
function(object, oracle=NULL, FDR.level=.05, average=FALSE, ...) {
    # find the oracle for each method
    tab = as.data.frame(object)

    tab = tab %>% mutate(method=as.character(method))

    if (is.null(oracle)) {
        # use the highest depth in each method
        oracles = tab %>% group_by(method) %>% filter(depth == max(depth))
    }
    else {
        # oracle is the same for all methods
        oracles = data.frame(method=unique(object$method)) %>%
            group_by(method) %>% do(oracle)
    }

    # calculate lfdr for each oracle
    lfdr1 = function(p) {
        # calculate for all p-values that aren't equal to 1 separately
        non1.lfdr = lfdr(p[p != 1])
        ret = rep(max(non1.lfdr), length(p))
        ret[p != 1] = non1.lfdr
        ret
    }
    oracles = oracles %>% group_by(method) %>% mutate(lfdr=lfdr1(pvalue))

    # combine with oracle
    sub.oracle = oracles %>% select(method, ID, o.pvalue=pvalue, o.qvalue=qvalue,
                                     o.coefficient=coefficient, o.lfdr=lfdr)
    tab = tab %>% inner_join(sub.oracle, by=c("method", "ID"))
    
    # summary operation
    ret = tab %>% group_by(depth, proportion, method, replication) %>%
        mutate(valid=(!is.na(coefficient) & !is.na(o.coefficient))) %>%
        summarize(significant=sum(qvalue < FDR.level),
                  pearson=cor(coefficient, o.coefficient, use="complete.obs"),
                  spearman=cor(coefficient, o.coefficient, use="complete.obs", method="spearman"),
                  MSE=mean((coefficient[valid] - o.coefficient[valid])^2),
                  estFDP=mean(o.lfdr[qvalue < FDR.level]),
                  rFDP=mean((o.qvalue > FDR.level)[qvalue < FDR.level]),
                  percent=mean(qvalue[o.qvalue < FDR.level] < FDR.level))

    # any case where none are significant, the estFDP/rFDP should be 0 (not NaN)
    # since technically there were no false discoveries
    ret = ret %>% mutate(estFDP=ifelse(significant == 0, 0, estFDP)) %>%
        mutate(rFDP=ifelse(significant == 0, 0, rFDP))

    if (average) {
        # average each metric within replications
        ret = ret %>% gather(metric, value, significant:percent) %>% group_by(proportion, method, metric) %>%
            summarize(value=mean(value)) %>% spread(metric, value)
    }
    
    ret = as.data.table(ret)
    class(ret) = c("summary.subsamples", "data.table", "data.frame")
    attr(ret, "seed") = attr(object, "seed")
    attr(ret, "FDR.level") = FDR.level
    ret
}
