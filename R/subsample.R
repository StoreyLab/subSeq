#' Subsample reads and perform statistical testing on each sample
#'
#' Perform subsampling at multiple proportions on a matrix of count
#' data representing mapped reads across multiple samples in many
#' genes. For each sample, perform some statistical operations.
#'
#' @param counts Matrix of unnormalized counts
#' @param proportions Vector of subsampling proportions in (0, 1]
#' @param method One or more methods to be performed at each subsample,
#' such as edgeR or DESeq (see Details)
#' @param replications Number of replications to perform at each depth
#' @param seed An initial seed, which will be stored in the output
#' so that any individual simulation can be reproduced.
#' @param qvalues Whether q-values should be calculated for multiple hypothesis
#' test correction at each subsample.
#' @param env Environment in which to find evaluate additional hander functions
#' that are given by name
#' @param ... Other arguments given to the handler, such as \code{treatment}
#'
#' @return A subsample S3 object, which is a data.table containing
#'
#' \item{pvalue}{A p-value calculated for each gene by the handler}
#' \item{coefficient}{An effect size (usually log fold change) calculated
#' for each gene by the handler}
#' \item{ID}{gene ID}
#' \item{count}{the number of reads to this specific gene in this subsample}
#' \item{depth}{the overall sequencing depth of this subsample}
#' \item{method}{the method used (the name of the handler)}
#'
#' @details Method represents the name of a handler function, which can be
#' custom-written by the user.
#'
#' If a gene has a count of 0 at a particular depth, we set the p-value to 1
#' and the coefficient to 0 to stay consistent between programs. If the gene has
#' a count that is not 0 but the p-value is NA, we set the p-value to 1 but
#' keep the estimated coefficient.
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
#' @import data.table
#' @importFrom dplyr group_by do mutate filter
#' @import magrittr
#' @importFrom qvalue qvalue
#' @importFrom Biobase exprs pData
#' @export
subsample <-
    function(counts, proportions, method="edgeR", replications=1,
             seed=NULL, qvalues = TRUE, env=parent.frame(), ...) {
        # error checking
        if (length(proportions) == 0) {
            stop("No proportions to sample")
        }
        if (any(proportions > 1 | proportions == 0)) {
            stop("Proportions must be in range (0, 1]")
        }

        # check that the counts is an unnormalized numeric matrix
        counts = as.matrix(counts)
        error = max(abs(round(counts) - counts))
        if (error > 1e-4 | any(counts < 0)) {
            stop("Counts should be unnormalized integer counts (not, e.g., RPKM)")
        }
        if (is.null(seed)) {
            # come up with a random initial seed (can't use current one, since
            # the point is to save it for later)
            seed = sample(.Machine$integer.max, 1)
        }

        # create a list mapping function name to function
        if (is.function(method)) {
            # if given a single function, make that into a list
            methods = list(method)
            names(methods) = deparse(substitute(method))
        }
        else {
            methods = lapply(method, function(m) {
                # function can be given directly: otherwise, check subSeq package
                # then caller's env
                if (is.function(m)) {
                    handler = m
                }
       #         else if (exists(m, mode="function", envir=as.environment("package:subSeq"))) {
        #            handler = get(m, mode="function", envir=as.environment("package:subSeq"))
         #       }
                else if (exists(m, mode="function", envir=env)) {
                    handler = get(m, mode="function", envir=env)
                }
                else if (m %in% c("edgeR", "voomLimma", "DESeq2", "edgeR.glm")) {
                    handler = get(m, mode="function")
                }
                else {
                    stop(paste("Could not find handler", m))
                }
                handler
            })
            names(methods) = method
        }

        # perform one for each method x proportion x replication
        params = expand.grid(method=names(methods), proportion=proportions, replication=1:replications)
        # don't need replications of full depth (will be identical)
        params = params %>% filter(!(proportion == 1 & replication > 1))

        # apply method to each subsample the specified number of times
        #m.ret = as.data.table(do.call(rbind, lapply(1:nrow(prop.reps),
        perform.subsampling <- function(method, proportion, replication) {
            # generate subcounts and use handler
            subcounts = generateSubsampledMatrix(counts, proportion, seed, replication)
            id = which(rowSums(subcounts) >= 5)
            subcounts = subcounts[id,] ### Filter counts at zero
            if (length(id) == 0) return("Error: counts too low at subsampling proportion")
            handler = methods[[method]]
            ret = handler(subcounts, ...)
            # add gene names (ID) and per-gene counts
            # If there is one handler row per gene, the remaining collumns of ret, "ID" and "counts", can be infered.
            infer.per.gene = dim( ret)[1] == dim( subcounts)[1]
            if ( !any( "ID" == colnames(ret))){
              if (infer.per.gene){
                ret$ID = rownames(subcounts)
              } else {
                stop("if a handler doesn't return one row per gene then it must specify an ID collumn")
              }
            }
            if ( !any( "count" == colnames(ret))){
              if (infer.per.gene){
                ret$count = as.integer(rowSums(subcounts))
              } else {
                ret$count = NA
              }
            }
            ret$depth = sum(subcounts)

            # in any cases of no reads, fix coefficient/pvalue to 0/1
            if ("count" %in% colnames(ret)) {
                if ("pvalue" %in% colnames(ret)) {
                    ret$pvalue[ret$count == 0 | is.na(ret$pvalue) | ret$pvalue == 1] = NA #new qvalue accepts NA values
                }
                if ("coefficient" %in% colnames(ret)) {
                    ret$coefficient[ret$count == 0 | is.infinite(ret$coefficient)] = 0
                }
            }
            ret
        }

        ret = params %>% group_by(method, proportion, replication) %>%
            do(perform.subsampling(.$method, .$proportion, .$replication))



        ## cleanup
        if (qvalues) {
            # calculate q-values
            ret0 = ret %>% filter(proportion == 1) %>% group_by(method) %>%
                  summarize(pi0=qvalue::qvalue(pvalue, lambda = seq(0.05,0.9, 0.05))$pi0) %>% group_by()
            ret = ret %>% inner_join(ret0, by = c("method"))
            ret = ret %>% group_by(proportion, method, replication) %>%
                mutate(qvalue=qvalue.fixedpi0(pvalue, pi0s=unique(pi0))) %>% group_by()

        }

        # turn into a subsamples object
        ret = as.data.table(as.data.frame(ret))
        class(ret) = c("subsamples", class(ret))
        attr(ret, "seed") = seed

        return(ret)
    }

# when calculating pi0, it is prudent to filter out cases where p-values are exactly 1.
# for example, methods of differential expression tend to give p-values of 1 when
# the read count is 0. This would lead to over-estimating pi0 and thus the resulting
# q-values, even though those genes have no impact on the rest of the genes.

qvalue.fixedpi0 <- function (p, fdr.level = NULL, pfdr = FALSE, pi0s=NULL, ...)
{
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  }
  else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level >
                                   1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  m <- length(p)
  u <- order(p)
  v <- rank(p, ties.method = "max")
  if (pfdr) {
    qvals <- (pi0s * m * p)/(v * (1 - (1 - p)^m))
  }
  else {
    qvals <- (pi0s * m * p)/v
  }
  qvals[u[m]] <- min(qvals[u[m]], 1)
  for (i in (m - 1):1) {
    qvals[u[i]] <- min(qvals[u[i]], qvals[u[i + 1]])
  }
  qvals_out[rm_na] <- qvals
  lfdr <- lfdr(p = p, pi0 = pi0s, ...)
  lfdr_out[rm_na] <- lfdr
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level))
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0s, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out)
  }
  return(qvals_out)
}
