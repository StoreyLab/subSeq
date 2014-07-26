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
#' @param progress whether to show a progress bar
#' @param seed An initial seed, which will be stored in the output
#' so that any individual simulation can be reproduced.
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
#' hammer.counts = hammer@assayData$exprs[, 1:4]
#' hammer.design = hammer@phenoData@data[1:4, ]
#' hammer.counts = hammer.counts[rowSums(hammer.counts) >= 5, ]
#' 
#' ss = subsample(hammer.counts, c(.01, .1, 1), treatment=hammer.design$protocol,
#'                  method=c("edgeR", "DESeq2", "voomLimma"))
#' 
#' @import data.table
#' @importFrom qvalue qvalue
#' 
#' @export
subsample <-
    function(counts, proportions, method="edgeR", replications=1, progress=FALSE,
             seed=NULL, env=parent.frame(), ...) {
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
        
        # perform one for each method and proportion, creating a data.table for each, and combine
        if (is.function(method)) {
            # if given a single function, make that into a list
            method.name = deparse(substitute(method))
            method = list(method)
        }
        results = lapply(method, function(m) {        
            if (is.function(m)) {
                # function can be given directly
                handler = m
                m = method.name  # make m a string
            }
            else if (exists(m, mode="function", envir=as.environment("package:subSeq"))) {
                # first look in subSeq package
                handler = get(m, mode="function", envir=as.environment("package:subSeq"))
            }
            else if (exists(m, mode="function", envir=env)) {
                # then look in the caller's environment
                handler = get(m, mode="function", envir=env)
            }
            else {
                stop(paste("Could not find handler", m))
            }

            prop.reps = expand.grid(proportions, 1:replications)
            # don't need replications of full depth (will be identical)
            prop.reps = prop.reps[!(prop.reps$Var1 == 1 & prop.reps$Var2 > 1), ]

            if (progress) {
                cat(paste0("Subsampling with method ", m, ":\n"))
                pb <- txtProgressBar(min=0, max=nrow(prop.reps))
            }
            
            # apply method to each subsample the specified number of times
            m.ret = as.data.table(do.call(rbind, lapply(1:nrow(prop.reps), function(i) {
                if (progress) {
                    setTxtProgressBar(pb, i)
                }
                
                proportion = prop.reps[i, 1]
                replication = prop.reps[i, 2]
                
                subcounts = generateSubsampledMatrix(counts, proportion, seed, replication)
                ret = handler(subcounts, ...)
                
                if (NROW(ret) == NROW(counts)) {
                    # if the output is the same size as the input, assume there
                    # is a one-to-one gene correspondence
                    # add gene names (ID) and per-gene counts
                    ret$ID = factor(rownames(counts))
                    ret$count = as.integer(rowSums(subcounts))
                }
                else {
                    if (is.null(ret$ID)) {
                        stop(paste("If handler does not return one row per gene,",
                                   "it must include an ID column."))
                    }
                }
                ret$depth = sum(subcounts)    
                ret$method = m
                ret$proportion = proportion
                ret$replication = replication
                
                # in any cases of no reads, fix coefficient/pvalue to 0/1
                if ("count" %in% colnames(ret)) {
                    ret$pvalue[ret$count == 0 | is.na(ret$pvalue)] = 1
                    ret$coefficient[ret$count == 0 | is.infinite(ret$coefficient)] = 0
                }
                ret
            })))
            
            ## cleanup
            # make sure method is a factor
            m.ret$method = factor(m.ret$method)
            
            if (progress) {
                close(pb)
                cat("Calculating q-values... ")
            }
            
            # calculate q-value
            m.ret[, qvalue:=qvalue.filtered1(pvalue), by=depth]
            
            class(m.ret) = c("subsamples", class(m.ret))
            
            attr(m.ret, "seed") = seed
            
            if (progress) {
                cat("done.\n")
            }
            
            m.ret
        })
        
        ret = do.call(combineSubsamples, results)
        return(ret)
    }

# when calculating pi0, it is prudent to filter out cases where p-values are exactly 1.
# for example, methods of differential expression tend to give p-values of 1 when 
# the read count is 0. This would lead to over-estimating pi0 and thus the resulting
# q-values, even though those genes have no impact on the rest of the genes.
qvalue.filtered1 = function(p) {
    # given a vector of p-values
    q.without1 = qvalue(p[p < 1])
    # when p == 1, the FDR is pi0
    q = rep(q.without1$pi0, length(p))
    q[p < 1] = q.without1$qvalue
    q
}
