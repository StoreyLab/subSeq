#' Subsample reads and perform statistical testing on each sample
#'
#' Perform subsampling at multiple proportions on a matrix of count
#' data representing mapped reads across multiple samples in many
#' genes. For each sample, perform some statistical operations.
#'
#' @param counts Matrix of unnormalized counts
#' @param treatments Vector (factor) of experimental treatments corresponding to
#' collumns of counts
#' @param proportions Vector of subsampling proportions in (0, 1]
#' @param bioReplicates Vector specifying number of samples from each treatment used in subsampling
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
#' ss = subsample(counts=        hammer.counts,
#'                treatments =   hammer.design$protocol,
#'                proportions=   c(.01, .1, 1),
#'                bioReplicates= c(2),
#'                method=        c("edgeR", "voomLimma"))
#'
#' @import data.table
#' @importFrom dplyr group_by do mutate filter
#' @import magrittr
#' @importFrom qvalue qvalue lfdr
#' @importFrom Biobase exprs pData
#' @export
subsample <-
  function(counts, treatments, proportions, bioReplicates, method="edgeR", replications=1,
           replacement= FALSE, ballanced.proportions= FALSE, seed=NULL, qvalues = TRUE, env=parent.frame(), ...) {
    # error checking
    if (length(proportions) == 0) {
      stop("No proportions to sample")
    }
    if (length( bioReplicates) == 0){
      stop("No numbers of bioReplicates to sample")
    }
    if (any(proportions > 1 | proportions == 0)) {
      stop("Proportions must be in range (0, 1]")
    }
    if ( length( treatments) != dim( counts)[2]){
      stop("Collumns of counts and elements of treatments should correspond
           so they should have the same length")
    }
    if ( any( !( bioReplicates %in% seq(2, dim(counts)[2])))){
      stop( "The number of biological replicates must be greater
            than one and less than the number of samples")
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
    } else {
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
    
    # perform one for each method x proportion x bioReplicates x replication
    params = expand.grid(method=names(methods), proportion=proportions, biological.replicate= bioReplicates, replication=1:replications)
    # don't need replications of full depth (will be identical)
    params = params %>% filter(!(proportion == 1 & replication > 1))
    # apply method to each subsample the specified number of times
    #m.ret = as.data.table(do.call(rbind, lapply(1:nrow(prop.reps),
    perform.ballanced.subsampling <- function(method, proportion, bioReplicates, replication) {
      # resample biological replicates
      inds <- getIndecesFromCatagoricalTreatment( treatments, bioReplicates, seed, replacement)
      treatment <- treatments[inds]
      #Get proportions for each index in inds
      #Calculating colsums once would help efficientcy, but it needs to be done in correct place to be readable
      if( ballanced.proportions == TRUE){
        total.counts <- colSums( counts)
        #ind.proportion * total.counts == min( total.counts)
        ind.proportions <- proportion * min( total.counts) / total.counts[ inds]
      } else {
        ind.proportions <- rep( proportion, length.out= length( inds))
      }
      # subsample reads and use handler
      subcounts = generateSubsampledMatrix(counts, inds, ind.proportions, seed, replication)
      id = which(rowSums(subcounts) >= 5)
      subcounts = subcounts[id,] ### Filter counts at zero
      if (length(id) == 0) return("Error: counts too low at subsampling proportion")
      handler = methods[[method]]
      ret = handler(subcounts, treatment, ...)
      ## add gene names (ID) and per-gene counts
      ## If there is one handler row per gene, the remaining collumns of ret, "ID" and "counts", can be infered.
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
    
    ret = params %>% group_by(method, proportion, biological.replicate, replication) %>%
      do(perform.ballanced.subsampling( .$method, .$proportion, .$biological.replicate, .$replication))
    
    ## cleanup
    if (qvalues) {
      # calculate q-values
      max.proportion <- max( ret$proportion)
      ret0 = ret %>% filter(proportion == max.proportion) %>% group_by(method) %>%
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
  lfdr <- qvalue::lfdr(p = p, pi0 = pi0s, ...)
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
