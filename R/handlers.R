## Handler functions for subSeq

# each of these takes a count matrix and any number of extra arguments (often
# a single vector representing the treatment variable), and returns a data
# frame with one row for each gene, and the columns "pvalue" and "coefficient" (a
# handler can include additional columns as well)

# one can also write custom handler functions and pass them to subsample.

#' handler for performing edgeR on subsampled RNA-Seq data
#' 
#' This encapsulates the edgeR \code{exactTest} function. It takes a count
#' matrix, potentially subsampled, a vector describing the treatment of each
#' sample, and the pair of conditions to compare, and performs a negative
#' binomial exact test for differential expression between each pair of
#' conditions.
#' 
#' @param count.matrix Matrix of unnormalized counts
#' @param treatment vector describing the treatment of each sample (column)
#' @param pair Pair of treatments to compare. If there are only two treatments,
#' this is optional.
#' 
#' @return A data frame with one row for each gene, with the columns
#' 
#' \item{coefficient}{log2 fold change between the conditions}
#' \item{pvalue}{P-value using a negative binomial test comparing the conditions}
#' 
#' @export
edgeR <-
    function(count.matrix, treatment, pair=NULL) {
        require("edgeR")
        treatment = factor(treatment)
        if (is.null(pair)) {
            if (length(levels(treatment)) > 2) {
                stop(paste("If there are more than two levels of treatment,",
                           "must be explicitly given pair to compare"))
            }
            pair = levels(treatment)
        }
        
        d = DGEList(counts=count.matrix, group=treatment)
        d = estimateCommonDisp(d)
        d = estimateTagwiseDisp(d)
        edgeR.result = exactTest(d, pair=pair)$table
        ret = data.frame(coefficient=edgeR.result$logFC, pvalue=edgeR.result$PValue)
        ret
    }


#' handler for performing an edgeR glm on subsampled RNA-Seq data
#' This encapsulates the edgeR \code{glmFit} and \code{glmLRT} 
#' functions. It takes a count matrix, potentially subsampled, a vector describing 
#' the group (usually treatment) of each sample, and a model matrix, and performs
#' a negative binomial regression to find each coefficient, followed by a
#' likelihood ratio test to test for significance of each
#' 
#' @param count.matrix Matrix of unnormalized counts
#' @param group vector describing the group of each sample
#' @param mm model matrix describing one or more regressors to fit
#' @param column on which column of the model matrix should the likelihood
#' ratio test should be performed
#' 
#' @return A data frame with one row for each gene, with the columns
#' 
#' \item{effect.size}{log fold change between the conditions}
#' \item{pvalue}{P-value using a likelihood ratio test comparing the conditions}
#' 
#' @export
edgeR.glm <-
    function(count.matrix, mm, column=2, group=NULL) {
        require("edgeR")
        if (is.null(group)) {
            group = rep(1, ncol(count.matrix))
        }
        d = DGEList(counts=count.matrix, group=group)
        d = estimateGLMTrendedDisp(d, mm)
        d = estimateGLMTagwiseDisp(d, mm)
        fit = glmFit(d, mm, dispersion=d$tagwise.dispersion)
        lrt = glmLRT(fit, column)
        ret = lrt$table[, c("logFC", "PValue")]
        colnames(ret) = c("coefficient", "pvalue")
        ret
    }


#' handler for performing voom normalization and then limma linear modeling
#' on subsampled RNA-Seq data
#' 
#' This encapsulates the voom normalization and limma differential expression
#' testing functions.
#' 
#' @param count.matrix Matrix of unnormalized counts
#' @param treatment vector describing the treatment of each sample (column)
#' @param ... Additional arguments to be passed to the voom normalization
#' 
#' @return A data frame with one row for each gene, with the columns
#' 
#' \item{coefficient}{log2 fold change between the conditions}
#' \item{pvalue}{P-value using a linear model of expression ~ treatment}
#' 
#' @export
voomLimma <-
    function(count.matrix, treatment, ...) {
        require("limma")
        v = voom(count.matrix, ...)
        f = lmFit(v, model.matrix(~ treatment))
        e = eBayes(f)
        data.frame(coefficient=f$coefficients[, 2], pvalue=e$p.value[, 2])
    }


#' handler for performing DESeq2 on subsampled RNA-Seq data
#' 
#' This encapsulates the DESeq2 \code{nbinomWaldTest} function. It takes a count
#' matrix, potentially subsampled, a vector describing the treatment of each
#' sample, and the pair of conditions to compare, and performs a negative
#' binomial Wald test for differential expression between each pair of
#' conditions.
#' 
#' @param count.matrix Matrix of unnormalized counts
#' @param treatment vector describing the treatment of each sample (column)
#' 
#' @return A data frame with one row for each gene, with the columns
#' 
#' \item{coefficient}{fold change between the conditions}
#' \item{pvalue}{P-value using a negative binomial test comparing the conditions}
#' 
#' @export
DESeq2 <-
    function(count.matrix, treatment) {
        require("DESeq2")
        dds = DESeqDataSetFromMatrix(countData = count.matrix,
                                     colData = data.frame(treatment=treatment),
                                     design = ~ treatment)
        dds = DESeq2::estimateSizeFactors(dds)
        dds = DESeq2::estimateDispersions(dds, quiet=TRUE)
        dds = DESeq2::nbinomWaldTest(dds, quiet=TRUE)
        # converting a DataFrame to a data frame requires a little trick:
        convert.df = selectMethod("as.data.frame", "DataFrame")
        ret = convert.df(results(dds)[c("log2FoldChange", "pvalue")])
        colnames(ret) = c("coefficient", "pvalue")
        ret
    }


#' handler for performing DEXSeq on subsampled RNA-Seq data
#' 
#' This encapsulates the DEXSeq dispersion estimation and \code{testForDEU}
#' functions. It takes a count matrix, potentially subsampled, and a data table
#' representing the study design, and tests for differential exon usage.
#' 
#' @param count.matrix Matrix of unnormalized counts. The row names should
#' represent the exon names
#' @param design Study design data frame for DEXSeq
#' @param geneIDs vector of gene IDs- one for each row of the count matrix
#' (repeated for each exon in a gene)
#' @param exonIDs vector of exon IDs- one for each row of the count matrix
#' 
#' @return A data frame with one row for each gene. This contains all columns
#' from the \link{DEXSeq::DEXSeqResults} function in DEXSeq, except that the 
#' \code{log2fold_untreated_treated} column is been changed to \code{coefficient}
#' and an \code{ID} column is added, with the format "geneID:exonID"
#' 
#' @examples
#' 
#' data( "pasillaExons", package="pasilla" )
#' count.matrix = pasillaExons@assayData$counts
#' 
#' design = pasillaExons@phenoData@data[c("condition", "type")]
#' geneIDs = pasillaExons@featureData@data$geneID
#' exonIDs = pasillaExons@featureData@data$exonID
#' 
#' proportions = c(.1, 1)
#' ss.exon = subsample(exon.counts, proportions, design=pasillaExons@phenoData@data, geneIDs)
#' 
#' DEXSeq(exon.counts, design, geneIDs, exonIDs)
#' 
#' @export
DEXSeq <-
    function(count.matrix, design, geneIDs, exonIDs) {
        require("DEXSeq")
        #formula.dispersion = count ~ sample + (exon + type) * condition
        ecs = DEXSeqDataSet(count.matrix, design, featureID=exonIDs, groupID=geneIDs)
        ecs = estimateSizeFactors(ecs)
        ecs = estimateDispersions(ecs, fitType='local')
        ecs = testForDEU(ecs)
        
        ecs = estimateExonFoldChanges(ecs)
        result = as.data.table(as.data.frame(DEXSeqResults(ecs)))

        # change some names
        #coefname = tail(colnames(result), 1)
        setnames(result, "log2fold_untreated_treated", "coefficient")
        result$ID = paste(result$groupID, result$featureID, sep=":")
        result
    }
