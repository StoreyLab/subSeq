edgeR <-
    function(count.matrix, treatment, pair=NULL) {
        treatment = factor(treatment)
        if (is.null(pair)) {
            if (length(levels(treatment)) > 2) {
                stop(paste("If there are more than two levels of treatment,",
                           "must be explicitly given pair to compare"))
            }
            pair = levels(treatment)
        }

        d = edgeR::DGEList(counts=count.matrix, group=treatment)
        d = edgeR::calcNormFactors(d)
        d = edgeR::estimateCommonDisp(d)
        d = edgeR::estimateTagwiseDisp(d)
        edgeR.result = edgeR::exactTest(d, pair=pair)$table
        ret = data.frame(coefficient=edgeR.result$logFC, pvalue=edgeR.result$PValue)
        ret
    }


edgeR.glm <-
    function(count.matrix, mm, column=2, group=NULL) {
        if (is.null(group)) {
            group = rep(1, ncol(count.matrix))
        }
        d = edgeR::DGEList(counts=count.matrix, group=group)
        d = edgeR::calcNormFactors(d)
        d = edgeR::estimateGLMTrendedDisp(d, mm)
        d = edgeR::estimateGLMTagwiseDisp(d, mm)
        fit = edgeR::glmFit(d, mm, dispersion=d$tagwise.dispersion)
        lrt = edgeR::glmLRT(fit, column)
        ret = lrt$table[, c("logFC", "PValue")]
        colnames(ret) = c("coefficient", "pvalue")
        ret
    }

voomLimma <-
    function(count.matrix, treatment, ...) {
        d = edgeR::DGEList(counts=count.matrix, group=treatment)
        d = edgeR::calcNormFactors(d)
        v = limma::voom(d, ...)
        f = limma::lmFit(v, model.matrix(~ treatment))
        e = limma::eBayes(f)
        data.frame(coefficient=f$coefficients[, 2], pvalue=e$p.value[, 2], t.test = e$t[,2])
    }

DESeq2 <-
    function(count.matrix, treatment) {
        dds = DESeq2::DESeqDataSetFromMatrix(countData = count.matrix,
                                     colData = data.frame(treatment=treatment),
                                     design = ~ treatment)
        dds = DESeq2::estimateSizeFactors(dds)
        dds = DESeq2::estimateDispersions(dds, quiet=TRUE)
        dds = DESeq2::nbinomWaldTest(dds, quiet=TRUE)
        # converting a DataFrame to a data frame requires a little trick:
        convert.df = selectMethod("as.data.frame", "DataFrame")
        ret = convert.df(DESeq2::results(dds)[c("log2FoldChange", "pvalue")])
        colnames(ret) = c("coefficient", "pvalue")
        ret
    }

DEXSeq <-
    function(count.matrix, design, geneIDs, exonIDs) {
        #formula.dispersion = count ~ sample + (exon + type) * condition
        ecs = DEXSeq::DEXSeqDataSet(count.matrix, design, featureID=exonIDs, groupID=geneIDs)
        ecs = DEXSeq::estimateSizeFactors(ecs)
        ecs = DEXSeq::estimateDispersions(ecs, fitType='local')
        ecs = DEXSeq::testForDEU(ecs)

        ecs = DEXSeq::estimateExonFoldChanges(ecs)
        result = as.data.table(as.data.frame(DEXSeq::DEXSeqResults(ecs)))

        # change some names
        #coefname = tail(colnames(result), 1)
        setnames(result, "log2fold_untreated_treated", "coefficient")
        result$ID = paste(result$groupID, result$featureID, sep=":")
        result
    }
