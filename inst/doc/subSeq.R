## ----foo,include=FALSE,echo=FALSE-------------
options(keep.source = TRUE, width = 48)
foo <- packageDescription("subSeq")

## ----load_subSeq------------------------------
library(subSeq)
data(hammer)

## ----download_hammer, eval=FALSE--------------
#  load(url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData"))
#  hammer = hammer.eset

## ----setup_hammer, dependson="load_subSeq"----
hammer.counts = hammer@assayData$exprs[, 1:4]
hammer.design = hammer@phenoData@data[1:4, ]
hammer.counts = hammer.counts[rowSums(hammer.counts) >= 5, ]

## ----proportions------------------------------
proportions = 10^seq(-2, 0, .5)
proportions

## ----subSeq_example, dependson=c("setup_hammer", "proportions")----
subsamples = subsample(hammer.counts, proportions, method=c("edgeR", "voomLimma"), treatment=hammer.design$protocol)

## ----show_subsamples, dependson="subSeq_edgeR"----
options(width=40)
subsamples

## ----summary_subsamples, dependson="subSeq_edgeR"----
subsamples.summary = summary(subsamples)
subsamples.summary

## ----plot_subsamples, dependson="summary_subsamples", echo=FALSE----
plot(subsamples.summary)

## ----custom_ggplot2, dependson="summary_subsamples", out.height="3in", out.width="3in"----
library(ggplot2)
ggplot(subsamples.summary, aes(x=depth, y=percent, col=method)) + geom_line()

## ----custom_ggplot2_2, dependson=c("summary_subsamples", "custom_ggplot2"), out.height="3in", out.width="3in"----
ggplot(subsamples.summary, aes(x=depth, y=pearson, col=method)) + geom_line()

## ----subsamples_myMethod, eval=FALSE----
#  subsamples = subsample(hammer.counts, proportions, method=c("edgeR", "DESeq2", "myMethod"), treatment=hammer.design$protocol)

## ----subsamples_more, dependson="subSeq_example"----
seed = getSeed(subsamples)

subsamples.more = subsample(hammer.counts, proportions, method=c("edgeR.glm"), mm=model.matrix(~ hammer.design$protocol), seed=seed)

## ----subsamples_more_combine, dependson="subsamples_more"----
subsamples.combined = combineSubsamples(subsamples, subsamples.more)
plot(summary(subsamples.combined))

## ----generate_subsampled, dependson="subsamples_more"----
submatrix = generateSubsampledMatrix(hammer.counts, .1, seed=seed)
dim(submatrix)
sum(submatrix)

## ----subSample_DEXSeq, dependson="load_SubSeq", results="hide"----
data("pasillaExons", package="pasilla")
exon.counts = pasillaExons@assayData$counts
design = pasillaExons@phenoData@data[c("condition", "type")]
geneIDs = pasillaExons@featureData@data$geneID
exonIDs = pasillaExons@featureData@data$exonID

subsamples.exon = subsample(exon.counts, proportions, method="DEXSeq", design=design, geneIDs=geneIDs, exonIDs=exonIDs)

## ----show_subsamples_exon, dependson="subSample_DEXSeq"----
head(subsamples.exon[proportion == 1])

