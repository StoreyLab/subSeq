#' @name hammer
#' 
#' @title ExpressionSet results from Hammer et al 2010
#' 
#' @docType data
#' 
#' @description An ExpressionSet containing the results of the Hammer et al 2010 RNA-Seq
#' study on the nervous system of rats (Hammer et al 2010). This dataset is used in the
#' examples and vignette for the subSeq package.
#' 
#' This was downloaded from the ReCount database of analysis-ready RNA-Seq datasets (Frazee et al 2011).
#' 
#' Hammer, P., Banck, M. S., Amberg, R., Wang, C., Petznick, G., Luo, S., Khrebtukova, I., Schroth, G. P.,
#' Beyerlein, P., and Beutler, A. S. (2010). mRNA-seq with agnostic splice site discovery for nervous system
#' transcriptomics tested in chronic pain. Genome research, 20(6), 847-860.
#' http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877581/
#' 
#' Frazee, A. C., Langmead, B., and Leek, J. T. (2011). ReCount: a multi-experiment resource of analysis-ready
#' RNA-seq gene count datasets. BMC Bioinformatics, 12, 449.
#' http://bowtie-bio.sourceforge.net/recount/
NULL

#' @name ss
#' 
#' @title Subsampling results using the hammer dataset
#' 
#' @docType data
#' 
#' @description The subsample object \code{ss} is the result from applying the \link{subsample} function to the \link{hammer} data set. The hypothesis test was a simple two-sample comparison (control vs. L5 SNL). Voom, DESeq2 and edgeR were used to test for differential expression at three different subsampling proportions: 0.01, 0.1 and 1. Genes with less than 5 counts across all replicates were filtered. For more details on how the object was generated, please see the \link{subsample} function.
#' 
#' The subsample object can then be used to determine whether an experiment has adequate read depth (see \link{plot} and \link{summary} functions). 
#'
#'  Hammer, P., Banck, M. S., Amberg, R., Wang, C., Petznick, G., Luo, S., Khrebtukova, I., Schroth, G. P.,
#' Beyerlein, P., and Beutler, A. S. (2010). mRNA-seq with agnostic splice site discovery for nervous system
#' transcriptomics tested in chronic pain. Genome research, 20(6), 847-860.
#' http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877581/
#' 
#' Frazee, A. C., Langmead, B., and Leek, J. T. (2011). ReCount: a multi-experiment resource of analysis-ready
#' RNA-seq gene count datasets. BMC Bioinformatics, 12, 449.
#' http://bowtie-bio.sourceforge.net/recount/
#' @examples 
#' # import the subsampling object (see ?subsample to see how ss is created)
#' data(ss)
#' 
#' # summarise object
#' sum_ss <- summary(ss)
#' #plot
#' if (interactive()) {
#'   plot(ss)
#' }
NULL
