subSeq: Subsampling of high-throughput sequencing count data
=======

When you use a RNA-Seq differential expression method, such as edgeR or DESeq2, you can answer a couple of biological questions:

* What genes are differentially expressed?
* What gene sets show differential expression?

However, what if we're interested in questions of experimental design:

* Do I have enough reads to detect most of the biologically relevant differences?
* If I run a similar experiment, should I run additional lanes to obtain more reads? Can I multiplex even more samples and work with *fewer* reads? 

One way to help answer these questions is to pretend you have *fewer* reads than you do, and to see how your results (the number of significant genes, your estimates of their effects, and so on) change. If you can achieve the same results with just 10% of your reads, it indicates that (when using your particular analysis method to answer your particular question) the remaining 90% of the reads added very little. In turn, if your conclusions changed considerably between 80% and 100% of your reads, it is likely they would change more if you added additional reads.

See also [subSeq: Determining appropriate sequencing depth through efficient read subsampling](http://bioinformatics.oxfordjournals.org/content/early/2014/09/03/bioinformatics.btu552.abstract).

Installation
-------------

First install the Bioconductor dependencies:

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("limma", "edgeR", "DESeq2", "DEXSeq", "pasilla"))

Then install the [devtools](https://github.com/hadley/devtools) package, and use it to install the [qvalue 2.0](https://github.com/StoreyLab/qvalue) and subSeq packages. 

    install.packages("devtools")
    library(devtools)
    install_github("jdstorey/qvalue")
    install_github("StoreyLab/subSeq", build_vignettes = TRUE)

Vignette
---------------------

Once you've installed the package, you can access the vignette with

    library(subSeq)
    vignette("subSeq")

You can also run the package's unit tests with

    library(testthat)
    test_package("subSeq")

If you run into a problem or have a question about the software's usage, please open a [GitHub issue](https://github.com/StoreyLab/subSeq/issues).
