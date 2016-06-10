This fork is under development and is not yet stable.

RNA-seq experiments provide the opportunity to test for differential expression across many genes in a single experiment. Variability in the results of such tests depends on the number of observations that are made. This fork examines such variability by picking RNA samples from the original experiment and using their associated reads to test for differential expression. This fork borrows heavily from the methods of [subSeq](https://github.com/StoreyLab/subSeq), which simulates data-sets with decreased numbers of reads. 

References:
[subSeq: Determining appropriate sequencing depth through efficient read subsampling](http://bioinformatics.oxfordjournals.org/content/early/2014/09/03/bioinformatics.btu552.abstract).