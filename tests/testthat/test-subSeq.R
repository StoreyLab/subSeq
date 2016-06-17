#TODO(riley) test replacement and ballanced.proportions features of subseq
### custom expectations ###

expect_na = function(x){
  eval(bquote(expect_equal(sum( !is.na(.(x))), 0)))
}
expect_not_na = function(x){
  eval(bquote(expect_equal(sum( is.na(.(x))), 0)))
}

### setup ###

data(hammer)

hammer.counts = hammer@assayData$exprs[, 1:4]
hammer.design = hammer@phenoData@data[1:4, ]

# create a smaller set of counts to perform (faster) subsampling on
counts = hammer.counts[rowSums(hammer.counts) > 5000, ]
proportions = c(.01, .1, 1)
treatment = hammer.design$protocol
biological.replicate <- 2


context("Individual handlers")

# for handlers, work with a small set that covers a range of sizes
smallcounts = hammer.counts[rowSums(hammer.counts) <= 100, ]
largecounts = hammer.counts[rowSums(hammer.counts) >= 100, ]

testcounts = rbind(smallcounts[sample(nrow(smallcounts), 20), ],
                   largecounts[sample(nrow(largecounts), 280), ])

test.output = function(output, numgenes=NULL) {
    # run tests on the output of a handler
    expect_that(output, is_a("data.frame"))
    expect_true("pvalue" %in% colnames(output))
    expect_true("coefficient" %in% colnames(output))
    expect_true(all(output$pvalue[!is.na(output$pvalue)] >= 0))
    expect_true(all(output$pvalue[!is.na(output$pvalue)] <= 1))

    if (!is.null(numgenes)) {
        expect_equal(nrow(output), numgenes)
    }
}

#test_that("edgeR.glm handler works", {
#    mm = model.matrix(~ treatment)
#    test.output(subSeq::edgeR.glm(testcounts, mm=mm), nrow(testcounts))
#    confounder = rnorm(ncol(testcounts))
#    mm2 = model.matrix(~ treatment + confounder)
#    test.output(subSeq::edgeR.glm(testcounts, mm=mm2), nrow(testcounts))

    # check that the edgeR applied to a random confounder is roughly uniform
    # (and thus that it actually is using that confounder)
#    conf.out = subSeq::edgeR.glm(counts, mm=mm2, column=3)
#    test.output(conf.out, nrow(counts))
#    expect_true(mean(conf.out$pvalue) < .7 & mean(conf.out$pvalue) > .3)
#})

#test_that("edgeR handler works", {
#    test.output(subSeq::edgeR(testcounts, treatment), nrow(testcounts))
#})

#test_that("DESeq2 handler works", {
#    test.output(DESeq2(testcounts, treatment), nrow(testcounts))
#})

#test_that("voomLimma handler works", {
#    test.output(voomLimma(testcounts, treatment), nrow(testcounts))
#})

# DEXSeq has problems; leaving it out

# test_that("DEXSeq handler works", {
#     data( "pasillaExons", package="pasilla" )
#     n = 200
#     exon.counts = head(pasillaExons@assayData$counts, n)
#     design = pasillaExons@phenoData@data[c("condition", "type")]
#     geneIDs = head(pasillaExons@featureData@data$geneID, n)
#     exonIDs = head(pasillaExons@featureData@data$exonID, n)
# 
#     ret.exon = DEXSeq(exon.counts, design, geneIDs, exonIDs)
#     test.output(ret.exon, n)
# })

context("Subsampling")

ss = subsample(counts, treatments= hammer.design$protocol, proportions,
               bioReplicates= biological.replicate, replications= 1,
               replacement= FALSE, method=c("edgeR", "voomLimma"))
ss.summ = summary(ss)

test_that("subsamples returns a data table with the right structure", {
    expect_that(ss, is_a("subsamples"))
    expect_that(ss, is_a("data.table"))

    # check that the proportions and depths have the properties we expect
    expect_equal(sort(unique(ss$proportion)), proportions)
    #this should be true only when all libraries do not have the same expected size
    expect_equal(max(ss$depth), sum(counts))
    expect_false(is.null(getSeed(ss)))
    
    # check the proportions are about what you expect (has some noise)
    proportion.changes = log2(ss$depth / max(ss$depth) / ss$proportion)
    expect_that(all(proportion.changes < .05), is_true())

    # check that the per-gene counts add up to the total depth
    countsums = ss[, list(total=sum(count)), by=c("depth", "method")]
    expect_equal(countsums$depth, countsums$total)
    
    # check that no replication was performed
    expect_true(all(ss$replication == 1))
})

test_that("quality metrics improve with increasing depth", {
    for (m in unique(ss.summ$method)) {
        sm = ss.summ[method == m]
        expect_equal(sm$proportion, proportions)
        expect_that(all(diff(sm$pearson) > 0), is_true())
        expect_that(all(diff(sm$percent) > 0), is_true())
        expect_that(all(diff(sm$MSE) < 0), is_true())
    }
})

test_that("summaries can be created with other p-value corrections", {
    ss.summ.BH = summary(ss, p.adjust.method="BH")
    ss.summ.bon = summary(ss, p.adjust.method="bonferroni")
    expect_that(all(ss.summ.BH$significant < ss.summ$significant), is_true())
    expect_that(all(ss.summ.bon$significant < ss.summ.BH$significant), is_true())
    
    # check that you can't give it a nonsense method
    expect_that(summary(ss, p.adjust.method="nomethod"), throws_error("should be one of"))
})

ss.rep = subsample(counts, proportions=  c(.1, 1), treatments=treatment, bioReplicates = biological.replicate,
                   method=c("edgeR", "voomLimma"), replications=2)

test_that("Replications (multiple at each proportion) works", {
    sumrep = summary(ss.rep)
    expect_equal(nrow(sumrep), 6)
    rep2 = sumrep[replication == 2]
    expect_equal(nrow(rep2), 2)
    expect_true(all(rep2$proportion == .1))

    # confirm there aren't replicates of 1.0
    expect_false(any(ss.rep$proportion == 1 & ss.rep$replication == 2))
    
    # confirm that averaging works
    av.rep = summary(ss.rep, average=TRUE)
    expect_equal(nrow(av.rep), 4)
    expect_null(av.rep$replication)
})

test_that("subSeq can handle low counts", {
    low.counts = hammer.counts[rowSums(hammer.counts) < 2000, ]
    low.counts = low.counts[sample(nrow(low.counts), 500), ]
    
    low.proportions = c(.01, .1, 1)
    ss.low = subsample(low.counts, treatment, low.proportions, 2,
                       method=c("edgeR", "voomLimma"))
    
    # test that plots still work
    expect_that(plot(summary(ss.low)), is_a("ggplot"))
    
    # significance might get wonky at this level, but correlations had better line up
    summ.low = summary(ss.low)
    for (m in unique(summ.low$method)) {
        sm.low = summ.low[method == m]
        expect_equal(sm.low$proportion, low.proportions)
        expect_that(all(diff(sm.low$pearson) > 0), is_true())
        expect_that(all(diff(sm.low$MSE) < 0), is_true())
    }
    
    # there should be no NAs or infinities
  #  expect_not_na(summ.low)
    expect_false(any(apply(summ.low, 2, is.nan)))    
    expect_false(any(apply(summ.low, 2, is.infinite)))    
})

test_that("Combining subsamples works", {
    seed = getSeed(ss)
    # try three other proportions
    more.proportions = c(.05, .3, .5, 1)
    ss2 = subsample(counts, treatment, more.proportions, biological.replicate,
                         method=c("edgeR", "voomLimma"), replacement = FALSE, seed=seed)
    
    combined = combineSubsamples(ss, ss2)
    expect_equal(getSeed(ss), getSeed(ss2))
})


context("Custom Handlers")

test_that("Can provide custom error handlers", {
    fake.pvalues = runif(nrow(counts))
    custom = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=-.2))
    }
    test.output(custom(counts, treatment), length(fake.pvalues))
    ss.custom = subsample(counts, treatment, proportions, biological.replicate,
                          method=custom, replacement = FALSE)
    
    expect_true(all(ss.custom$method == "custom"))
    expect_equal(fake.pvalues, ss.custom[depth == min(ss.custom$depth)]$pvalue)
    
    # check it can be given as a string as well
    ss.custom2 = subsample(counts, treatment, proportions, bioReplicates = biological.replicate,
                           method="custom", replacement= TRUE)
    expect_true(all(ss.custom2$method == "custom"))
    expect_equal(fake.pvalues, ss.custom2[depth == min(ss.custom2$depth)]$pvalue)
})

test_that("Handlers can have columns that others don't", {
    fake.pvalues = runif(nrow(counts))
    othercols = replicate(3, runif(nrow(counts)))
    custom1 = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=othercols[, 3], other=othercols[, 1]))
    }
    custom2 = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=othercols[, 2], other=othercols[, 2]))
    }
    custom3 = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=othercols[, 1], other3=othercols[, 3]))
    }

    ss.custom = subsample(counts, treatment, proportions, biological.replicate,
                          method=c("edgeR", "custom1", "custom2", "custom3"))

    # we expect it to fill in missing values with NAs
    expect_true(all(c("other", "other3") %in% colnames(ss.custom)))
    expect_na(ss.custom[method == "edgeR", other]) #this is currently not working
    expect_na(ss.custom[method == "custom1", other3])
    expect_na(ss.custom[method == "custom2", other3])
    expect_na(ss.custom[method == "custom3", other])
    for (p in proportions) {
        expect_equal(ss.custom[method == "custom1" & proportion == p, other], othercols[, 1])
        expect_equal(ss.custom[method == "custom2" & proportion == p, other], othercols[, 2])
        expect_equal(ss.custom[method == "custom3" & proportion == p, other3], othercols[, 3])
    }
})

test_that("Handlers don't have to return one row per gene", {
    # some handlers, such as for gene set enrichment, don't necessarily return
    # one row per gene. Check it can return more and less
    for (n in c(nrow(counts) / 2, nrow(counts) * 2)) {
        fake.pvalues = runif(n)
        othercol = runif(n)
        coefs = rnorm(n)
        custom.different = function(counts, treatment) {
            return(data.frame(pvalue=fake.pvalues, coefficient=coefs, other=othercol,
                              ID=as.character(1:n)))
        }
        
        ss.custom = subsample(counts, treatment, proportions, biological.replicate, method=c("edgeR", "custom.different"))
        ss.edgeR = ss.custom[method == "edgeR"]
        ss.custom.different = ss.custom[method == "custom.different"]
        
        expect_equal(nrow(ss.edgeR), nrow(counts) * length(proportions))
        expect_equal(nrow(ss.custom.different), length(coefs) * length(proportions))

        expect_not_na(ss.edgeR$ID)
        expect_not_na(ss.custom.different$ID)
        expect_not_na(ss.custom.different$other)
        expect_na(ss.edgeR$other)
        expect_na(ss.custom.different$count)
        
        expect_not_na(ss.custom[method == "custom.different", other])
        
        # check that it doesn't affect the summary of the data
        custom.summ = summary(ss.custom)
        expect_not_na(custom.summ$pearson)
        expect_not_na(custom.summ$spearman)
    }
    
    # check that the handler does need to return an ID column in that case
    custom.noID = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=coefs, other=othercol))
    }
    expect_that(subsample(counts, treatments=treatment, proportions,
                          bioReplicates = biological.replicate, method="custom.noID"),
                throws_error("if a handler doesn't return one row per gene then it must specify an ID collumn"))
})


context("Reproducibility from seeds")

s.edgeR = ss[method == "edgeR"]
s.voom = ss[method == "voomLimma"]

test_that("seeds are reproducible between methods", {
    summ.edgeR = ss.summ[method == "edgeR"]
    summ.voom = ss.summ[method == "voomLimma"]

    expect_equal(s.edgeR$depth, s.voom$depth)
    expect_equal(s.edgeR$count, s.voom$count)
})

test_that("seeds are reproducible between runs", {
    # perform it again with the same seed, and see that it matches the
    # first replication all the way through
    ss2 = subsample(counts, proportions, treatments=treatment,
                    bioReplicates = biological.replicate,
                         method=c("edgeR", "voomLimma"),
                         seed=getSeed(ss))
    
    expect_equal(ss$count, ss2$count)
    expect_equal(ss$pvalue, ss2$pvalue)
    expect_equal(ss$coefficient, ss2$coefficient)

    # for sanity check, check the same with the summary
    ss.summ2 = summary(ss2)
    expect_equal(ss.summ$significant, ss.summ2$significant)
    expect_equal(ss.summ$MSE, ss.summ2$MSE)
})

test_that("generateSubsampledMatrix retrieves the correct subsampled matrices", {
    treatments= hammer.design$protocol
    p = min(ss$proportion)
    seed = getSeed(ss)
    indeces <- getIndecesFromCatagoricalTreatment( treatments, rep= 2, seed, replacement= FALSE)
    subm = generateSubsampledMatrix(counts, indeces, rep_len(p, length(indeces)), seed)
    expect_that(subm, is_a("matrix"))
    expect_equal(dim(subm), dim(counts))
    expect_equal(sum(subm), min(ss$depth))

    expect_equal(rowSums(subm), s.edgeR[proportion == p]$count, check.names = FALSE)
    expect_equal(rowSums(subm), s.voom[proportion == p]$count, check.names = FALSE)

    # confirm that edgeR on that matrix gives the same results
    subm.results = edgeR(subm, treatment=treatment)
    expect_equal(subm.results$pvalue, s.edgeR[proportion == p]$pvalue)
    expect_equal(subm.results$coefficient, s.edgeR[proportion == p]$coefficient)
    
    # confirm summary object also contains correct seed
    summ = summary(ss)
    expect_equal(subm, generateSubsampledMatrix(counts, indeces, rep_len(p, length(indeces)), getSeed(summ)))
    
    # confirm generateSubsampledMatrix works if you explicitly
    # tell it replication is 1
    subm.rep1 = generateSubsampledMatrix(counts, indeces, rep_len(p, length(indeces)), seed, replication=1)
    expect_equal(subm, subm.rep1)
})

test_that("Performing multiple replicates is reproducible", {
    ss.rep.2 = subsample(counts, treatment=treatment, c(.1, 1), biological.replicate,
                    method=c("edgeR", "voomLimma"), replications=2,
                    seed=getSeed(ss.rep))
    expect_equal(ss.rep$pvalue, ss.rep.2$pvalue)
    
})

context("Plotting")

test_that("plotting is possible without errors", {
    expect_that(plot(ss.summ), is_a("ggplot"))
})


context("Error handling")

test_that("Raises an error on edgeR if there are >2 treatments", {
    new.treatment = c("A", "A", "B", "C")
    # check with multiple handlers
    expect_that(subsample(counts, treatments=new.treatment, proportions, biological.replicate,
                          method="edgeR"), throws_error("more than two levels"))
})

test_that("Raises an error if it cannot find the handler", {
    expect_that(subsample(counts, treatment=treatment, proportions, biological.replicate,
                          method="nomethod"),
                throws_error("Could not find handler nomethod"))
})

test_that("error messages are thrown when proportions are incorrect", {    
    expect_that(subsample(counts, treatments=treatment, c(), biological.replicate, method="edgeR"),
                throws_error("No proportions"))
    expect_that(subsample(counts, treatments=treatment, c(.1, 1, 2), biological.replicate, method="edgeR"),
                throws_error("Proportions must be in range"))
    expect_that(subsample(counts, treatment=treatment, c(0, 1), biological.replicate, method="edgeR"),
                throws_error("Proportions must be in range"))
})

test_that("error message is thrown if counts were normalized", {
    # confirming that there was no normalization
    expect_that(subsample(scale(counts), treatments=treatment, proportions, bioReplicates = biological.replicate,
                          method="edgeR"), throws_error("unnormalized"))
    sc.counts = scale(counts, center=FALSE)
    expect_that(subsample(sc.counts, treatment=treatment, proportions, bioReplicates = biological.replicate,
                          method="edgeR"), throws_error("unnormalized"))
})

test_that("combineSubsamples raises an error when combining different seeds", {
    more.proportions = c(.05, .3, .5)
    ss2 = subsample(counts, treatment=treatment, more.proportions, bioReplicates = biological.replicate,
                    method=c("edgeR", "voomLimma"))

    expect_false(getSeed(ss) == getSeed(ss2))
    expect_that(combineSubsamples(ss, ss2), throws_error("different.*seed"))
})

# handlers:
# -write baySeq, DEXSeq

# plots:
# -volcano plots
# -per-gene plots
