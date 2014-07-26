#' plot metrics as a function of subsampled read depth
#' 
#' Plot the number of genes found significant, the Spearman correlation of the
#' effect size estimates with the full experiment, and the empirical false
#' discovery rate as a function of the subsampled read depth. This determines
#' whether these metrics saturate, which indicates that the experiment has an
#' appropriate sequencing depth.
#' 
#' @param x a \code{summary.subsamples} object
#' @param ... further arguments passed to or from other methods.
#' 
#' @import ggplot2
#' @import reshape2
#' 
#' @export
plot.summary.subsamples <-
function(x, ...) {
    vars = c("significant", "estFDP", "spearman", "MSE")
    id.vars = c("method", "depth", "replication", "proportion")
    metrics = as.data.table(melt(as.data.frame(x), id=id.vars))

    metrics = metrics[variable %in% vars]

    # change order and appearance of levels
    metrics$variable = factor(metrics$variable, levels=vars)
    levels(metrics$variable) = c("# significant genes at 5% FDR", "Estimated FDP",
                                 "Spearman corr of estimates", "MSE of estimates")

    # average within replications
    metrics[, average.depth:=mean(depth, na.rm=TRUE), by=c("method", "proportion", "variable")]
    metrics[, average.value:=mean(value, na.rm=TRUE), by=c("method", "proportion", "variable")]
    
    g = (ggplot(metrics, aes(col=method)) + geom_line(aes(x=average.depth, y=average.value)) +
             facet_wrap("variable", nrow=2, scale="free_y") +
             theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="top") +
             xlab("Depth") + ylab("Metric"))
    if (any(metrics$replication > 1)) {
        g = g + geom_point(aes(x=depth, y=value))
    }
    return(g)
}

## helper functions

subsample.plot.helper <-
    function(summary.ss, stat="significant") {
        (ggplot(summary.ss, aes_string(x="depth", y=stat, col="method")) +
             geom_line())
    }

g_legend <-
    function(a.gplot) {
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }
