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
#' @importFrom dplyr filter group_by mutate
#' @importFrom tidyr gather
#' 
#' @export
plot.summary.subsamples <-
function(x, ...) {
    vars = c("significant", "estFDP", "spearman", "MSE")
    
    FDR.level = attr(x, "FDR.level")
    if (is.null(FDR.level) || FDR.level == .05) {
        FDR.level = "5%"
    }
    
    metrics <- as.data.frame(x) %>% gather(metric, value, significant:percent) %>%
                    filter(metric %in% vars)

    # change order and appearance of levels
    metrics$metric = factor(metrics$metric, levels=vars)
    levels(metrics$metric) = c(paste("# significant genes at", FDR.level, "FDR"), "Estimated FDP",
                                 "Spearman corr of estimates", "MSE of estimates")

    # average within replications
    metrics <- metrics %>% group_by(method, depth, metric) %>%
        mutate(average.depth=mean(depth, na.rm=TRUE), average.value=mean(value, na.rm=TRUE))
    
    g = (ggplot(metrics, aes(col=method)) + geom_line(aes(x=average.depth, y=average.value)) +
             facet_wrap(~ metric, nrow=2, scales="free_y") +
             theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="top") +
             xlab("Depth") + ylab("Metric"))
    if (any(metrics$replication > 1)) {
        g = g + geom_point(aes(x=depth, y=value))
    }
    return(g)
}

#' plot metrics as a function of subsampled read depth
#' 
#' Plot the number of genes found significant, the Spearman correlation of the
#' effect size estimates with the full experiment, and the empirical false
#' discovery rate as a function of the subsampled read depth. This determines
#' whether these metrics saturate, which indicates that the experiment has an
#' appropriate sequencing depth.
#' 
#' This is an alias for the \link{plot.summary.subsamples} function, so that
#' plotting can be done directly on the subsamples object. We recommend using
#' \code{summary(ss)} first, so that the summary operation does not have to
#' be performed each time the figure is plotted, and so the summary object
#' can be examined on its own.
#' 
#' @param x a \code{subsamples} object
#' @param ... further arguments passed to or from other methods.
#' 
#' @export
plot.subsamples <-
    function(x, ...) {
        summ = summary(x)
        plot(summ)
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
