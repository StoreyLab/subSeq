#' combine multiple subsamples objects
#' 
#' Given two or more subsamples objects, combine them into one larger object, on
#' which we can perform all the usual analyses and plots.
#' 
#' @param ... Two or more subsamples objects
#' 
#' @details If there are columns in some subsamples objects that are not in others,
#' the missing values will be filled with NA
#' 
#' @export
combineSubsamples <-
    function(...) {
        lst = list(...)

        # confirm that the global seed is the same
        seeds = sapply(lst, function(s) attr(s, "seed"))
        if (length(unique(seeds)) != 1) {
            stop(paste("Cannot combine subsamples objects that have",
                 "different random seeds: methods are not comparable"))
        }

        all.colnames = Reduce(union, lapply(lst, colnames))
        with.NA = lapply(lst, function(s) {
            ifnot = lapply(all.colnames, function(i) rep(NA, nrow(s)))
            as.data.frame(mget(all.colnames, envir=as.environment(s),
                                ifnotfound=ifnot))
        })

        ret = as.data.table(do.call(rbind, with.NA))

        class(ret) = c("subsamples", "data.table", "data.frame")
        # save the seed
        attr(ret, "seed") = attr(lst[[1]], "seed")
        ret
    }
