##' Database sampling tools
##'
##' Tools for designing database sampling schemes.  This can be done with
##' \code{gsDesign}, but this package is notably faster for large
##' databases.
##' @docType package
##' @useDynLib dbSamplr
##' @importFrom Rcpp sourceCpp
##' @name dbSamplr
NULL

## unload the shared object
.onUnload <- function(libpath) {
    library.dynam.unload("dbSamplr", libpath)
}
