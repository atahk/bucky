.onLoad <- function(libname, pkgname) {
    if (suppressWarnings(requireNamespace("texreg", quietly=TRUE))) {
        setGeneric("extract", function(model, ...) standardGeneric("extract"),
                   package = "texreg")
        setMethod("extract",
                  signature = className("robustified", pkgname),
                  definition = extract.robustified)
        setMethod("extract",
                  signature = className("mi.estimates", pkgname),
                  definition = extract.mi.estimates)
    }
    invisible()
}
