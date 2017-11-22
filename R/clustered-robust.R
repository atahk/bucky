vcovCR <- function(x, cluster = NULL, type = c("CR", "CR0", "CR1")) {
    mf <- x$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$cluster <- match.call()$cluster
    cluster <- eval(mf, parent.frame())$"(cluster)"
    if (is.null(cluster))
        stop("Must specify cluster variable.")
    if (any(is.na(cluster)))
        stop("Cluster cannot have missing values.")
    nclust <- length(unique(cluster))
    type <- match.arg(type)
    switch(type,
           CR = omega <- crossprod(apply(estfun(x), 2, function(x) tapply(x, cluster, sum))) * nclust/(nclust-1),
           CR0 = omega <- crossprod(apply(estfun(x), 2, function(x) tapply(x, cluster, sum))),
           CR1 = omega <- crossprod(apply(estfun(x), 2, function(x) tapply(x, cluster, sum))) * (nclust/(nclust-1))*((nobs(x)-1)/(nobs(x)-x$rank)))
    vcov. <- sandwich(x, meat.=omega/nobs(x))
    attr(vcov., "type") <- type
    return(vcov.)
}

robust.summary <- function(x, cluster, type, omega, ...) {
    mf <- match.call()
    m <- match(c("x", "cluster", "type", "omega"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1]] <- quote(robustify)
    return(summary(eval(mf), ...))
}

robustify <- function (x, cluster, type, omega, ...) 
{
    mf <- match.call()
    if ("cluster" %in% names(mf)) {
        m <- match(c("x", "cluster", "type"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1]] <- quote(vcovCR)
        mthd.plus <- paste("robust standard errors clustered on", 
            as.character(enquote(mf$cluster)[-1]))
    }
    else {
        m <- match(c("x", "type", "omega"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1]] <- quote(vcovHC)
        mthd.plus <- "robust standard errors"
    }
    if (is.list(x)) {
        ro <- x
        class(ro) <- c("robustified", class(x))
    } else {
        ro <- list(coefficients = coef(x))
        class(ro) <- "robustified"
    }
    ro$cov.robust <- eval(mf)
    mf[names(mf) == "x"] <- list(as.name("..."))
    ro$robust.call <- mf
    attr(ro, "method") <- mthd.plus
    return(ro)
}

summary.robustified <- function(object, ...) {
    s <- try(NextMethod())
    if (!is.list(s))
        s <- list()
    s$robust.call <- object$robust.call
    s$cov.scaled <- object$cov.robust
    s$coefficients <- as.matrix(coeftest(x=object, vcov.=s$cov.scaled, ...))
    if ("fstatistic" %in% names(s))
        s <- s[names(s) != "fstatistic"]
    attr(s, "method") <- paste(attr(s$coefficients, "method"), "with", attr(object, "method"))
    attr(s$coefficients, "method") <- NULL
    if (all(class(s) == "list"))
        class(s) <- "summary.robustified"
    else
        class(s) <- c("summary.robustified", class(s))
    return(s)
}

print.robustified <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (!is.null(x$call)) {
        cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
            "\n", sep = "")
    }
    if (!is.null(x$robust.call)) {
        cat("\nCall (vcov):  ", paste(deparse(x$robust.call), sep = "\n", collapse = "\n"), 
            "\n", sep = "")
    }
    cat("\n")
    if (length(coef(x))) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts)) 
            cat("  [contrasts: ", apply(cbind(names(co), co), 
                1L, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
                      print.gap = 2, quote = FALSE)
        cat("\n")
    }
    if (!is.null(x$na.action)) {
        if (nzchar(mess <- naprint(x$na.action))) 
            cat("  (", mess, ")\n\n", sep = "")
    }
    invisible(x)
}

print.summary.robustified <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      signif.stars = getOption("show.signif.stars"), ...) {
    if (!is.null(x$call)) {
        cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
            "\n", sep = "")
    }
    if (!is.null(x$robust.call)) {
        cat("\nCall (vcov):  ", paste(deparse(x$robust.call), sep = "\n", collapse = "\n"), 
            "\n", sep = "")
    }
    cat("", strwrap(paste0(attr(x, "method"), ":")), "", sep="\n")
    printCoefmat(coef(x), digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\n")
    if (!is.null(x$na.action)) {
        if (nzchar(mess <- naprint(x$na.action))) 
            cat("  (", mess, ")\n\n", sep = "")
    }
    invisible(x)
}

vcov.robustified <- function(object, ...) object$cov.robust
vcov.summary.robustified <- function(object, ...) object$cov.scaled

extract.robustified <- function(model, include.aic = TRUE, include.bic = TRUE,
                                include.loglik = TRUE, include.deviance = TRUE,
                                include.rsquared = TRUE, include.adjrs = TRUE,
                                include.nobs = TRUE, include.rmse = FALSE) {
    s <- summary(model)
    gof <- list()
    if (is.list(s)) {
        if (include.rsquared)
            gof[["R$^2$"]] <- try(s$r.squared)
        if (include.adjrs)
            gof[["Adj. R$^2$"]] <- try(s$adj.r.squared)
        if (include.rmse)
            gof[["RMSE"]] <- try(s$sigma[[1]])
    }
    if (include.aic)
        gof[["AIC"]] <- try(AIC(model))
    if (include.bic)
        gof[["BIC"]] <- try(BIC(model))
    if (include.loglik)
        gof[["Log Likelihood"]] <- try(logLik(model))
    if (include.deviance)
        gof[["Deviance"]] <- try(deviance(model))
    if (include.nobs)
        gof[["Num. obs."]] <- try(nobs(model))
    gof <- gof[!(sapply(gof, inherits, what="try-error") | sapply(gof, is.null))]
    return(texreg::createTexreg(
        coef.names = rownames(coef(s)),
        coef = coef(s)[, 1],
        se = coef(s)[, 2],
        pvalues = coef(s)[, 4],
        gof.names = names(gof),
        gof = unlist(gof),
        gof.decimal = sapply(gof, is.numeric) & !sapply(gof, is.integer)
    ))
}
