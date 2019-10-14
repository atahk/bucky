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
        if (!("type" %in% names(mf))) {
            if (inherits(x, "glm"))
                mf[["type"]] <- "CR"
            else if (inherits(x, "lm"))
                mf[["type"]] <- "CR1"
            else
                mf[["type"]] <- "CR"
        }
        mthd.plus <- paste("robust standard errors clustered on", 
            as.character(enquote(mf$cluster)[-1]))
    }
    else {
        m <- match(c("x", "type", "omega"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1]] <- quote(vcovHC)
        if (!("type" %in% names(mf))) {
            if (inherits(x, "glm"))
                mf[["type"]] <- "HC0"
            else if (inherits(x, "lm"))
                mf[["type"]] <- "HC1"
            else
                mf[["type"]] <- "HC0"
        }
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

vcovHC.polr <- function (x, type = c("HC1", "HC0", "HC"), omega = NULL, sandwich = TRUE, ...)
{
    type <- match.arg(type)
    if (type == "HC") 
        type <- "HC0"
    if (is.null(omega)) {
        switch(type,
               HC0 = omega <- function(residuals, diaghat, df) 1/length(residuals),
               HC1 = omega <- function(residuals, diaghat, df) 1/df)
    }
    estf <- estfun(x)
    out <- crossprod(estf, omega(rep(1,nrow(estf)), NULL, df.residual(x)) * estf)
    if (sandwich) 
        out <- sandwich(x, meat. = out, ...)
    return(out)
}

predict.robustified <- function(object, newdata = NULL, se.fit = FALSE,
                                interval = c("none", "confidence", "prediction"),
                                level = 0.95,
				na.action = na.pass,
                                ...) {
    interval <- match.arg(interval)
    if ((interval == "none") && !se.fit) {
        return(NextMethod())
    } else {
        if (inherits(object, "lm")) {
            call <- match.call()
            if ("type" %in% names(call))
                type <- eval(call[[which("type" == names(call))]])
            else if (inherits(object, "glm"))
                type <- "link"
            else
                type <- "response"
            if (!(type %in% c("link", "response", "terms")))
                stop("invalid type argument")
            else if (type == "terms")
                stop("type \"terms\" not yet implemented")
            if (missing(newdata)) {
                X <- model.matrix(object)
                offset <- object$offset
            } else {
                tt <- terms(object)
                Terms <- delete.response(tt)
                m <- model.frame(Terms, newdata, na.action = na.action,
                                 xlev = object$xlevels)
                if (!is.null(cl <- attr(Terms, "dataClasses"))) 
                    .checkMFClasses(cl, m)
                X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
                offset <- rep(0, nrow(X))
                if (!is.null(off.num <- attr(tt, "offset"))) 
                    for (i in off.num)
                        offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
                if (!is.null(object$call$offset)) 
                    offset <- offset + eval(object$call$offset, newdata)
            }
            beta <- object$coefficients
            fit <- as.vector(X %*% beta)
            if (!is.null(offset)) 
                fit <- fit + as.vector(offset)
            stderr <- sqrt(apply((X %*% object$cov.robust) * X, 1, sum))
            switch(interval, confidence = {
                fit <- outer(stderr, qnorm(0.5 + level * c(fit=0, lwr=-0.5, upr=0.5))) + fit
            }, prediction = {
                stop("prediction confidence intervals not available for robustified objects")
            })
            switch(type, response = {
                stderr <- stderr * abs(family(object)$mu.eta(as.matrix(fit)[,1]))
                fit <- family(object)$linkinv(fit)
            }, link = , terms = )
            if (is.matrix(fit) && is.null(rownames(fit)) && !is.null(names(stderr)))
                rownames(fit) <- names(stderr)
            else if ((!is.matrix(fit)) && is.null(names(fit)) && !is.null(names(stderr)))
                names(fit) <- names(stderr)
            if (se.fit)
                return(list(fit = fit, se.fit = stderr))
            else
                return(fit)
        }
        if (interval != "none")
            stop("confidence intervals not yet implemented for robustified objects of this type")
        else
            stop("standard errors not yet implemented for robustified objects of this type")
    }
}

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
