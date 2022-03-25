mi.eval <- function(EXPR, robust, cluster, coef., vcov., df.=NULL, parallel=FALSE, lazy=NULL, ...) {
    mf.raw <- match.call()
    if (!("EXPR" %in% names(mf.raw)))
        stop("Must specify R command to apply across multiply imputed datasets.")
    if (!("robust" %in% names(mf.raw)))
        robust <- ("cluster" %in% names(mf.raw))
    if (is.null(robust))
        robust <- ("cluster" %in% names(mf.raw))
    if ((!robust) & ("cluster" %in% names(mf.raw)))
        stop("Cannot use \"cluster\" option with \"robust=FALSE\".")
    coef.expr <- mf.raw[match(c("coef.","EXPR"), names(mf.raw), 1L)]
    if (!("coef." %in% names(mf.raw)))
        coef.expr[[1L]] <- quote(stats::coef)
    m <- c(match(c("vcov.", "EXPR"), names(mf.raw), 1L), match("cluster", names(mf.raw), 0L))
    vcov.expr <- mf.raw[m]
    mthd.plus <- NULL
    if (!("vcov." %in% names(mf.raw))) {
        if ("cluster" %in% names(vcov.expr)) {
            vcov.expr[[1L]] <- quote(bucky::vcovCR)
                mthd.plus <- as.character(enquote(mf.raw$cluster)[-1])
            if (nchar(mthd.plus) >= 40)
                mthd.plus <- c("robust standard errors clustered on:",mthd.plus)
            else
                mthd.plus <- paste("robust standard errors clustered on",mthd.plus)
        }
        else if (robust) {
            vcov.expr[[1L]] <- quote(sandwich::vcovHC)
            mthd.plus <- "robust standard errors"
        }
        else
            vcov.expr[[1L]] <- quote(stats::vcov)
    }
    else {
        mthd.plus <- c("variance-covariance matrix computed using:", deparse(vcov.expr))
    }
    coef.expr[[2]] <- quote(imp.list[[i]])
    vcov.expr[[2]] <- 
    names(coef.expr)[2] <- names(vcov.expr)[2] <- ""
    mf <- mf.raw[["EXPR"]]
    m <- match(c("data"), names(mf), 0L)
    if (m==0)
        stop("Must use \"data\" option for multiple imputation.")
    if (is.null(lazy))
        lazy <- is.name(mf[[m]])
    if (!lazy) {
        imputations <- eval(mf[[m]])
        mf[[m]] <- as.name("imputations")
    }
    mfmclass <- eval(substitute(class(x), list(x=mf[[m]])))
    if ("mids" %in% mfmclass)
        imp.info <- list(size=eval(substitute((x)$m, list(x=mf[[m]]))),
                    names=NULL,
                    sub=function(i, mf, m) { mf[[m]] <- substitute(mice::complete(x, y), list(x=mf[[m]], y=i)); return(eval(mf)) })
    else if ("amelia" %in% mfmclass || "imputationList" %in% mfmclass)
        imp.info <- list(size=eval(substitute(length((x)$imputations), list(x=mf[[m]]))),
                    names=eval(substitute(names((x)$imputations), list(x=mf[[m]]))),
                    sub=function(i, mf, m) { mf[[m]] <- substitute((x)$imputations[[y]], list(x=mf[[m]], y=i)); return(eval(mf)) })
    else if ("data.frame" %in% mfmclass)
        imp.info <- list(size=1L,
                    names=NULL,
                    sub=function(i, mf, m) { return(eval(mf)) })
    else if (eval(substitute(is.list(x) && is.data.frame((x)[[1L]]), list(x=mf[[m]]))))
        imp.info <- list(size=eval(substitute(length(x), list(x=mf[[m]]))),
                         names=eval(substitute(names(x), list(x=mf[[m]]))),
                         sub=function(i, mf, m) { mf[[m]] <- substitute((x)[[y]], list(x=mf[[m]], y=i)); return(eval(mf)) })
    else
        stop("No imputations found.")
    num.imp <- imp.info$size
    if (num.imp < 1L)
        stop("Empty list of imputations.")
    if (num.imp == 1L) {
        warning("Only one data set found. Returning regular model estimates.")
        return(eval(mf))
    }
    if (is.null(parallel)) {
        parallel <- requireNamespace("parallel", quietly=TRUE)
        if (parallel)
            parallel <- getOption("mc.cores", parallel::detectCores() - 1L) > 1L
        if (is.na(parallel))
            parallel <- FALSE
    }
    if (parallel) {
        parallel <- requireNamespace("parallel", quietly=TRUE)
        if (!parallel)
            warning("Can't load \"parallel\" parallel. Using serial computation.")
    }
    if (parallel) {
        imp.list <- parallel::mclapply(1:num.imp, imp.info$sub, mf=mf, m=m, ..., mc.silent=TRUE, mc.allow.recursive=FALSE)
        if (any(sapply(imp.list, inherits, what="try-error"))) {
            stop(imp.list[which(sapply(imp.list, inherits, what="try-error"))[1]],
                 call.=FALSE)
        }
    } else {
        imp.list <- lapply(1:num.imp, imp.info$sub, mf=mf, m=m, ...)
    }
    names(imp.list) <- imp.info$names
    if (is.null(df.)) {
        if (all(sapply(imp.list, inherits, what="glm") | sapply(imp.list, inherits, what="glmerMod")))
            df. <- Inf
        else
            df. <- df.residual
    }
    if (is.function(df.)) {
        df. <- try(min(sapply(imp.list, df.)), silent=TRUE)
        if (inherits(df., "try-error") || !is.numeric(df.))
            df. <- Inf
    }
    if (is.null(df.))
        df. <- Inf
    coeflist <- suppressWarnings(lapply(1:length(imp.list), function(i, expr) { expr[[2]] <- substitute(imp.list[[i]], list(i=i)); eval(expr) }, expr=coef.expr))
    vcovlist <- suppressWarnings(lapply(1:length(imp.list), function(i, expr) { expr[[2]] <- substitute(imp.list[[i]], list(i=i)); eval(expr) }, expr=vcov.expr))
    vcovnames <- rownames(vcovlist[[1]])
    coefnames <- names(coeflist[[1]])
    k <- min(length(coeflist[[1]]),nrow(vcovlist[[1]]))
    ind.c <- 1:k
    ind.v <- 1:k
    if ((length(names(coeflist[[1]])) > 0) & all(names(coeflist[[1]]) %in% rownames(vcovlist[[1]])))
        ind.v <- match(names(coeflist[[1]]), rownames(vcovlist[[1]]))
    else if ((length(rownames(vcovlist[[1]])) > 0) & all(rownames(vcovlist[[1]]) %in% names(coeflist[[1]])))
        ind.c <- match(rownames(vcovlist[[1]]),names(coeflist[[1]]))
    coeflist <- vapply(coeflist, function(x, ind) return(as.vector(x)[ind]), rep(0,k), ind=ind.c)
    vcovlist <- vapply(vcovlist, function(x, ind) return(as.matrix(x)[ind,ind,drop=FALSE]), diag(k), ind=ind.v)
    if (is.null(coefnames))
        coefnames <- vcovnames[ind.v]
    else
        coefnames <- coefnames[ind.c]
    if (k==1) {
        vcov.within <- as.matrix(mean(vcovlist))
        coef.mean <- mean(coeflist)
        vcov.between <- as.matrix(var(as.vector(coeflist)))
    }
    else {
        vcov.within <- suppressWarnings(apply(vcovlist,1:2,mean))
        coef.mean <- apply(coeflist,1,mean)
        vcov.between <- var(t(coeflist))
    }
    names(coef.mean) <- rownames(vcov.between) <- colnames(vcov.between) <- coefnames
    df.model <- (num.imp - 1) * (1 + (num.imp/(num.imp+1)) * diag(vcov.within) / diag(vcov.between))^2
    if (any(is.finite(df.))) {
        dfobs <- ((1 - 2/(df. + 3)) * df.) * (diag(vcov.within)/(diag(vcov.within) + diag(vcov.between)))
        df.model <- 1/(1/dfobs + 1/df.model)
    }
    m.out <- list(coefficients = coef.mean, vcov = vcov.within + vcov.between * (1 + 1/num.imp),
                  num.imp = num.imp, df.residual = df.model)
    imp.nobs <- try(sapply(imp.list, nobs), silent=TRUE)
    if (!inherits(imp.nobs, "try-error")) {        
        if (length(unique(imp.nobs)) == 1) {
            m.out$nobs <- unique(imp.nobs)
            names(m.out$nobs) <- NULL
        }
        else if (!all(is.na(imp.nobs))) {
            warning("Models run on data sets of different sizes.")
            m.out$nobs <- imp.nobs
        }
    }
    imp.family <- try(lapply(imp.list, family), silent=TRUE)
    if (!inherits(imp.family, "try-error") && length(unique(imp.family)) == 1)
        m.out$family <- imp.family[[1L]]
    imp.class <- try(lapply(imp.list, class), silent=TRUE)
    if (!inherits(imp.class, "try-error") && length(unique(imp.class)) == 1)
        m.out$model.class <- imp.class[[1L]]
    imp.terms <- try(lapply(imp.list, terms), silent=TRUE)
    if (!inherits(imp.terms, "try-error")) {
        for (i in 1:length(imp.terms))
            attr(imp.terms[[i]], ".Environment") <- NULL
        if (length(unique(imp.terms)) == 1)
            m.out$terms <- imp.terms[[1L]]
    }
    imp.resid <- try(sapply(imp.list, function(x) x$residuals), silent=TRUE)
    if (!inherits(imp.resid, "try-error")) {
        if (is.matrix(imp.resid))
            m.out$residuals <- apply(imp.resid, 1, mean, na.rm=TRUE)
        else if (is.numeric(imp.resid))
            m.out$residuals <- imp.resid
        else if (is.list(imp.resid) && all(sapply(imp.resid, is.numeric)))
            m.out$residuals <- rep(NA, min(sapply(imp.resid, length)))
    }
    if (any(is.finite(df.)))
        m.out$df.complete <- df.
    if ("missMatrix" %in% eval(substitute(names(x), list(x=mf[[m]]))))
        m.out$miss <- eval(substitute(apply((x)$missMatrix,2,mean), list(x=mf[[m]])))
    if (!is.null(mthd.plus))
        m.out$note <- mthd.plus
    class(m.out) <- "mi.estimates"
    m.out$mi.call <- mf.raw
    m.out$call <- mf
    return(m.out)
}

vcov.mi.estimates <- function(object, ...) return(object$vcov)
print.mi.estimates <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n", paste(deparse(x$mi.call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}
summary.mi.estimates <- function(object, ...) {
    m.out <- list(coefficients = coeftest(object, ...), df.residual=object$df.residual, mi.call=object$mi.call, call=object$call, num.imp=object$num.imp, note=object$note)
    attr(m.out$coefficients, "method") <- paste("Coefficients")
    class(m.out) <- "summary.mi.estimates"
    return(m.out)
}
print.summary.mi.estimates <- function(x, ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Coefficients:\n")
    printCoefmat(x$coefficients, ...)
    if (is.null(x$note))
        cat("\nComputed using multiple imputation with",x$num.imp,"imputed data sets.")
    else
        cat("\nComputed using multiple imputation with",x$num.imp,"imputed data sets\nand",paste(x$note,collapse="\n"))
    cat("\n\n")
    invisible(x)
}
extract.mi.estimates <- function(model, include.nobs = TRUE, include.imp = TRUE) {
    s <- summary(model)
    gof <- list()
    if (include.nobs)
        gof[["Num. obs."]] <- try(nobs(model))
    if (include.imp)
        gof[["Imputations"]] <- try(as.integer(model$num.imp))
    gof <- gof[!(sapply(gof, inherits, what="try-error") | sapply(gof, is.null))]
    return(texreg::createTexreg(
        coef.names = rownames(s$coef),
        coef = s$coef[, 1],
        se = s$coef[, 2],
        pvalues = s$coef[, 4],
        gof.names = names(gof),
        gof = unlist(gof),
        gof.decimal = sapply(gof, is.numeric) & !sapply(gof, is.integer)
    ))
}

predict.mi.estimates <- function (object, newdata = NULL, type = c("link", "response", "terms"), se.fit = FALSE, terms = NULL, na.action = na.pass, ...) {
    type <- match.arg(type)    
    if (missing(newdata))
        stop("The predict method for mi.estimates objects is only implemented when newdata is specified")
    if (type == "terms")
        stop("The predict method for mi.estimates objects is not yet implemented for type=\"terms\"")
    if (!any(c("glm","lm") %in% object$model.class))
        stop("The predict method for mi.estimates objects is only implemented for those based on glm and lm objects")
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, newdata, na.action = na.action)
    X <- model.matrix(Terms, m)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(object$terms, "offset"))) 
        for (i in off.num)
            offset <- offset + eval(attr(object$terms, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
        offset <- offset + eval(object$call$offset, newdata)
    predictor <- as.vector(X %*% coef(object))
    if (!is.null(offset)) 
        predictor <- predictor + offset
    if (type == "response" && !is.null(object$family))
        predictor <- object$family$linkinv(predictor)
    names(predictor) <- rownames(X)
    if (se.fit) {
        predictor.se <- as.vector(sqrt(apply(X * (X %*% vcov(object)), 1, sum)))
        if (type == "response" && !is.null(object$family))
            predictor.se <- predictor.se * abs(object$family$mu.eta(predictor))
        names(predictor.se) <- rownames(X)
        return(list(fit = predictor, se.fit = predictor.se))
    } else {
        return(predictor)
    }
}
