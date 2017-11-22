mi.eval <- function(EXPR, robust, cluster, coef., vcov., df.=NULL, parallel=FALSE, ...) {
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
                mthd.plus <- as.character(enquote(mf$cluster)[-1])
            if (nchar(mthd.plus) >= 40)
                mthd.plus <- c("robust standard errors clustered on:",mthd.plus)
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
    if (is.null(eval(mf[[m]])$imputations))
        stop("No imputations found.")
    num.imp <- length(eval(mf[[m]])$imputations)
    if (num.imp < 1)
        stop("Empty list of imputations.")
    if (num.imp == 1) {
        mf[[m]] <- substitute((x)$imputations[[y]], list(x=mf[[m]], y=1L))
        warning("Only one data set found. Returning regular model estimates.")
        return(eval(mf))
    }
    if (parallel) {
        parallel <- requireNamespace("parallel", quietly=TRUE)
        if (!parallel)
            warning("Can't load \"parallel\" parallel. Using serial computation.")
    }
    if (parallel)
        imp.list <- parallel::mclapply(1:num.imp, function(i, mf) { mf[[m]] <- substitute((x)$imputations[[y]], list(x=mf[[m]], y=i)); return(eval(mf)) }, mf=mf, ...)
    else
        imp.list <- lapply(1:num.imp, function(i, mf) { mf[[m]] <- substitute((x)$imputations[[y]], list(x=mf[[m]], y=i)); return(eval(mf)) }, mf=mf, ...)
    names(imp.list) <- names(eval(mf[[m]])$imputations)
    if (is.null(df.)){
        if (all(sapply(imp.list, inherits, what="glm")))
            df. <- Inf
        else
            df. <- df.residual
    }
    if (is.function(df.))
        df. <- min(sapply(imp.list, df.))
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
    else if ((length(rownames(vcovlist[[1]])) > 0) & all(rownames(vcovlist[[1]])) %in% names(coeflist[[1]]))
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
    if (any(is.finite(df.)))
        m.out$df.complete <- df.
    if ("missMatrix" %in% eval(substitute(names(x), list(x=mf[[m]]))))
        m.out$miss <- eval(substitute(apply((x)$missMatrix,2,mean), list(x=mf[[m]])))
    if (!is.null(mthd.plus))
        m.out$note <- mthd.plus
    class(m.out) <- "mi.estimates"
    m.out$call <- mf.raw
    m.out$subcall <- mf
    return(m.out)
}

vcov.mi.estimates <- function(object, ...) return(object$vcov)
print.mi.estimates <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
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
    m.out <- list(coefficients = coeftest(object, ...), df.residual=object$df.residual, call=object$call, subcall=object$subcall, num.imp=object$num.imp, note=object$note)
    attr(m.out$coefficients, "method") <- paste("Coefficients")
    class(m.out) <- "summary.mi.estimates"
    return(m.out)
}
print.summary.mi.estimates <- function(x, ...) {
    cat("\nCall:\n", paste(deparse(x$subcall), sep = "\n", collapse = "\n"), 
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
extract.mi.estimates <- function(model, include.nobs = TRUE) {
    s <- summary(model)
    gof <- list()
    if (include.nobs)
        gof[["Num. obs."]] <- try(nobs(model))
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
