### Skewness
skewness <- function(x, normal = FALSE) {
    if (!is.null(ncol(x))) ans <- apply(x, 2, skewness, normal)
    else if (normal) {
        x <- x[!is.na(x)]
        ans <- sum((x - mean(x))^3) / ((length(x) - 1) * var(x)^(3/2))
    } else {
        x <- x[!is.na(x)]
        ans <- sum((x - mean(x))^3) * sqrt(length(x)) / sum((x - mean(x))^2)^(3/2)
    }
    return(ans)
}

### Kurtosis
kurtosis <- function(x, excess = FALSE, normal = FALSE) {
    if (!is.null(ncol(x))) ans <- apply(x, 2, kurtosis, excess, normal)
    else if (normal) {
        x <- x[!is.na(x)]
        ans <- sum((x - mean(x))^4)/((length(x) - 1) * var(x)^2)
    } else {
        x <- x[!is.na(x)]
        ans <- sum((x - mean(x))^4) * length(x) / sum((x - mean(x))^2)^2
    }
    if (excess)	ans - 3
    return(ans)
}

### My lag function

lag.default <- function(x, k = -1, intersect = FALSE, ...) {
    if (isTRUE(all.equal(length(k), 1))) {
        ans <- stats:::lag.default(x, k)
    } else {
        ans <- NULL
        for (i in k) {
            if (isTRUE(all.equal(intersect, TRUE))) {
                ans <- ts.intersect(ans, stats:::lag.default(x, i))
            } else ans <- ts.union(ans, stats:::lag.default(x, i))
        }
        if (is.matrix(x)) {
            colnames(ans) <- paste(rep(colnames(x), each = length(k)), '.t', k, sep = '')
        } else colnames(ans) <- paste(deparse(substitute(x)),'.t',k, sep = '')
    }
    return(ans)
}

fracdiff <- function(x, k = 1) {
    if (all(abs(k %% 1) < 1e-12)) {
        if (abs(k) < 1e-12) {
            x
        } else diff(x, diff = k)
    } else {
        ## Note that initial values are not deleted for fracionally differenced series
        ## as any values is initial for some lagged value.
        if (length(k) != 1) stop('Not handled yet')
        ans <- .Fortran('fracdiff', as.numeric(x), as.integer(dim(x)), k)[[1]]
        dim(ans) <- dim(x)
        class(ans) <- class(x)
        tsp(ans) <- tsp(hasTsp(x))
        return(ans)
    }
}

### My LaTeX functions
LaTeX <- function(obj, ...) {
    if (missing(obj)) stop('No Object to tranform to LaTeX code')
    UseMethod('LaTeX')
}

LaTeX.default <- function(obj, ...) {
    
}

LaTeX.data.frame <- function(obj, na.sub = '', digits = 3, format = 'fg', ...) {
    ## Construct index of factor columns
    fac.id <- sapply(obj, is.factor) | sapply(obj, is.character)
    
    ## Construct index of unmric columns
    num.id <- sapply(obj, is.numeric)
	
    ## Coerce columns into character strings with the rigt formatting
    obj[fac.id] <- sapply(obj[fac.id], as.character)
    obj[num.id] <- sapply(obj[num.id], formatC, digits = digits, format = format)
    
    ## Substitute NAs
    obj[is.na(obj)] <- na.sub
    obj[num.id] <- sapply(obj[num.id], sub, pattern = 'NA', replacement = na.sub)
    
    ## Cheating solution untill time is in excess by using memisc package to print
    toLatex(obj, digit = digits)
}

LaTeX.matrix <- function(obj, na.sub = '', digits = 3, format = 'fg', se.symb1 = '(', se.symb2 = '[', ...) {
    NextMethod('LaTeX')
}

LaTeX.autotest.I1 <- function(obj, digits = 3, ...) {
    out <- coefmatrix(obj$value, obj$p.value, digits = digits)
    lout <- matrix(, 0, 3)
    colnames(lout) <- c('r', 'df', '5pct.c.v.')
    for (i in seq(length(obj$df))) {
        lout <- rbind(lout, c(i, obj$df[i], formatC(1 - pchisq(0.95, obj$df[i]), digits = digits, format = 'f')), "")
    }
    out <- cbind(lout, out)
    toLatex(out, digit = digits)
}

### Overwriting ts.union and ts.intersect to get better names
ts.union <- function(...) {
    ans <- stats:::.cbind.ts(list(...), "", dframe = FALSE, union = TRUE)
    colnames(ans) <- unlist(sapply(colnames(ans), function(x) sub('^.', '', x)))
    return(ans)
}

ts.intersect <- function(...) {
    ans <- stats:::.cbind.ts(list(...), "", dframe = FALSE, union = FALSE)
    colnames(ans) <- unlist(sapply(colnames(ans), function(x) sub('^.', '', x)))
    return(ans)
}
##
