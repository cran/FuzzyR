### FuzzyR - Utilities


## @export
meshgrid <- function(a, b) {
    list(x = outer(b * 0, a, FUN = "+"), y = outer(b, a * 0, FUN = "+"))
}

## Function: fun.exists
##  Description:
##      to check whether a function exists in a given function list
##  Input:
##      fun: function to be checked
##      where: function list
##  Output:
##      TRUE/FALSE
## @export
fun.exists <- function(fun, where) {
    sum(sapply(unlist(where), function(x) identical(match.fun(x), match.fun(fun))))
}


#' @title fuzzyr.match.fun
#' @description
#' This is a modification of the original match.fun, where parent.frame(2) is changed to parent.env(environment()).
#' @param FUN item to match as function: a function, symbol or character string.
#' @param descend logical; control whether to search past non-function objects.
#' @details See \code{\link{match.fun}}.
#' @export
fuzzyr.match.fun <- function(FUN, descend = TRUE)
{
    if (is.function(FUN))
        return(FUN)
    if (!(is.character(FUN) && length(FUN) == 1L || is.symbol(FUN))) {
        FUN <- eval.parent(substitute(substitute(FUN)))
        if (!is.symbol(FUN))
            stop(gettextf("'%s' is not a function, character or symbol",
                          deparse(FUN)), domain = NA)
    }
    envir <- parent.env(environment())
    if (descend)
        FUN <- get(as.character(FUN), mode = "function", envir = envir)
    else {
        FUN <- get(as.character(FUN), mode = "any", envir = envir)
        if (!is.function(FUN))
            stop(gettextf("found non-function '%s'", FUN), domain = NA)
    }
    return(FUN)
}

match.fun <- fuzzyr.match.fun


## Function: in.range
##  Description:
##      to check whether the value x is in a given range
##  Input:
##      x: numeric or vector
##      range: range(a,b)
##  Output:
##      a vector of TRUE/FALSE
## @export
in.range <- function(x, range) {
    ifelse(x >= range[1] & x <= range[2], TRUE, FALSE)
}

## Function: multiply
##  Description:
##      An alternative function of `*` which allows in-compatible matrix multiplication
##  Input:
##      x,y: numeric or complex vectors or objects which can be coerced to such, or other objects for which methods have been written.
##  Output:
##      the same as the function `*`
## @export
multiply <- function(x, y) {

    if (is.matrix(x) && is.matrix(y)) {
        if (ncol(x) == 1 && ncol(y) != 1) {
            x <- c(x)
        } else if (ncol(y) == 1 && ncol(x) != 1) {
            y <- c(y)
        }
    }

    x * y
}


## Function: init.params.gbell
##  Description:
##      To initiate the parameters for gbell membership function
##  Input:
##      x: range of input, or a vector of numerical input
##      n: number of membership functions.
## @export
init.params.gbell <- function(x, n = 2) {

    params <- NULL

    x.range <- range(x)

    a <- diff(x.range) / ((n - 1) * 2)
    if (a == 0) a = 1
    b <- 1
    c <- seq(x.range[1], x.range[2], length.out = n)

    params <- cbind(a, b, c)
    params
}


## Function: init.params.it2gbell
##  Description:
##      To initiate the parameters for it2gbell membership function
##  Input:
##      x: range of input, or a vector of numerical input
##      n: number of membership functions.
## @export
init.params.it2gbell <- function(x, n = 2) {

    params <- NULL

    x.range <- range(x)

    a.lower <- diff(x.range) / 5
    if (a.lower == 0) a.lower = 1
    a.upper = a.lower * 1.2
    b = 1
    c <- seq(x.range[1], x.range[2], length.out = n)

    params <- cbind(a.lower, a.upper, b, c)
    params
}


## This function can be replaced by the function expand.grid().
## e.g. expand.grid(c(list(1:3), list(1:2), list(1:3)))
## @export
genrule <- function(in_n, in_mf_n, pos = 1) {

    if (pos < in_n) {
        sub_rule = genRule(in_n, in_mf_n, pos + 1)
    } else if (pos == in_n) {
        return(matrix(c(1:in_mf_n[pos]), ncol = 1))
    } else {
        stop("error")
    }

    rule = NULL
    for (i in 1:in_mf_n[pos]) {
        tmp_rule = cbind(rep(i, nrow(sub_rule)), sub_rule)
        if (length(rule) != 0) {
            rule = rbind(rule, tmp_rule)
        } else {
            rule = tmp_rule
        }
    }
    rule
}


rspe <- function(f, y, ref = 0) {

    sign <- ifelse(y >= 0,
                   ifelse((f - y) >= 0, 1, -1),
                   ifelse((f - y) <= 0, 1, -1)
    )

    rspe <- sign * abs(f - y) / (abs(f - y) + abs(y - ref))
    rspe[f == y] = 0
    rspe[f == ref] = sign[f == ref] * 0.5

    rspe
}

#' @title Fuzzy Accuracy
#' @description
#' This function is to provide performance indicators by using eight different accuracy measures including a new measure UMBRAE.
#' @param f A vector of forecasting values produced by a model to be evaluated.
#' @param y A vector of observed values.
#' @param f.ref A vector of forecasting values produced by a benchmark method to be compared.
#' @param scale.mase A single value which is the scaling factor of the measure MASE.
#' @return A vector of results by each measure.
#' @examples
#' f <- rnorm(10)
#' y <- rnorm(10)
#' fuzzyr.accuracy(f, y)
#' @author Chao Chen
#' @references
#' [1] C. Chen, J. Twycross, and J. M. Garibaldi, “A new accuracy measure based on bounded relative error for time series forecasting,” PLOS ONE, vol. 12, no. 3, pp. 1–23, 2017. \cr
#' \doi{10.1371/journal.pone.0174202}
#' @export

fuzzyr.accuracy <- function(f, y, f.ref = 0, scale.mase = NULL) {

    e <- y - f

    mae <- mean(abs(e), na.rm = TRUE)

    mse <- mean(e^2, na.rm = TRUE)
    rmse <- sqrt(mse)

    mase <- NA
    if (!is.null(scale.mase)) {
        mase <- mean(abs(e / scale.mase), na.rm = TRUE)
    }

    e.ref <- y - f.ref
    mrae <- mean(abs(e / e.ref), na.rm = TRUE)
    gmrae <- exp(mean(log(abs(e / e.ref)), na.rm = TRUE))

    pe <- e / y * 100
    mape <- mean(abs(pe), na.rm = TRUE)

    spe <- e / (abs(y) + abs(f)) * 200
    smape <- mean(abs(spe), na.rm = TRUE)

    mbrae <- mean(abs(rspe(f, y, f.ref)))
    umbrae <- mbrae / (1 - mbrae)

    out <- c(mae, rmse, mase, mrae, gmrae, mape, smape, umbrae)
    names(out) <- c("MAE", "RMSE", "MASE", "MRAE", "GMRAE", "MAPE", "sMAPE", "uMbRAE")

    out
}


#' Plot firing strength with different inference method
#'
#' Plots a 2D graph of the firing strength for a antecedent produced by different inference method
#'
#' @param IP A matrix representing the input stack, number of inputs (columns) by number of outputs (rows).
#' @param mfType The type of fuzzy membership function\cr
#' @param mfPara The parameters for the given type of membership function\cr
#' @param fuzMethod The type of fuzzy membership function for non-singleton fuzzification
#' @param fuzPara The parameters for the given fuz.type of membership function
#' @param SFLS When TRUE, shows the firing strength produced by SFLS
#' @param STD When TRUE, shows the firing strength produced by std-NSFLS
#' @param CEN When TRUE, shows the firing strength produced by cen-NSFLS
#' @param SIM When TRUE, shows the firing strength produced by sim-NSFLS
#' @param step For discrete fuzzification
#' @param fisRange Field of definition, for example, c(1,10)
#' @return A two dimensional graph displaying all the firing strength produced by given method.
#' @examples
#' cmp.firing(1, 'gaussmf', c(1, 2.5, 1), 'gbell', c(0.4, 2), step=100)
#'
#' @author Yu Zhao
#' @export
cmp.firing <- function(IP, mfType, mfPara, fuzMethod, fuzPara, SFLS = TRUE, STD = TRUE, CEN = FALSE, SIM = FALSE, step = 100, fisRange = NULL) {
    x <- IP
    if (!is.null(fisRange)) {
        Rmin <- fisRange[1]
        Rmax <- fisRange[2]
    }
    for (i in 1:length(x)) {
        if (typeof(fuzMethod) == "character") {
            if (fuzMethod == "gauss") {
                if (!is.null(fisRange)) {
                    xx <- seq(Rmin, Rmax, (Rmax - Rmin) / step)
                }
                else {
                    Range <- fuzPara * 4
                    xx <- seq(x[i] - Range, x[i] + Range, 2 * Range / step)
                    fisRange <- c(x[i] - Range, x[i] + Range)
                }
                MF <- gauss.fuzzification(x[i], fuzPara)
                y1 <- sapply(c(MF), function(F) F(xx))
                Input <- xx
                Degree <- y1
                plot(Input, Degree, type = "l", lwd = 1)
            } else if (fuzMethod == "gbell") {
                if (!is.null(fisRange)) {
                    xx <- seq(Rmin, Rmax, (Rmax - Rmin) / step)
                }
                else {
                    Range <- fuzPara[1] * 2 + (fuzPara[1] * 2 * (2 / fuzPara[2]))
                    xx <- seq(x[i] - Range, x[i] + Range, 2 * Range / step)
                    fisRange <- c(x[i] - Range, x[i] + Range)
                }
                MF <- gbell.fuzzification(x[i], fuzPara)
                y1 <- sapply(c(MF), function(F) F(xx))
                Input <- xx
                Degree <- y1
                plot(Input, Degree, type = "l", lwd = 1)
            } else if (fuzMethod == "tri") {
                para <- fuzPara
                if (length(para) == 1) {
                    para <- c(para, para)
                }
                if (!is.null(fisRange)) {
                    xx <- seq(Rmin, Rmax, (Rmax - Rmin) / step)
                }
                else {
                    Range <- sum(para)
                    xx <- seq(x[i] - (para[1] * 2), x[i] + (para[2] * 2), Range / step)
                    fisRange <- c(x[i] - (para[1] * 2), x[i] + (para[2] * 2))
                }
                MF <- tri.fuzzification(x[i], fuzPara)
                y1 <- sapply(c(MF), function(F) F(xx))
                Input <- xx
                Degree <- y1
                plot(Input, Degree, type = "l", lwd = 1)
            }
        }
        else {
            if (!is.null(fisRange)) {
                xx <- seq(Rmin, Rmax, (Rmax - Rmin) / step)
            }
            else {
                stop("Field of definition (fisRange) must be given.")
            }
            MF <- fuzMethod(x[i], fuzPara)
            y1 <- sapply(c(MF), function(F) F(xx))
            Input <- xx
            Degree <- y1
            plot(Input, Degree, type = "l", lwd = 1)
        }

        MF <- genmf(mf.type = mfType, mf.params = mfPara)
        y2 <- sapply(c(MF), function(F) F(xx))
        points(xx, y2, type = "l", lwd = 1)
        y <- apply(cbind(y1, y2), 1, min)
        points(xx, y, type = "l", col = "tomato", lwd = 3)
        if (SFLS) {
            abline(v = x[i], lty = 3, col = "forestgreen", lwd = 2)
            yy <- evalmf(IP, mfType, mfPara)
            legend("topleft", pch = 16, legend = paste(round(yy, 3), "(SFLS)"),
                   col = "forestgreen", cex = 0.75, bty = "n")
        }
        std <- mean(xx[y == max(y)])
        if (STD) {
            abline(v = std, lty = 3, col = "orange", lwd = 2)
            zz <- evalmf(IP, mfType, mfPara, fuzMethod, fuzPara, "tnorm.min.max", fisRange)
            legend("topleft", inset = c(0, 0.1), pch = 16, legend = paste(round(zz, 3), "(Std-NSFLS)"),
                   col = "orange", cex = 0.75, bty = "n")
        }
        if (!max(y) == 0) {
            y <- y / sum(y)
            x[i] <- sum(xx * y)
        }
    }
    if (CEN) {
        abline(v = x[i], lty = 3, col = "blue4", lwd = 2)
        xx <- evalmf(IP, mfType, mfPara, fuzMethod, fuzPara, "tnorm.min.defuzz.centroid", fisRange)
        legend("topright", pch = 16, legend = paste(round(xx, 3), "(Cen-NSFLS)"),
               col = "blue4", cex = 0.75, bty = "n")
    }
    if (SIM) {
        aa <- evalmf(IP, mfType, mfPara, fuzMethod, fuzPara, "similarity.set", fisRange)
        abline(h = aa, lty = 3, col = "brown", lwd = 2)
        legend("topright", inset = c(0, 0.1), pch = 16, legend = paste(round(aa, 3), "(Sim-NSFLS)"),
               col = "brown", cex = 0.75, bty = "n")
    }
    title(main = paste("Firing strength ( input =", IP, ")"))
}


