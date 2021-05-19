### FuzzyR - Fuzzy (Set) Operation


#' @title Fuzzification
#' @description
#' To convert the crisp input x to a fuzzy membership function with specified fuzzification method
#' @param fuzzification.method The fuzzification method
#' @param x The required parameters for a fuzzification method
#' @param mf.params The parameters for a membership function
#' @return The corresponding fuzzy membership function
#' @examples
#' x <- 3
#' mf <- x.fuzzification(gbell.fuzzification, x, c(1,2))
#' # This is the same as:
#' mf <- genmf(gbellmf, c(1,2,x))
#'
#' evalmf(1:10, mf)
#' @author Chao Chen
#' @export

x.fuzzification <- function(fuzzification.method, x, mf.params) {
    F <- match.fun(fuzzification.method)
    F(x, mf.params)
}


#' @title Fuzzy tnorm
#' @description
#' To conduct t-norm operation for given fuzzy member functions
#' @param operator The t-norm operator such as min, prod
#' @param ... fuzzy membership functions
#' @return A membership function, which is the t-norm of membership functions
#' @examples
#' mf1 <- genmf(gbellmf, c(1,2,3))
#' mf2 <- genmf(gbellmf, c(4,5,6))
#' mf3 <- fuzzy.tnorm(prod, mf1, mf2)
#' tmp1 <- evalmf(1:10, mf1)
#' tmp2 <- evalmf(1:10, mf2)
#' tmp3 <- evalmf(1:10, mf3)
#' identical(tmp3, tmp1*tmp2)
#' tmp3
#' @author Chao Chen
#' @export

fuzzy.tnorm <- function(operator, ...) {

    operator <- match.fun(operator)
    operator.list <- c(min, prod)

    if (!fun.exists(operator, operator.list)) {
        stop("unsupported t-norm operator")
    }

    fuzzy.t(operator, ...)
}


#' @title Fuzzy t-conorm
#' @description
#' To conduct t-conorm operation for given fuzzy member functions
#' @param operator The t-conorm operator such as max
#' @param ... fuzzy membership functions
#' @return A membership function, which is the t-conorm of membership functions
#' @examples
#' mf1 <- genmf(gbellmf, c(1,2,3))
#' mf2 <- genmf(gbellmf, c(4,5,6))
#' mf3 <- fuzzy.tconorm(max, mf1, mf2)
#' tmp1 <- evalmf(1:10, mf1)
#' tmp2 <- evalmf(1:10, mf2)
#' tmp3 <- evalmf(1:10, mf3)
#' identical(tmp3, pmax(tmp1, tmp2))
#' tmp3
#' @author Chao Chen
#' @export

fuzzy.tconorm <- function(operator, ...) {

    operator <- match.fun(operator)
    operator.list <- c(max)

    if (!fun.exists(operator, operator.list)) {
        stop("unsupported t-conorm operator")
    }

    fuzzy.t(operator, ...)
}


#' @title Fuzzy t-norm/t-conorm operation
#' @description
#' To conduct t-norm or t-conorm operation for given fuzzy member functions
#' @param operator The supported t-norm/t-conorm operators are min, prod, max
#' @param ... fuzzy membership functions
#' @return A membership function, which is the t-norm/t-conorm of membership functions
#' @examples
#' mf1 <- genmf(gbellmf, c(1,2,3))
#' mf2 <- genmf(gbellmf, c(4,5,6))
#' mf3 <- fuzzy.t(max, mf1, mf2)
#' tmp1 <- evalmf(1:10, mf1)
#' tmp2 <- evalmf(1:10, mf2)
#' tmp3 <- evalmf(1:10, mf3)
#' identical(tmp3, pmax(tmp1, tmp2))
#' tmp3
#' @author Chao Chen
#' @export

fuzzy.t <- function(operator, ...) {

    operator <- match.fun(operator)
    mf.list <- list(...)

    fuzzy.t <- function(x) {
        mfx <- sapply(mf.list, function(mf) { FUN <- match.fun(mf); FUN(x) })
        mfx <- matrix(mfx, nrow = length(x))
        apply(mfx, 1, operator)
    }
}


#' @title Fuzzy rule firing
#' @description
#' To get the firing strength for the given input fuzzification membership function and the antecedent membership function in the domain of [lower, upper]
#' @param operator t-norm operator
#' @param x.mf the fuzzy input membership function
#' @param ante.mf the antecedent membership function
#' @param lower lower bound of the input
#' @param upper upper bound of the input
#' @return the rule firing strenth
#' @examples
#' x.mf <- x.fuzzification(gbell.fuzzification, 3, c(1,2))
#' ante.mf <- genmf(gbellmf, c(1,2,6))
#' firing.strength <- fuzzy.firing(min, x.mf, ante.mf, lower=0, upper=10)
#' firing.strength
#' @author Chao Chen
#' @export

fuzzy.firing <- function(operator, x.mf, ante.mf, lower, upper) {

    x.mf <- c(x.mf)
    ante.mf <- c(ante.mf)

    firing.mf <- fuzzy.tnorm(operator, x.mf[[1]], ante.mf[[1]])
    firing.strength <- fuzzy.optimise(firing.mf, lower, upper)

    if (length(x.mf) > 1 || length(ante.mf) > 1) {
        firing.mf <- fuzzy.tnorm(operator, x.mf[[length(x.mf)]], ante.mf[[length(ante.mf)]])
        #firing.mf <- fuzzy.tnorm(operator, tail(x.mf, 1)[[1]], tail(ante.mf,1)[[1]])
        firing.strength <- c(firing.strength, fuzzy.optimise(firing.mf, lower, upper))
    }

    firing.strength
}


fuzzy.firing.defuzz <- function(operator, x.mf, ante.mf, lower, upper, defuzz.type) {

    x.mf <- c(x.mf)
    ante.mf <- c(ante.mf)

    firing.mf <- fuzzy.tnorm(operator, x.mf[[1]], ante.mf[[1]])
    x <- seq(lower, upper, len = 100)
    mu <- evalmf(x, firing.mf)
    x.defuzz <- defuzz(x, mu, defuzz.type)
    firing.strength <- evalmf(x.defuzz, firing.mf)

    if (length(x.mf) > 1 || length(ante.mf) > 1) {
        firing.mf <- fuzzy.tnorm(operator, x.mf[[length(x.mf)]], ante.mf[[length(ante.mf)]])
        mu <- evalmf(x, firing.mf)
        x.defuzz <- defuzz(x, mu, defuzz.type)
        firing.strength <- c(firing.strength, evalmf(x.defuzz, firing.mf))
    }

    firing.strength
}


fuzzy.firing.similarity.set <- function(x.mf, ante.mf, lower, upper) {
    x.mf <- c(x.mf)
    ante.mf <- c(ante.mf)

    mf.intersection <- fuzzy.tnorm('min', x.mf[[1]], ante.mf[[1]])
    mf.union <- fuzzy.tconorm('max', x.mf[[1]], ante.mf[[1]])

    integration.intersection <- integrate(mf.intersection, lower, upper)$value
    integration.union <- integrate(mf.union, lower, upper)$value
    firing.strength <- integration.intersection / integration.union

    if (length(x.mf) > 1 || length(ante.mf) > 1) {
        mf.intersection <- fuzzy.tnorm('min', x.mf[[length(x.mf)]], ante.mf[[length(ante.mf)]])
        mf.union <- fuzzy.tconorm('max', x.mf[[length(x.mf)]], ante.mf[[length(ante.mf)]])

        integration.intersection <- integrate(mf.intersection, lower, upper)$value
        integration.union <- integrate(mf.union, lower, upper)$value
        firing.strength <- c(firing.strength, integration.intersection / integration.union)
    }

    firing.strength
}


#' @title Fuzzy optimisation
#' @description
#' to get an approximation of the maximum membership grade for a given membership function in the domain of [lower, upper]
#' @param fuzzy.mf fuzzy member function
#' @param lower lower bound of the input
#' @param upper upper bound of the input
#' @return an approximation of the maximum membership grade in the given domain
#' @examples
#' mf <- genmf(gbellmf, c(1,2,3))
#' x <- seq(4, 5, by=0.01)
#' max(evalmf(x, mf))
#' fuzzy.optimise(mf, 4, 5)
#' @author Chao Chen
#' @export

fuzzy.optimise <- function(fuzzy.mf, lower, upper) {
    if (lower == upper) {
        # This is for the firing of singleton fuzzification, where lower==upper should be x' to get f(x')
        fuzzy.mf(lower)
    } else {
        max1 <- optimise(fuzzy.mf, c(lower, upper), maximum = T, tol = 6.32e-10)$objective
        max2 <- fuzzy.mf(seq(lower, upper, length.out=100))
        max(max1, max2)
    }
}


## Function: ekm
##  Description:
##          to optimize with the EKM algorithm
##      OR
##          to obtain which wl and wr for optimizing the weigthed mean with the EKM algorithm
##  Input:
##      wl: the lower bound vector of the firing intervals (wl >= 0)
##      wr: the upper bound vector of the firing intervals (wr >= 0)
##      f: the lower/upper bound vector of the tsk outputs
##      maximum: T for maximum, and F for minimum (default)
##      w.which: if TRUE, then return the list of wl and wr which are used to achieve the optimum
##  Output:
##          y: the optimized output
##      OR
##          wl.which: a list of which wl to be used
##          wr.which: a list of which wr to be used
## @author Chao Chen

ekm <- function(wl, wr, f, maximum = F, w.which = F, sorted = F, k.which = F) {

    if (sum(wl > wr)) {
        stop("the lower bound should be no larger than upper bound")
    } else if (sum(wl < 0) || sum(wr < 0)) {
        stop("the firing strength should not be negative number")
    }

    y <- NULL
    wr.which <- rep(TRUE, length(wr))
    k <- NULL

    s1 = sum(wl)
    if (identical(wl, wr)) {
        if (s1 != 0) {
            y <- wl %*% f / s1
        } else {
            y <- mean(f)
        }
    } else if (s1 == 0) {
        if (maximum) {
            y <- max(f[wr != 0])
        } else {
            y <- min(f[wr != 0])
        }
        wr.which[-which(f == y)] = FALSE
    }

    if (is.null(y)) {

        idx.trim <- which(wr != 0)
        wl <- wl[idx.trim]
        wr <- wr[idx.trim]
        f <- f[idx.trim]

        if (length(unique(f)) == 1) {
            y <- unique(f)
        } else {
            n <- length(f)

            if (!sorted) {
                idx.order <- order(f)
                idx.trim <- idx.trim[idx.order]
                f <- f[idx.order]
                wl <- wl[idx.order]
                wr <- wr[idx.order]
            }

            if (maximum) {
                mflag <- 1
                k <- round(n / 1.7)
                w <- if (k < n) c(wl[1:k], wr[(k + 1):n]) else wl[1:n]
            } else {
                mflag <- -1
                k <- round(n / 2.4)
                w <- if (k < n) c(wr[1:k], wl[(k + 1):n]) else wr[1:n]
            }

            a <- w %*% f
            b <- sum(w)
            y <- as.numeric(a / b)
            ## The following 'complicated' operation is because the precision limit of R or CPU!
            ## This may lead to wrong results! We suggest an alternative function km.da.
            ## In fact, y should be within [min(f), max(f)], but R does not give this fact in some case!
            #k.tmp <- if(sum(round(f,15) >= round(y,15)) != 0) which.max(round(f,15) >= round(y,15)) else n
            #            k.tmp <- which.max(round(f,15) > round(y,15) - 2^-50)
            #            k.new <- if(k.tmp > 1) k.tmp - 1 else 1
            k.new <- length(which(f <= y))
            if (k.new == n) k.new <- n - 1
            k.hist <- k

            #            cnt <- 1
            while (k.new != k && !(k.new %in% k.hist)) {
                #                cnt <- cnt + 1
                s <- sign(k.new - k)

                idx <- (min(k, k.new) + 1):max(k, k.new)

                a <- a - mflag * s * (f[idx] %*% (wr[idx] - wl[idx]))
                b <- b - mflag * s * sum(wr[idx] - wl[idx])
                y <- as.numeric(a / b)

                k <- k.new
                k.hist <- c(k.hist, k)
                #k.tmp <- if(sum(round(f,15) >= round(y,15)) != 0) which.max(round(f,15) >= round(y,15)) else n
                #                k.tmp <- which.max(round(f,15) > round(y,15) - 2^-50)
                #                k.new <- if(k.tmp > 1) k.tmp - 1 else 1
                k.new <- length(which(f < y + 2^-50))
                if (k.new == n) k.new <- n - 1
            }
            #            cat("cnt:", cnt, "\n")
            #            cat("k.hist: [", k.hist, "]\n")
            #cat("k: ", k, "\n")

            if (maximum) {
                wr.which[idx.trim[1:k]] = FALSE
            } else {
                if (k < n) wr.which[idx.trim[(k + 1):n]] = FALSE
            }
        }
    }

    wl.which <- !wr.which

    if (w.which) {
        cbind(wl.which, wr.which)
    } else {
        if (k.which) {
            c(y, k)
        } else {
            y
        }
    }
}


#' @title km.da
#' @description
#' A Direct Approach for Determining the Switch Points in the Karnik-Mendel Algorithm.
#' @param wl A vector of lower membership grades.
#' @param wr A vector of upper membership grades.
#' @param f A vector of the primary values in the discrete universe of discourse X.
#' @param maximum T, to calculate the maximum centroid; F, to calulate the minimum centroid.
#' @param w.which T, to show which membership grade to be used to calculate maximum/minimum centroid for each primary value.
#' @param sorted T, to indicate that the primary values have already been put in ascending order.
#' @param k.which T, to show the index of the switch point selected by the algorithm.
#' @return w.which=T, a two-column matrix indicating which membership grades to be used;
#' w.which=F and k.which=T, a vector of the centroid and the switch point;
#' w.which=F and k.which=F, a single value of the centroid.
#' @examples
#' wr <- runif(100, 0, 1)
#' wl <- wr * runif(100, 0, 1)
#' f <- abs(runif(100, 0, 1))
#' f <- sort(f)
#' km.da(wl, wr, f)
#' @author Chao Chen
#' @references
#' [1] C. Chen, R. John, J. Twycross, and J. M. Garibaldi, “A Direct Approach for Determining the Switch Points in the Karnik–Mendel Algorithm,” IEEE Transactions on Fuzzy Systems, vol. 26, no. 2, pp. 1079–1085, Apr. 2018. \cr
#' \doi{10.1109/TFUZZ.2017.2699168}
#'
#' [2] C. Chen, D. Wu, J. M. Garibaldi, R. John, J. Twycross, and J. M. Mendel, “A Comment on ‘A Direct Approach for Determining the Switch Points in the Karnik-Mendel Algorithm,’” IEEE Transactions on Fuzzy Systems, vol. 26, no. 6, pp. 3905–3907, 2018. \cr
#' \doi{10.1109/TFUZZ.2018.2865134}
#' @export

km.da <- function(wl, wr, f, maximum = F, w.which = F, sorted = F, k.which = F) {

    if (sum(wl > wr)) {
        stop("the lower bound should be no larger than upper bound")
    } else if (sum(wl < 0) || sum(wr < 0)) {
        stop("the firing strength should not be negative number")
    }

    y <- NULL
    wr.which <- rep(TRUE, length(wr))
    k <- NULL

    s1 = sum(wl)
    if (identical(wl, wr)) {
        if (s1 != 0) {
            y <- wl %*% f / s1
        } else {
            y <- mean(f)
        }
    } else if (s1 == 0) {
        if (maximum) {
            y <- max(f[wr != 0])
        } else {
            y <- min(f[wr != 0])
        }
        wr.which[-which(f == y)] = FALSE
    }

    if (is.null(y)) {

        idx.trim <- which(wr != 0)
        wl <- wl[idx.trim]
        wr <- wr[idx.trim]
        f <- f[idx.trim]

        if (length(unique(f)) == 1) {
            y <- unique(f)
        } else {
            n <- length(f)

            if (!sorted) {
                idx.order <- order(f)
                idx.trim <- idx.trim[idx.order]
                f <- f[idx.order]
                wl <- wl[idx.order]
                wr <- wr[idx.order]
            }

            if (maximum) {
                positive.w.cusum <- cumsum(wl[1:(n - 1)])
                negative.w.cusum <- cumsum(wr[n:2])
            } else {
                positive.w.cusum <- cumsum(wr[1:(n - 1)])
                negative.w.cusum <- cumsum(wl[n:2])
            }

            delta.f <- diff(f)
            positive <- delta.f * positive.w.cusum
            negative <- delta.f[(n - 1):1] * negative.w.cusum

            positive.cusum <- c(0, cumsum(positive))
            negative.cusum <- c(cumsum(negative)[(n - 1):1], 0)

            derivatives <- positive.cusum - negative.cusum
            #k <- which.max(derivatives>=0)-1
            #k <- sum(derivatives<0)
            k <- length(which(derivatives < 0))
            if (k == 0) k = 1

            if (maximum) {
                w <- c(wl[1:k], wr[(k + 1):n])
            } else {
                w <- c(wr[1:k], wl[(k + 1):n])
            }

            a <- w %*% f
            b <- sum(w)
            y <- as.numeric(a / b)

            #cat("k: ", k, "\n")

            if (maximum) {
                wr.which[idx.trim[1:k]] = FALSE
            } else {
                if (k < n) wr.which[idx.trim[(k + 1):n]] = FALSE
            }
        }
    }

    wl.which <- !wr.which

    if (w.which) {
        cbind(wl.which, wr.which)
    } else {
        if (k.which) {
            c(y, k)
        } else {
            c(y)
        }
    }
}

