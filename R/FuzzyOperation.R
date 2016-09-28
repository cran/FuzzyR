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

    if(!fun.exists(operator, operator.list)) {
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

    if(!fun.exists(operator, operator.list)) {
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

    fuzzy.t <- function(x) {
        mf.list <- list(...)

        mfx <- sapply(mf.list, function(mf) {FUN <- match.fun(mf); FUN(x)})
        mfx <- matrix(mfx, nrow=length(x))
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

    if(length(x.mf) > 1 || length(ante.mf) > 1) {
        firing.mf <- fuzzy.tnorm(operator, x.mf[[length(x.mf)]], ante.mf[[length(ante.mf)]])
        #firing.mf <- fuzzy.tnorm(operator, tail(x.mf, 1)[[1]], tail(ante.mf,1)[[1]])
        firing.strength <- c(firing.strength, fuzzy.optimise(firing.mf, lower, upper))
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

    if(lower == upper) {
        # This is for the firing of singleton fuzzification, where lower==upper should be x' to get f(x')
        fuzzy.mf(lower)
    } else {
        optimise(fuzzy.mf, c(lower, upper), maximum=T, tol=10^-6)$objective
    }
}
