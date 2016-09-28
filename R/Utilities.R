### FuzzyR - Utilities

## @export
meshgrid <- function(a,b) {
  list(x=outer(b*0, a, FUN="+"), y=outer(b, a*0, FUN="+"))
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

    if(is.matrix(x) && is.matrix(y)) {
        if(ncol(x) == 1 && ncol(y) != 1) {
            x <- c(x)
        } else if(ncol(y) == 1 && ncol(x) != 1) {
            y <- c(y)
        }
    }

    x*y
}

## This function can be replaced by the function expand.grid().
## e.g. expand.grid(c(list(1:3), list(1:2), list(1:3)))
## @export
genrule <- function(in_n, in_mf_n, pos=1) {

	if(pos < in_n) {
		sub_rule = genRule(in_n, in_mf_n, pos + 1)
	} else if(pos == in_n) {
		return(matrix(c(1:in_mf_n[pos]), ncol=1))
	} else {
		stop("error")
	}

	rule = NULL
	for(i in 1:in_mf_n[pos]) {
		tmp_rule = cbind(rep(i, nrow(sub_rule)), sub_rule)
		if(length(rule) != 0) {
			rule = rbind(rule, tmp_rule)
		} else {
			rule = tmp_rule
		}
	}
	rule
}



## Function: init.params.gbell
##  Description:
##      To initiate the parameters for gbell membership function
##  Input:
##      x: range of input, or a vector of numerical input
##      n: number of membership functions.
## @export
init.params.gbell <- function(x, n=2) {

	params <- NULL

	x.range <- range(x)

    a <- diff(x.range) / 5
    b <- 1
    c <- seq(x.range[1], x.range[2], length.out=n)

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
init.params.it2gbell <- function(x, n=2) {

	params <- NULL

	x.range <- range(x)

    a.lower <- diff(x.range) / 5
    a.upper = a.lower * 1.2
    b = 1
    c <- seq(x.range[1], x.range[2], length.out=n)

    params <- cbind(a.lower, a.upper, b, c)
	params
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

ekm <- function(wl, wr, f, maximum=F, w.which=F, sorted=F, k.which=F) {

    if(sum(wl > wr)) {
        stop("the lower bound should be no larger than upper bound")
    } else if(sum(wl < 0) || sum(wr <= 0)) {
        stop("the firing strength should not be negative number")
    }

    y <- NULL
    wr.which <- rep(TRUE, length(wr))

    #if(identical(wl, wr) && length(unique(wl)) == 1) {
    if(identical(wl, wr)) {
        #y <- mean(f)
        y <- sum(wl * f) / sum(wl)
    } else if(sum(wl) == 0) {
        if(maximum) {
            y <- max(f)
            wr.which[-which.max(f)] = FALSE
        } else {
            y <- min(f)
            wr.which[-which.min(f)] = FALSE
        }
    }

    if(is.null(y)) {

        idx.trim <- which(wr!=0)
        wl <- wl[idx.trim]
        wr <- wr[idx.trim]
        f <- f[idx.trim]

        if(length(unique(f)) == 1) {
            y <- unique(f)
        } else {
            n <- length(f)

            if(!sorted) {
                idx.order <- order(f)
                idx.trim <- idx.trim[idx.order]
                f <- f[idx.order]
                wl <- wl[idx.order]
                wr <- wr[idx.order]
            }

            if(maximum) {
                mflag <- 1
                k <- round(n / 1.7)
                w <- if(k < n) c(wl[1:k], wr[(k+1):n]) else wl[1:n]
            } else {
                mflag <- -1
                k <- round(n / 2.4)
                w <- if(k < n) c(wr[1:k], wl[(k+1):n]) else wr[1:n]
            }

            a <- w %*% f
            b <- sum(w)
            y <- as.numeric(a / b)
            ## The following 'complicated' operation is because the precision of R is terrible!
            ## In fact, y should be within [min(f), max(f)], but R does not give this fact!
            #k.tmp <- if(sum(round(f,15) >= round(y,15)) != 0) which.max(round(f,15) >= round(y,15)) else n
            k.tmp <- which.max(round(f,15) > round(y,15) - 2^-50)
            k.new <- if(k.tmp > 1) k.tmp - 1 else 1
            #k.new <- length(which(f<y+2^-50))
            #if(k.new == n) k.new <- n - 1
            k.hist <- k

            counter <- 1
            while(k.new != k && !(k.new %in% k.hist)) {
                counter <- counter + 1
                s <- sign(k.new - k)

                idx <- (min(k, k.new) + 1) : max(k, k.new)

                a <- a - mflag * s * (f[idx] %*% (wr[idx] - wl[idx]))
                b <- b - mflag * s * sum(wr[idx] - wl[idx])
                y <- as.numeric(a / b)

                k <- k.new
                k.hist <- c(k.hist, k)
                #k.tmp <- if(sum(round(f,15) >= round(y,15)) != 0) which.max(round(f,15) >= round(y,15)) else n
                k.tmp <- which.max(round(f,15) > round(y,15) - 2^-50)
                k.new <- if(k.tmp > 1) k.tmp - 1 else 1
                #k.new <- length(which(f<y+2^-50))
                #if(k.new == n) k.new <- n - 1
            }
            #cat("counter:", counter, "\n")
            #cat("k: ", k, "\n")

            if(maximum) {
                wr.which[idx.trim[1:k]] = FALSE
            } else {
                if(k < n) wr.which[idx.trim[(k+1):n]] = FALSE
            }
        }
    }

    wl.which <- !wr.which

    if(w.which) {
        cbind(wl.which, wr.which)
    } else {
        if(k.which) {
            c(y, k) 
        } else {
            y
        }
    }
}
