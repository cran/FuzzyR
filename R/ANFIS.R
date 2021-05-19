### FuzzyR - ANFIS Toolkit

# require(plyr)

#' @import plyr
NULL

#' @title ANFIS model builder
#' @description
#' To build an ANFIS model from an existing FIS model
#' @param fis A fuzzy inference system model initialised by \code{\link{newfis}}.
#' @return An ANFIS model
#' @examples
#' fis <- anfis.tipper()
#' anfis <- anfis.builder(fis)
#' @author Chao Chen
#' @references
#' [1] C. Chen, R. John, J. Twycross, and J. M. Garibaldi, “An extended ANFIS architecture and its learning properties for type-1 and interval type-2 models,” in Proceedings IEEE International Conference on Fuzzy Systems, 2016, pp. 602–609. \cr
#' \doi{10.1109/FUZZ-IEEE.2016.7737742}
#'
#' [2] C. Chen, R. John, J. Twycross, and J. M. Garibaldi, “Type-1 and interval type-2 ANFIS: a comparison,” in Proceedings IEEE International Conference on Fuzzy Systems, 2017, pp. 1–6.  \cr
#' \doi{10.1109/FUZZ-IEEE.2017.8015555}
#' @export

anfis.builder <- function(fis) {

    LI <- 1
    L1 <- 2
    L2 <- 3
    L3 <- 4
    L4 <- 5
    L5 <- 6

    input.num <- length(fis$input)
    input.mf.num <- NULL
    for (i in 1:input.num)
        input.mf.num[i] <- length(fis$input[[i]]$mf)

    output.num <- length(fis$output)
    output.mf.num <- NULL
    for (i in 1:output.num)
        output.mf.num[i] <- length(fis$output[[i]]$mf)

    if (output.num != 1) {
        stop("the fis with multiple outputs is not supported")
    }

    rule.num <- nrow(fis$rule)
    rule.antecedent <- matrix(fis$rule[, 1:input.num], nrow = rule.num)
    rule.consequent <- matrix(fis$rule[, (input.num + 1):(input.num + output.num)], nrow = rule.num)
    rule.weight <- matrix(fis$rule[, ncol(fis$rule) - 1], nrow = rule.num)
    rule.operator <- matrix(fis$rule[, ncol(fis$rule)], nrow = rule.num)

    ## initiate the structure of anfis
    # TODO the main properties of fis such as type, methods are also considered to be added
    anfis <- list(andMethod = fis$andMethod, orMethod = fis$orMethod, impMethod = fis$impMethod, aggMethod = fis$aggMethod,
                  defuzzMethod = fis$defuzzMethod, rule = fis$rule, layer = NULL)

    ## to build layer I, the input layer
    anfis <- anfis.addlayer(anfis, "inputs")

    offset <- 0
    for (i in 1:input.num) {
        fuzzification.method <- fis$input[[i]]$fuzzification.method
        if (grepl(".fuzzification", fuzzification.method, fixed = T) == FALSE) {
            fuzzification.method <- paste0(fuzzification.method, ".fuzzification")
        }
        anfis <- anfis.addnode(anfis, LI, NULL, NULL, NULL, c((1 + offset):(input.mf.num[i] + offset)),
                               fis$input[[i]]$fuzzification.method, fis$input[[i]]$fuzzification.params,
                               fis$input[[i]]$range, fis$input[[i]]$name)
        offset = offset + input.mf.num[i]
    }

    ## build layer 1 with nodes of rule antecedents (input member functions)
    anfis <- anfis.addlayer(anfis, "antecedent")

    for (i in 1:input.num) {
        for (j in 1:input.mf.num[i]) {
            fan.in <- c(i)
            fan.out <- which(rule.antecedent[, i] == j)
            anfis <- anfis.addnode(anfis, L1, NULL, NULL, fan.in,
                                   fan.out, fis$input[[i]]$mf[[j]]$type, fis$input[[i]]$mf[[j]]$params,
                                   fis$input[[i]]$range, fis$input[[i]]$mf[[j]]$name)
        }
    }


    #build layer 2 which contains nodes of rule firing strength
    #build layer 3 which contains nodes of normalised firing strength
    #build layer 4 which contains nodes of rule consequents
    anfis <- anfis.addlayer(anfis, "firing strength")
    anfis <- anfis.addlayer(anfis, "normalised firing strength")
    anfis <- anfis.addlayer(anfis, "consequent")

    for (i in 1:rule.num) {
        ## add nodes for layer 2
        offset <- cumsum(c(0, input.mf.num))[1:input.num]
        fan.in <- rule.antecedent[i,] + offset
        fan.in <- fan.in[which(rule.antecedent[i,] != 0)]
        anfis <- anfis.addnode(anfis, L2, NULL, NULL, fan.in, c(1:rule.num), NULL, NULL)

        ## add nodes for layer 3
        anfis <- anfis.addnode(anfis, L3, NULL, NULL, c(1:rule.num), c(i), NULL, NULL)

        ## TODO: multiple-output ANFIS is under consideration.
        ##       Also, currently, the number of nodes in layer 4 is the same as the number of rules.
        ##       Less number of nodes need to be considered.
        ## add nodes for layer 4
        fan.in <- c(i)
        fan.out <- which(rule.consequent[i,] != 0)
        #offset <- cumsum(c(0, output.mf.num))[1:output.num]
        idx <- rule.consequent[i, 1]
        mf.type <- fis$output[[1]]$mf[[idx]]$type
        mf.params <- fis$output[[1]]$mf[[idx]]$params
        anfis <- anfis.addnode(anfis, L4, NULL, NULL, fan.in, fan.out, mf.type, mf.params)
    }


    #build layer 5 which contains nodes of final output with defuzziation
    anfis <- anfis.addlayer(anfis, "output")
    anfis <- anfis.addnode(anfis, L5, NULL, NULL, c(1:rule.num), NULL, NULL, NULL)
    anfis
}


## @title ANFIS layer adder
## @description
## To add a layer for the given ANFIS model
## @param anfis The given ANFIS model
## @param name The layer name
## @return The ANFIS model with added layer
## @details This function is not designed for external use.
## @examples
## fis <- anfis.tipper()
## anfis <- list(andMethod=fis$andMethod, orMethod=fis$orMethod, impMethod=fis$impMethod, aggMethod=fis$aggMethod,
##                        defuzzMethod=fis$defuzzMethod, rule=fis$rule, layer=NULL)
## anfis <- anfis.addlayer(anfis, "inputs")
## @author Chao Chen
## @export

anfis.addlayer <- function(anfis, name) {
    anfis$layer <- append(anfis$layer, list(list(name = name, node = NULL)))
    anfis
}


## @title ANFIS node adder
## @description
## To add a node to the specified layer of a given ANFIS model.
## @param anfis The given ANFIS model
## @param idx The index of the layer where node is going to be added.
## @param value node value
## @param dedo The derivative of the error to the output of this node.
## @param fan.in The indexes of input nodes from previous layer.
## @param fan.out The indexes of output nodes in next layer.
## @param mf.type membership function type, if this is a node of fuzzy membership function.
## @param mf.params The parameters for the membership function
## @param range The universe of discourse for the membership function
## @param mf.name The name of membership function
## @return The ANFIS model with added node.
## @details This function is not designed for external use.
## @author Chao Chen

anfis.addnode <- function(anfis, idx, value = NULL, dedo = NULL, fan.in = NULL, fan.out = NULL, mf.type = NULL, mf.params = NULL, range = NULL, mf.name = NULL) {

    anfis$layer[[idx]]$node <- append(anfis$layer[[idx]]$node,
                                      list(list(value = value, dedo = dedo, fan.in = fan.in, fan.out = fan.out, mf.type = mf.type, mf.params = mf.params, range = range, mf.name = mf.name)))

    anfis
}


#' @title ANFIS evaluator
#' @description
#' To evaluate a ANFIS model with input data
#' @param anfis The given ANFIS model
#' @param input.stack The input data
#' @return The output of the anfis for given input data.
#' @examples
#' fis <- anfis.tipper()
#' anfis <- anfis.builder(fis)
#' data.num <- 5
#' input.num <- length(fis$input)
#' input.stack <- matrix(rnorm(data.num*input.num), ncol=input.num)
#' y <- matrix(rnorm(data.num))
#' data.trn <- cbind(input.stack, y)
#' anfis.eval(anfis, input.stack)
#' @author Chao Chen
#' @export

anfis.eval <- function(anfis, input.stack) {

    output.LI <- anfis.LI.eval(anfis, input.stack)
    output.L1 <- anfis.L1.eval(anfis, output.LI, input.stack)
    output.L2 <- anfis.L2.eval(anfis, output.L1)
    output.L4.mf <- anfis.L4.mf.eval(anfis, input.stack)
    output.L2.which <- anfis.L2.which(anfis, output.L2, output.L4.mf)
    output.L3 <- anfis.L3.eval(anfis, output.L2, output.L2.which)
    output.L4 <- anfis.L4.eval(output.L3, output.L4.mf)
    output.L5 <- anfis.L5.eval(output.L4)

    output.anfis <- output.L5

    output.anfis
}


#' @title ANFIS optimiser
#' @description
#' To optimise the performance of a given ANFIS model by learning the parameters in L1 and L4.
#' @param anfis The given ANFIS model
#' @param data.trn The input and output data pairs as training data
#' @param data.chk The input and output data pairs as checking (validation) data
#' @param epoch.total The total training epochs.
#' @param stepsize The initial stepsize
#' @param rate.inc increasing rate of the stepsize
#' @param rate.dec decrasing rate of the stepsize
#' @param method The learning algorithms for Layer 1 and Layer 4 respectively. default method=c("gradient", "lse")
#' @param err.log T or F, the flag indicate whether to save the error log.
#' @param online 0 -- batch; 1 -- online; 2 -- semi-online
#' @param lambda The forgetting rate for the LSE algorithm
#' @param opt.by To optimise the ANFIS model by: err.opt -- optimisation error; err.trn -- training error; err.chk -- checking (validation) error.
#' @param err.trn.fix T or F. When KM defuzzification is used for IT2 ANFIS, err.trn is not equal to err.opt. Hence, this flag is used for users to choose whether to fix this issue. The default value is set to T for the compatibility with previous built IT2 models. For T1 ANFIS, this flag can be set to F for speed improvement.
#' @return The optimised ANFIS model.
#' @examples
#' fis <- anfis.tipper()
#' anfis <- anfis.builder(fis)
#' data.num <- 5
#' input.num <- length(fis$input)
#' input.stack <- matrix(rnorm(data.num*input.num), ncol=input.num)
#' y <- matrix(rnorm(data.num))
#' data.trn <- cbind(input.stack, y)
#' anfis.eval(anfis, input.stack)
#' anfis.final <- anfis.optimise(anfis, data.trn, epoch.total=500,
#'                                  stepsize=0.01, rate.inc=1.1, rate.dec=0.9)
#' @author Chao Chen
#' @references
#' [1] C. Chen, R. John, J. Twycross, and J. M. Garibaldi, “An extended ANFIS architecture and its learning properties for type-1 and interval type-2 models,” in Proceedings IEEE International Conference on Fuzzy Systems, 2016, pp. 602–609. \cr
#' \doi{10.1109/FUZZ-IEEE.2016.7737742}
#'
#' [2] C. Chen, R. John, J. Twycross, and J. M. Garibaldi, “Type-1 and interval type-2 ANFIS: a comparison,” in Proceedings IEEE International Conference on Fuzzy Systems, 2017, pp. 1–6.  \cr
#' \doi{10.1109/FUZZ-IEEE.2017.8015555}
#' @export

anfis.optimise <- function(anfis, data.trn, data.chk = NULL, epoch.total = 100, stepsize = 0.1, rate.inc = 1.1, rate.dec = 0.9, method = c("gradient", "lse"), err.log = F, online = 0, lambda = 1, opt.by = "err.opt", err.trn.fix = T) {

    LI <- 1
    L1 <- 1 + 1
    L2 <- 1 + 2
    L3 <- 1 + 3
    L4 <- 1 + 4
    L5 <- 1 + 5

    if (method[2] == "lse") {
        for (i in 1:length(anfis$layer[[L4]]$node)) {
            if (anfis$layer[[L4]]$node[[i]]$mf.type != "linearmf") {
                stop("currently, the LSE method can only be applied to 'linearmf' for nodes in Layer 4")
            }
        }
    }

    if (length(online) != 1 || !is.element(online, c(0, 1, 2))) {
        stop("parameter online: 0 -- batch; 1 -- online; 2 -- semi-online")
    }

    input.num <- length(anfis$layer[[LI]]$node)
    output.num <- length(anfis$layer[[L5]]$node)

    data.num <- nrow(data.trn)
    input.stack <- input.stack.all <- data.trn[, 1:input.num]
    target <- target.all <- as.matrix(data.trn[, -(1:input.num)])
    scale.mase <- mean(abs(data.trn[, input.num] - data.trn[, (input.num + 1)]))

    anfis.optimum <- anfis
    err.min <- Inf
    err.min.rec <- NULL
    err.min.rec.tmp <- NULL
    err.opt.rec <- Inf
    epoch.last <- 1

    if (err.log) {
        # Add '.GlobalEnv$' to replace super assignment '<<-'
        .GlobalEnv$err.opt <- matrix(0, nrow = epoch.total, ncol = 8)
        .GlobalEnv$err.trn <- matrix(0, nrow = epoch.total, ncol = 8)
        .GlobalEnv$err.chk <- matrix(0, nrow = epoch.total, ncol = 8)
        .GlobalEnv$theta.L1 <- NULL
        .GlobalEnv$theta.L4 <- NULL
        colnames(.GlobalEnv$err.opt) <- colnames(.GlobalEnv$err.trn) <- colnames(.GlobalEnv$err.chk) <- c("MAE", "RMSE", "MASE", "MRAE", "GMRAE", "MAPE", "sMAPE", "uMbRAE")
    }

    epoch <- p <- 1
    if (online != 1) {
        p <- data.num
    } else {
        anfis.output.tmp <- target
        output.L3.tmp <- NULL
    }

    while (epoch <= epoch.total) {
        if (online == 1) {
            input.stack <- matrix(input.stack.all[p,], nrow = 1)
            target <- matrix(target.all[p,], nrow = 1)
        }

        output.LI <- anfis.LI.eval(anfis, input.stack)
        output.L1 <- anfis.L1.eval(anfis, output.LI, input.stack)
        output.L2 <- anfis.L2.eval(anfis, output.L1)
        output.L4.mf <- anfis.L4.mf.eval(anfis, input.stack)
        output.L2.which <- anfis.L2.which(anfis, output.L2, output.L4.mf)
        output.L3 <- anfis.L3.eval(anfis, output.L2, output.L2.which)

        if (method[2] == "lse") {
            anfis <- anfis.optimise.lse(anfis, output.L3, input.stack, target, online, lambda)
            output.L4.mf <- anfis.L4.mf.eval(anfis, input.stack)
            #            if(ncol(output.L3[[1]]) != 1) {
            #                output.L2.which <- anfis.L2.which(anfis, output.L2, output.L4.mf)
            #                output.L3 <- anfis.L3.eval(anfis, output.L2, output.L2.which)
            #            }
        }

        output.L4 <- anfis.L4.eval(output.L3, output.L4.mf)
        output.L5 <- anfis.L5.eval(output.L4)

        if (online == 1) {
            ## option 2
            #anfis.output.tmp[p] <- output.L5

            ## option 3
            output.L3.tmp <- lapply(1:length(output.L3), function(i) rbind(output.L3.tmp[[i]], output.L3[[i]]))
        }

        if (p == data.num) {
            if (online == 1) {
                ## option 1
                #anfis.output <- anfis.eval(anfis, input.stack.all)

                ## option 2
                #anfis.output <- anfis.output.tmp

                ## option 3
                output.L4.mf.tmp <- anfis.L4.mf.eval(anfis, input.stack.all)
                output.L4.tmp <- anfis.L4.eval(output.L3.tmp, output.L4.mf.tmp)
                output.L5.tmp <- anfis.L5.eval(output.L4.tmp)
                anfis.output <- output.L5.tmp
                output.L3.tmp <- NULL
            } else {
                anfis.output <- output.L5
            }

            if (err.log) {
                err.tmp <- fuzzyr.accuracy(anfis.output, data.trn[, (input.num + 1)], data.trn[, input.num], scale.mase)
                .GlobalEnv$err.opt[epoch,] <- err.tmp

                if (err.trn.fix) {
                    err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.trn[, 1:input.num]), data.trn[, (input.num + 1)], data.trn[, input.num], scale.mase)
                }
                .GlobalEnv$err.trn[epoch,] <- err.tmp

                if (!is.null(data.chk)) {
                    err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.chk[, 1:input.num]), data.chk[, (input.num + 1)], data.chk[, input.num], scale.mase)
                    .GlobalEnv$err.chk[epoch,] <- err.tmp
                }

                theta.L1.tmp <- NULL
                for (i in 1:length(anfis$layer[[L1]]$node)) {
                    theta.L1.tmp <- c(theta.L1.tmp, anfis$layer[[L1]]$node[[i]]$mf.params)
                }

                theta.L4.tmp <- NULL
                for (i in 1:length(anfis$layer[[L4]]$node)) {
                    theta.L4.tmp <- c(theta.L4.tmp, anfis$layer[[L4]]$node[[i]]$mf.params)
                }

                .GlobalEnv$theta.L1 <- rbind(.GlobalEnv$theta.L1, theta.L1.tmp)
                .GlobalEnv$theta.L4 <- rbind(.GlobalEnv$theta.L4, theta.L4.tmp)
            }

            err.opt.rmse <- sqrt(mean((anfis.output - target.all)^2))
            if (opt.by == "err.trn") {
                if (err.log) {
                    err.epoch <- .GlobalEnv$err.trn[epoch, 2]
                } else {
                    if (err.trn.fix) {
                        err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.trn[, 1:input.num]), data.trn[, (input.num + 1)], data.trn[, input.num], scale.mase)
                        err.epoch <- err.tmp[2]
                    } else {
                        err.epoch <- err.opt.rmse
                    }
                }
            } else if (opt.by == "err.chk") {
                if (err.log) {
                    err.epoch <- .GlobalEnv$err.chk[epoch, 2]
                } else {
                    err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.chk[, 1:input.num]), data.chk[, (input.num + 1)], data.chk[, input.num], scale.mase)
                    err.epoch <- err.tmp[2]
                }
            } else {
                err.epoch <- err.opt.rmse
            }

            #cat(" epoch[", epoch, "], err=", err.epoch, "\n")
            if (err.epoch < err.min) {
                err.min <- err.epoch
                err.min.rec.tmp <- c(err.min.rec.tmp, err.min)
                if (length(err.min.rec.tmp) == 1000) {
                    err.min.rec <- c(err.min.rec, err.min.rec.tmp)
                    err.min.rec.tmp <- NULL
                }
                anfis.optimum <- anfis
            }

            ## to update step_size

            #cat("epoch - epoch.last: ", epoch - epoch.last, "\n")
            #cat("stepsize: ", stepsize, "\n")

            err.opt.rec <- c(err.opt.rec, err.opt.rmse)
            if (epoch - epoch.last >= 4) {
                err.opt.rec <- tail(err.opt.rec, 5)
                decflag <- -sign(diff(err.opt.rec))

                decflag.sum <- sum(decflag)
                if (decflag.sum == 4) {
                    stepsize <- stepsize * rate.inc
                    epoch.last <- epoch
                } else if (decflag.sum <= 0) {
                    #} else if( identical(decflag, c(-1,1,-1,1)) ) {
                    stepsize = stepsize * rate.dec
                    epoch.last = epoch
                }
            }
        }


        ## to optimise params in L1 and/or L4 by gradient decent method
        de.do5 <- anfis.dE.dO5(output.L5, target)
        do5.do4 <- anfis.dO5.dO4(output.L4)
        de.do4 <- anfis.dE.dO4(anfis, de.do5, do5.do4)
        de.dp4 <- anfis.dE.dP4(anfis, de.do4, output.L3, input.stack)
        do4.do3 <- anfis.dO4.dO3(output.L4, output.L4.mf)
        de.do3 <- anfis.dE.dO3(de.do4, do4.do3, output.L3)
        do3.do2 <- anfis.dO3.dO2(anfis, output.L2, output.L2.which)
        de.do2 <- anfis.dE.dO2(de.do3, do3.do2)
        do2.do1 <- anfis.dO2.dO1(anfis, output.L2, output.L1)
        de.do1 <- anfis.dE.dO1(anfis, output.L1, de.do2, do2.do1)
        de.dp1 <- anfis.dE.dP1(anfis, de.do1, input.stack)

        ## TODO: to decide whether to update params according to the configuration
        if (method[1] == "gradient") {
            anfis <- anfis.optimise.gradient(anfis, L1, de.dp1, stepsize)
        }

        if (method[2] == "gradient" || method[2] == "both") {
            anfis <- anfis.optimise.gradient(anfis, L4, de.dp4, stepsize)
        }

        if (online == 1 && p < data.num) {
            p <- p + 1
        } else {
            if (online == 1) { p <- 1 }
            epoch <- epoch + 1
        }
    }
    err.min.rec <- c(err.min.rec, err.min.rec.tmp)
    anfis.optimum <- append(anfis.optimum, list(RMSE_MIN_ROUTE = err.min.rec))
    show(err.min.rec)
    anfis.optimum
}


## @title The gradient algorithm for ANFIS optimisation
## @description
## To update the mf.params in Layer 1 or Layer 4 by the gradient decent methed
## @param anfis The given ANFIS model
## @param L The index of the layer in which the parameters to be optimised
## @param de.dp The derivatives of the output error to the parameters
## @param stepsize the stepsize for updating the parameters
## @return The ANFIS model with updated parameters
## @details This function is not designed for external use.
## @author Chao Chen
## @references

anfis.optimise.gradient <- function(anfis, L, de.dp, stepsize) {

    de.dp.sum <- lapply(de.dp, apply, 2, sum)
    len = sqrt(sum(unlist(de.dp.sum)^2))

    if (len != 0) {
        eta <- stepsize / len
        for (i in 1:length(de.dp.sum)) {
            mf.params <- anfis$layer[[L]]$node[[i]]$mf.params - eta * de.dp.sum[[i]]
            mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
            if ((mf.type == "it2gbellmf" && sign(mf.params[3]) * abs(mf.params[1]) > sign(mf.params[3]) * abs(mf.params[2]))
                || (mf.type == "itlinearmf" && mf.params[1] > mf.params[2])) {
                #                    mf.params[1] <- mf.params[2]
                mf.params[c(1, 2)] <- mf.params[c(2, 1)]
            }
            anfis$layer[[L]]$node[[i]]$mf.params <- mf.params
        }
    }

    anfis
}


## @title The LSE algorithm for ANFIS optimisation
## @description
## to update the mf.params in L4 by the LSE algorithm
## @param anfis The given ANFIS model
## @param output.L3 The output of nodes from Layer 3
## @param input.stack The input data (training data)
## @param target The output data (training data)
## @param online 0 -- batch; 1,2 -- online
## @param lambda the forgetting rate
## @return The ANFIS model with updated parameters
## @details This function is not designed for external use.
## @author Chao Chen
## @references

anfis.optimise.lse <- function(anfis, output.L3, input.stack, target, online = 0, lambda = 1) {

    #    if(mean(sapply(output.L3, ncol)) != 1) {
    #        stop("lse method here can only be used with type-1 anfis")
    #    }

    L <- 4 + 1
    node.num <- length(anfis$layer[[L]]$node)

    input.x <- cbind(1, input.stack)
    data.num <- nrow(input.x)

    if (ncol(output.L3[[1]]) == 1) {
        A <- matrix(sapply(output.L3, function(w) c(w) * input.x), nrow = data.num)
    } else {
        A <- matrix(sapply(output.L3, function(w) apply(w, 1, mean) * input.x), nrow = data.num)
    }

    B <- target

    if (!online || is.null(anfis$S)) {
        #if(is.null(anfis$S)) {
        S <- 999999 * diag(ncol(A))
        if (!online) {
            theta <- matrix(rep(0, node.num * ncol(input.x)), ncol = 1)
        } else {
            theta <- matrix(sapply(anfis$layer[[L]]$node, function(x) x$mf.params), ncol = 1)
        }
    } else {
        S <- anfis$S
        theta <- matrix(sapply(anfis$layer[[L]]$node, function(x) x$mf.params), ncol = 1)
    }

    for (i in 1:data.num) {
        ai = A[i,]
        bi = B[i,]

        ## (S %*% ai) %*% (t(ai) %*% S) is much more efficient than (S %*% ai %*% t(ai) %*% S)
        S <- (1 / lambda) * (S - (S %*% ai) %*% (t(ai) %*% S) / c(lambda + (t(ai) %*% S %*% ai)))
        theta <- theta + S %*% ai %*% (t(bi) - t(ai) %*% theta)
    }

    theta <- matrix(theta, ncol = node.num)
    for (i in 1:node.num) {
        anfis$layer[[L]]$node[[i]]$mf.params <- theta[, i]
    }

    if (online) {
        anfis$S <- S
    }

    anfis
}


#' @title The evaluator for nodes in Layer I
#' @description
#' To evaluate the input Layer (LI) of anfis
#' @param anfis The given ANFIS model
#' @param input.stack The input data
#' @return The output of nodes in Layer I
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.LI.eval <- function(anfis, input.stack) {

    L <- 1
    node.num <- length(anfis$layer[[L]]$node)

    output.LI <- list()
    for (i in 1:node.num) {
        mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
        mf.params <- anfis$layer[[L]]$node[[i]]$mf.params
        FUN <- match.fun(mf.type)
        x <- matrix(input.stack[, i], ncol = 1)
        output.LI <- append(output.LI, list(apply(x, 1, FUN, mf.params)))
    }

    output.LI
}


#' @title The evaluator for nodes in Layer 1
#' @description
#' To evaluate the antecedent layer (L1) of anfis
#' @param anfis The given ANFIS model
#' @param output.LI The output of nodes in Layer I
#' @param input.stack The input data
#' @return The output of nodes in Layer 1
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.L1.eval <- function(anfis, output.LI, input.stack) {

    L <- 1 + 1
    node.num <- length(anfis$layer[[L]]$node)

    output.L1 <- list()
    operator <- anfis$andMethod
    for (i in 1:node.num) {
        ante.mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
        ante.mf.params <- anfis$layer[[L]]$node[[i]]$mf.params
        ante.mf <- genmf(ante.mf.type, ante.mf.params)

        fan.in <- anfis$layer[[L]]$node[[i]]$fan.in
        x.mf <- output.LI[[fan.in]]
        x.mf.type <- anfis$layer[[L - 1]]$node[[fan.in]]$mf.type
        if (x.mf.type == "singleton.fuzzification") {
            x.lower <- input.stack[, fan.in]
            x.upper <- x.lower
        } else {
            x.range <- anfis$layer[[L - 1]]$node[[fan.in]]$range
            x.lower <- x.range[1]
            x.upper <- x.range[2]
        }

        firing.strength <- t(mapply(fuzzy.firing, operator, x.mf, list(ante.mf), x.lower, x.upper, SIMPLIFY = T))
        firing.strength <- matrix(firing.strength, nrow = nrow(input.stack))
        rownames(firing.strength) <- NULL
        output.L1 <- append(output.L1, list(firing.strength))
    }

    output.L1
}


#' @title The evaluator for nodes in Layer 2
#' @description
#' To evaluate the nodes in Layer 2 of the given ANFIS model
#' @param anfis The given ANFIS model
#' @param output.L1 The output of nodes in Layer 1
#' @return The output of nodes in Layer 2
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.L2.eval <- function(anfis, output.L1) {

    L <- 2 + 1
    node.num <- length(anfis$layer[[L]]$node)
    n <- nrow(output.L1[[1]])

    output.L2 <- list()
    connectors <- c(anfis$andMethod, anfis$orMethod)
    ## This mapvalues is for the alternative calculation below
    connectors <- mapvalues(connectors, from = c("prod", "max", "min"), to = c("multiply", "pmax", "pmin"), warn_missing = F)
    rule.connector <- connectors[anfis$rule[, ncol(anfis$rule)]]

    for (i in 1:node.num) {
        FUN <- match.fun(rule.connector[i])
        fan.in <- anfis$layer[[L]]$node[[i]]$fan.in

        #rule.strength <- output.L1[[fan.in[1]]]
        #for( j in fan.in[-1] ) {
        #    rule.strength <- matrix(mapply(FUN, rule.strength, output.L1[[j]]), nrow=n)
        #}

        ## an alternative approach for calculating the above rule.strength
        rule.strength <- Reduce(FUN, output.L1[fan.in])
        if (is.vector(rule.strength)) {
            ## required when pmin or pmax is the operator
            rule.strength <- matrix(rule.strength, nrow = n)
        }

        ## to use NieTan as an alternative defuzzification approach of KM algorithms
        if (anfis$defuzzMethod == "NT" && ncol(rule.strength) == 2) {
            rule.strength <- matrix(apply(rule.strength, 1, mean))
        }

        output.L2 <- append(output.L2, list(rule.strength))
    }

    output.L2
}


#' @title L2.which
#' @description
#' To determin which output (w.lower, w.upper) to be used by the ekm algorithm
#' @param anfis The given ANFIS model
#' @param output.L2 The output of nodes in Layer 2
#' @param output.L4.mf The linear membership grades of nodes in Layer 4
#' @return A list of matrix indicating which output (w.lower, w.upper) in layer 2 should be used by the ekm algorithm
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.L2.which <- function(anfis, output.L2, output.L4.mf) {

    L <- 2 + 1
    node.num <- length(anfis$layer[[L]]$node)
    n <- nrow(output.L2[[1]])

    col.num <- sapply(output.L2, ncol)

    if (mean(col.num) == 1) {
        output.L2.which <- list(lapply(1:node.num, function(x) matrix(rep(TRUE, n))))
    } else {
        output.L2.which.lower <- vector("list", node.num)
        output.L2.which.upper <- vector("list", node.num)
        for (i in 1:n) {
            w.L2 <- t(sapply(output.L2, function(x) c(x[i, 1], x[i, ncol(x)])))
            f.L4 <- t(sapply(output.L4.mf, function(x) c(x[i, 1], x[i, ncol(x)])))
            w.lower <- w.L2[, 1]
            w.upper <- w.L2[, 2]
            f.lower <- f.L4[, 1]
            f.upper <- f.L4[, 2]
            w.which.lower <- km.da(w.lower, w.upper, f.lower, maximum = F, w.which = T)
            w.which.upper <- km.da(w.lower, w.upper, f.upper, maximum = T, w.which = T)

            output.L2.which.lower <- lapply(1:node.num, function(idx) rbind(output.L2.which.lower[[idx]], if (col.num[idx] == 1) TRUE else w.which.lower[idx,]))
            output.L2.which.upper <- lapply(1:node.num, function(idx) rbind(output.L2.which.upper[[idx]], if (col.num[idx] == 1) TRUE else w.which.upper[idx,]))
        }
        output.L2.which <- list(output.L2.which.lower, output.L2.which.upper)
    }

    output.L2.which
}


#' @title The evaluator for nodes in Layer 3
#' @description
#' To evaluate the nodes in Layer 3 of the given ANFIS model
#' @param anfis The given ANFIS model
#' @param output.L2 The output of nodes in Layer 2
#' @param output.L2.which A list of matrix indicating which output (w.lower, w.upper) in layer 2 should be used by the ekm algorithm
#' @return The output of nodes in Layer 3
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.L3.eval <- function(anfis, output.L2, output.L2.which) {

    L <- 3 + 1
    node.num <- length(anfis$layer[[L]]$node)

    output.L3 <- list()

    if (mean(sapply(output.L2, ncol)) == 1) {
        w.tmp <- matrix(unlist(output.L2), nrow = nrow(output.L2[[1]]))
        w.sum <- as.matrix(rowSums(w.tmp))
        for (i in 1:node.num) {
            w.normal <- ifelse(w.sum != 0, w.tmp[, i] / w.sum, 1 / node.num)
            output.L3 <- append(output.L3, list(w.normal))
        }
    } else {
        n <- nrow(output.L2[[1]])
        output.L3 <- vector("list", node.num)
        for (i in 1:n) {
            w.L2 <- unlist(sapply(output.L2, function(x) x[i,]))
            w.which.lower <- unlist(sapply(output.L2.which[[1]], function(x) x[i,]))
            w.which.upper <- unlist(sapply(output.L2.which[[2]], function(x) x[i,]))

            w.normal.lower <- if (sum(w.L2[w.which.lower]) != 0) w.L2[w.which.lower] / sum(w.L2[w.which.lower]) else 1 / node.num
            w.normal.upper <- if (sum(w.L2[w.which.upper]) != 0) w.L2[w.which.upper] / sum(w.L2[w.which.upper]) else 1 / node.num
            output.L3 <- lapply(1:node.num, function(idx) rbind(output.L3[[idx]], c(w.normal.lower[[idx]], w.normal.upper[[idx]])))
        }
    }

    output.L3
}


#' @title The evaluator for membership functions of nodes in Layer 1
#' @description
#' To evaluate the membership functions of nodes in Layer 4
#' @param anfis The given ANFIS model
#' @param input.stack The input data
#' @return The membership grades of the membership functions of nodes in Layer 4
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.L4.mf.eval <- function(anfis, input.stack) {

    L <- 4 + 1
    node.num <- length(anfis$layer[[L]]$node)

    data.num <- nrow(input.stack)
    input.x <- cbind(1, input.stack)
    output.L4.mf <- list()

    for (i in 1:node.num) {
        conse.mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
        conse.mf.params <- anfis$layer[[L]]$node[[i]]$mf.params
        conse.mf <- genmf(conse.mf.type, conse.mf.params)
        output.L4.mf <- append(output.L4.mf, list(matrix(evalmf(input.x, conse.mf), nrow = data.num)))
    }

    output.L4.mf
}


#' @title The evaluator for nodes in Layer 4
#' @description
#' To evaluate the nodes in Layer 4
#' @param output.L3 The output of nodes in Layer 3
#' @param output.L4.mf The membership grades of the membership functions of nodes in Layer 4
#' @return The output of nodes in Layer 4
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export

anfis.L4.eval <- function(output.L3, output.L4.mf) {

    output.L4 <- mapply(function(x, y) matrix(mapply("*", x, y), nrow = nrow(x)), output.L3, output.L4.mf, SIMPLIFY = F)
    output.L4
}


#' @title The evaluator for nodes in Layer 5
#' @description
#' To evaluate the nodes in Layer 5
#' @param output.L4 The output of nodes in Layer 4
#' @return The output of nodes in Layer 5
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' See the source code of \code{\link{anfis.eval}} for usage.
#' @author Chao Chen
#' @export
anfis.L5.eval <- function(output.L4) {

    col.L4 <- mean(sapply(output.L4, ncol))

    if (col.L4 != 1) {
        output.L4 <- lapply(output.L4, function(x) { cbind(x[, 1], x[, ncol(x)]) })
    }

    n.row <- nrow(output.L4[[1]])
    n.col <- ncol(output.L4[[1]])
    n <- length(output.L4)

    output.L4 <- array(unlist(output.L4), dim = c(n.row, n.col, n))
    output.L5 <- apply(output.L4, c(1, 2), sum)

    output.L5 <- matrix(apply(output.L5, 1, mean), nrow = n.row)

    output.L5
}


#' @title anfis.dE.dO5
#' @description
#' To calculate the derivatives of output error with respect to output.L5.
#' NOTE: currently, only single output in L5 is supported
#' @param output.L5 the model outputs
#' @param y the target outputs
#' @return The derivatives of output error with respect to output.L5
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dO5 <- function(output.L5, y) {
    de.do5 <- -2 * (y - output.L5)
    de.do5
}


#' @title anfis.dO5.dO4
#' @description
#' To calculate the derivatives of output.L5 with respect to output.L4.
#' NOTE: currently, only single output in L5 is supported
#' @param output.L4 The output of nodes in Layer 4.
#' @return The derivatives of output.L5 with respect to output.L4.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dO5.dO4 <- function(output.L4) {
    do5.do4 <- lapply(output.L4, function(x) matrix(1 / ncol(x), nrow(x), ncol(x)))
    do5.do4
}


#' @title anfis.dE.dO4
#' @description
#' to calculate the derivatives of output error with respect to output.L4.
#' @param anfis The given ANFIS model
#' @param de.do5 The derivatives of output error with respect to output.L5
#' @param do5.do4 The derivatives of output.L5 with respect to output.L4.
#' @return The derivatives of output error with respect to output.L4.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dO4 <- function(anfis, de.do5, do5.do4) {
    de.do4 <- lapply(do5.do4, function(x) x * c(de.do5))
    de.do4
}


#' @title anfis.dE.dP4
#' @description
#' To calculate the derivatives of output error with respect to parameters in Layer 4.
#' @param anfis The given ANFIS model
#' @param de.do4 The derivatives of output error with respect to output.L4
#' @param output.L3 The output of nodes in Layer 3
#' @param input.stack The input data pairs.
#' @return The derivatives of output error with respect to parameters in Layer 4.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dP4 <- function(anfis, de.do4, output.L3, input.stack) {

    L <- 4 + 1
    node.num <- length(anfis$layer[[L]]$node)

    n <- nrow(input.stack)
    input.x <- cbind(1, input.stack)

    de.dp4 <- list()
    for (i in 1:node.num) {
        conse.mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
        #conse.mf.params <- anfis$layer[[L]]$node[[i]]$mf.params

        tmp <- matrix(mapply("*", output.L3[[i]], de.do4[[i]]), nrow = n)

        if (conse.mf.type == "linearmf") {
            de.dp4 <- append(de.dp4, list(apply(tmp, 1, sum) * input.x))
        } else if (conse.mf.type == "itlinearmf") {
            de.dp4 <- append(de.dp4, list(cbind(tmp, apply(tmp, 1, sum) * input.stack)))
        } else {
            stop("dE.dP4 error: consequent member function type not supported")
        }
    }

    de.dp4
}


#' @title anfis.dO4.dO3
#' @description
#' To calculate the derivatives of output.L4 with respect to output.L3.
#' @param output.L4 The output of nodes in Layer 4
#' @param output.L4.mf The membership grades of the membership functions of nodes in Layer 4
#' @return The derivatives of output.L4 with respect to output.L3.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dO4.dO3 <- function(output.L4, output.L4.mf) {
    do4.do3 <- mapply(function(x, y) matrix(mapply("*", matrix(1, nrow(x), ncol(x)), y), nrow = nrow(x)), output.L4, output.L4.mf, SIMPLIFY = F)
    do4.do3
}


#' @title anfis.dE.dO3
#' @description
#' to calculate the derivatives of output error with respect to output.L3.
#' @param de.do4 The derivatives of output error with respect to output.L4
#' @param do4.do3 The derivatives of output.L4 with respect to output.L3.
#' @param output.L3 The output of nodes in Layer 3.
#' @return The derivatives of output error with respect to output.L3.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dO3 <- function(de.do4, do4.do3, output.L3) {
    de.do3 <- mapply(function(x, y, z) { if (ncol(z) == 1) matrix(apply(x * y, 1, sum)) else x * y }, de.do4, do4.do3, output.L3, SIMPLIFY = F)
    de.do3
}


#' @title anfis.dO3.dO2
#' @description
#' To calculate the derivatives of output.L3 with respect to output.L2.
#' @param anfis The given ANFIS model
#' @param output.L2 The output of nodes in Layer 2
#' @param output.L2.which A list of matrix indicating which output (w.lower, w.upper) in layer 2 should be used by the ekm algorithm
#' @return The derivatives of output.L3 with respect to output.L2. do3.left[j].do2[i] <- do3.do2[[i]][[1]][[j]]
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dO3.dO2 <- function(anfis, output.L2, output.L2.which) {

    L <- 2 + 1
    node.num <- length(anfis$layer[[L]]$node)
    row.num <- nrow(output.L2[[1]])

    output.L2.sum <- list()
    output.L2.w <- list()
    for (i in 1:length(output.L2.which)) {
        tmp1 <- matrix(unlist(output.L2), nrow = row.num)
        tmp2 <- matrix(unlist(output.L2.which[[i]]), nrow = row.num)
        tmp3 <- t(sapply(1:row.num, function(idx) tmp1[idx,][tmp2[idx,]]))
        output.L2.w <- append(output.L2.w, list(tmp3))
        output.L2.sum <- append(output.L2.sum, list(matrix(apply(tmp3, 1, sum))))
    }

    do3.do2 <- list()
    for (i in 1:node.num) {
        #fan.out <- anfis$layer[[L]]$node[[i]]$fan.out
        do3.do2.tmp.k <- list()
        for (k in 1:length(output.L2.which)) {
            do3.do2.tmp.j <- list()
            for (j in 1:node.num) {
                if (i == j) {
                    # output.L2.sum[[k]]^2 may be zero if output.L2.sum[[k]] is too small
                    #tmp <- (output.L2.sum[[k]] - output.L2.w[[k]][,i]) / output.L2.sum[[k]]^2
                    tmp <- (output.L2.sum[[k]] - output.L2.w[[k]][, i]) /
                        output.L2.sum[[k]] /
                        output.L2.sum[[k]]
                    tmp[is.nan(tmp)] = 1 / node.num
                } else {
                    #tmp <- -output.L2.w[[k]][,j] / output.L2.sum[[k]]^2
                    tmp <- -output.L2.w[[k]][, j] /
                        output.L2.sum[[k]] /
                        output.L2.sum[[k]]
                    tmp[is.nan(tmp)] = 1 / node.num
                }
                do3.do2.tmp.j <- append(do3.do2.tmp.j, list(c(tmp) * output.L2.which[[k]][[i]]))
            }
            do3.do2.tmp.k <- append(do3.do2.tmp.k, list(do3.do2.tmp.j))
        }

        do3.do2 <- append(do3.do2, list(do3.do2.tmp.k))
    }

    do3.do2
}


#' @title anfis.dE.dO2
#' @description
#' to calculate the derivatives of output error with respect to output.L2.
#' @param de.do3 The derivatives of output error with respect to output.L3
#' @param do3.do2 The derivatives of output.L3 with respect to output.L2.
#' @return The derivatives of output error with respect to output.L2.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dO2 <- function(de.do3, do3.do2) {

    node.num <- length(do3.do2)

    de.do2 <- list()
    for (i in 1:node.num) {
        row.num <- nrow(do3.do2[[i]][[1]][[1]])
        col.num <- ncol(do3.do2[[i]][[1]][[1]])
        de.do2.tmp.i <- matrix(0, row.num, col.num)

        for (j in 1:length(do3.do2[[i]])) {
            for (k in 1:length(do3.do2[[i]][[j]])) {
                de.do2.tmp.i <- de.do2.tmp.i + de.do3[[k]][, j] * do3.do2[[i]][[j]][[k]]
            }
        }

        colnames(de.do2.tmp.i) <- NULL
        de.do2 <- append(de.do2, list(de.do2.tmp.i))
    }

    de.do2
}


#' @title anfis.dO2.dO1
#' @description
#' To calculate the derivatives of output.L2 with respect to output.L1.
#' @param anfis The given ANFIS model
#' @param output.L2 The output of nodes in Layer 2
#' @param output.L1 The output of nodes in Layer 1
#' @return The derivatives of output.L2 with respect to output.L1. do2[j].do1[i] <- do2.do1[[i]][[which(fan.out==j)]]
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dO2.dO1 <- function(anfis, output.L2, output.L1) {

    L <- 1 + 1
    node.num <- length(output.L1)
    row.num <- nrow(output.L1[[1]])

    do2.do1 <- list()
    for (i in 1:node.num) {
        fan.out <- anfis$layer[[L]]$node[[i]]$fan.out

        do2.do1.i <- list()
        for (j in fan.out) {
            #do2.do1.j <- matrix(mapply(function(x, y) ifelse(y!=0, x/y, 0), output.L2[[j]], output.L1[[i]]), nrow=row.num)

            ## An alternative approach to calculate do2.do1.j above
            fan.in <- anfis$layer[[L + 1]]$node[[j]]$fan.in
            do2.do1.in <- setdiff(fan.in, i)
            if (length(do2.do1.in) == 0) {
                do2.do1.j <- matrix(1, nrow = row.num, ncol = ncol(output.L1[[i]]))
            } else {
                do2.do1.j <- Reduce(multiply, output.L1[do2.do1.in])
                if (ncol(output.L1[[i]]) == 2 && ncol(do2.do1.j) == 1) {
                    do2.do1.j <- do2.do1.j[, c(1, 1)]
                }
            }

            ## Added for IT2 ANFIS based on Nie-Tan's approach as the defuzzification method
            if (anfis$defuzzMethod == "NT" && ncol(output.L1[[i]]) == 2) {
                do2.do1.j <- do2.do1.j / 2
            }

            do2.do1.i <- append(do2.do1.i, list(do2.do1.j))
        }
        do2.do1 <- append(do2.do1, list(do2.do1.i))
    }

    do2.do1
}


#' @title anfis.dE.dO1
#' @description
#' to calculate the derivatives of output error with respect to output.L1.
#' @param anfis The given ANFIS model
#' @param output.L1 The output of nodes in Layer 1
#' @param de.do2 The derivatives of output error with respect to output.L2
#' @param do2.do1 The derivatives of output.L2 with respect to output.L1.
#' @return The derivatives of output error with respect to output.L1.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dO1 <- function(anfis, output.L1, de.do2, do2.do1) {

    L <- 1 + 1
    node.num <- length(do2.do1)
    row.num <- nrow(output.L1[[1]])

    de.do1 <- list()
    for (i in 1:node.num) {
        col.num <- ncol(output.L1[[i]])
        fan.out <- anfis$layer[[L]]$node[[i]]$fan.out
        de.do1.tmp.i <- matrix(0, row.num, col.num)

        if (length(fan.out) > 0) {
            for (j in 1:length(fan.out)) {
                de.do1.tmp.j <- multiply(de.do2[[fan.out[j]]], do2.do1[[i]][[j]])
                if (ncol(de.do1.tmp.j) > col.num) {
                    de.do1.tmp.j <- apply(de.do1.tmp.j, 1, sum)
                }
                de.do1.tmp.i <- de.do1.tmp.i + de.do1.tmp.j
            }
        }

        de.do1 <- append(de.do1, list(de.do1.tmp.i))
    }

    de.do1
}


#' @title anfis.dE.dP1
#' @description
#' To calculate the derivatives of output error with respect to parameters in Layer 1.
#' @param anfis The given ANFIS model
#' @param de.do1 The derivatives of output error with respect to output.L1
#' @param input.stack The input data pairs.
#' @return The derivatives of output error with respect to parameters in Layer 1.
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dP1 <- function(anfis, de.do1, input.stack) {

    L <- 1 + 1
    node.num <- length(anfis$layer[[L]]$node)

    de.dp1 <- list()
    for (i in 1:node.num) {
        fan.in <- anfis$layer[[L]]$node[[i]]$fan.in
        x.mf.type <- anfis$layer[[L - 1]]$node[[fan.in]]$mf.type
        if (x.mf.type != "singleton.fuzzification") {
            stop("dE.dP1: only singleton fuzzification is supported")
        }

        ante.mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
        ante.mf.params <- anfis$layer[[L]]$node[[i]]$mf.params

        FUN <- match.fun(paste("anfis.dE.dP1.", ante.mf.type, sep = ""))
        de.dp1.tmp <- FUN(de.do1[[i]], input.stack[, fan.in], ante.mf.params)
        de.dp1 <- append(de.dp1, list(de.dp1.tmp))
    }

    de.dp1
}


#' @title anfis.dE.dP1.gbellmf
#' @description
#' To calculate the derivatives of E versus mf.params.L1 for gbellmf: 1 / ( 1 + (((x - c)/a)^2)^b)
#' NOTE: only singleton fuzzification is supported
#' @param de.do1 The derivatives of output error with respect to output.L1
#' @param x The crisp input
#' @param mf.params parameters for membership functions
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dP1.gbellmf <- function(de.do1, x, mf.params) {
    de.dp1.gbellmf <- c(de.do1) * anfis.dMF.dP.gbellmf(x, mf.params)
}


#' @title anfis.dE.dP1.it2gbellmf
#' @description
#' to calculate the derivatives of E versus mf.params.L1 for it2gbellmf
#' NOTE: only singleton fuzzification is supported
#' @param de.do1 The derivatives of output error with respect to output.L1
#' @param x The crisp input
#' @param mf.params parameters for membership functions
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dE.dP1.it2gbellmf <- function(de.do1, x, mf.params) {
    de.dp1.it2gbellmf.lower <- c(de.do1[, 1]) * anfis.dMF.dP.gbellmf(x, mf.params[-2])
    de.dp1.it2gbellmf.upper <- c(de.do1[, 2]) * anfis.dMF.dP.gbellmf(x, mf.params[-1])
    de.dp1.gbellmf <- cbind(as.matrix(de.dp1.it2gbellmf.lower[, 1]), as.matrix(de.dp1.it2gbellmf.upper[, 1]),
                            matrix(de.dp1.it2gbellmf.lower[, -1] + de.dp1.it2gbellmf.upper[, -1], ncol = 2))
}


#' @title anfis.dMF.dP.gbellmf
#' @description
#' to calculate the derivatives of membership grades with respect to its parameters
#' @param x The crisp input
#' @param mf.params parameters for membership functions
#' @details This function is not recommended for external use, but can be used for debugging or learning.
#' @author Chao Chen
#' @export

anfis.dMF.dP.gbellmf <- function(x, mf.params) {

    a = mf.params[1]
    b = mf.params[2]
    c = mf.params[3]

    denom = 0

    if (a == 0) {
        stop("anfis.dE.dP1.gbellmf: a == 0 founded!")
    }

    tmp1 = ((x - c) / a)^2
    tmp2 = tmp1^b
    denom = (1 + tmp2)^2

    dmf.dp = NULL
    # partial mf to partial a
    #tmp = (2*b*tmp2)/(a*denom)
    tmp = (2 * b * tmp2) / (1 + tmp2) / a / (1 + tmp2)
    tmp[tmp2 == Inf] = 0
    tmp[which(is.nan(tmp))] = 0
    tmp[tmp == Inf] = 0
    dmf.dp = cbind(dmf.dp, tmp)

    # partial mf to partial b
    #tmp = -1*log(tmp1)*tmp2/denom
    tmp = -1 * log(tmp1) * tmp2 / (1 + tmp2) / (1 + tmp2)
    tmp[tmp1 == 0] = 0
    tmp[which(is.nan(tmp))] = 0
    tmp[tmp == Inf] = 0
    dmf.dp = cbind(dmf.dp, tmp)

    # partial mf to partial c
    #tmp = (2*b*tmp2)/((x - c)*(denom))
    tmp = (2 * b * tmp2) /
        (x - c) /
        (1 + tmp2) /
        (1 + tmp2)
    tmp[x == c] = 0
    tmp[which(is.nan(tmp))] = 0
    tmp[tmp == Inf] = 0
    dmf.dp = cbind(dmf.dp, tmp)

    colnames(dmf.dp) <- NULL
    dmf.dp
}
