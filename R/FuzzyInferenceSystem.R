### FuzzyR - Fuzzy Inference System

# add by Tajul (to handle 'no visible binding for global variable')
globalVariables(c('D_x', 'D_y', 'D_out', 'genRule', 'perturb', 'NumInputs', 'Range', 'NumMFs',
                  'mfParams', 'NumOutputs', 'NumRules', 'ruleList', 'err.opt', 'err.trn', 'err.chk', 'myaccuracy',
                  'theta.L1', 'ekm', 'theta.L4', 'optimise', 'persp', 'rainbow', 'plot', 'axis', 'mtext', 'polygon',
                  'gray', 'lines', 'text', 'mtext', 'pdf', 'dev.off', 'show', 'tail', 'OUT_COUPLE', 'pushViewport',
                  'plotViewport', 'popViewport', 'RULE_CONS', 'IN_COUPLE', 'drawAnteEval', 'drawConsEval',
                  'drawRuleEval', 'viewport', 'drawRuleAggr', 'trCentroid', 'plot_it2', 'trCOS', 'iterative_method',
                  'drawAllSteps', 'rangex', 'OUT_RULE_CONS'))
#end

#' @import grid
#' @importFrom stats integrate
#' @importFrom graphics abline legend points title par
NULL


#' Create a fis using newfis function
#'
#' Creates a fis object.
#'
#' @param fisName String representing the fis name.
#' @param fisType Type of the fis, default is 'mamdani'.
#' @param mfType Type of membership functions, 't1' or 'it2'
#' @param andMethod The AND method for the fis, default is 'min'.
#' @param orMethod The OR method for the fis, default is 'max'.
#' @param impMethod The implication method for the fis, default is 'min'.
#' @param aggMethod The aggregation method for the fis, default is 'max'.
#' @param defuzzMethod The defuzzification method for the fis, default is 'centroid'.
#' @return A new fis structure.
#' @examples
#' fis <- newfis("fisName")
#' @export
newfis <- function(fisName, fisType = "mamdani", mfType = "t1",
                   andMethod = "min", orMethod = "max", impMethod = "min", aggMethod = "max",
                   defuzzMethod = "centroid") {
    fis <- list(name = fisName, type = fisType, mfType = mfType,
                andMethod = andMethod, orMethod = orMethod,
                impMethod = impMethod, aggMethod = aggMethod,
                defuzzMethod = defuzzMethod,
                input = NULL, output = NULL, rule = NULL)
}


#' @title Convert a fis
#' @description
#' Convert a fis object from one type to another (e.g. from singleton to non-singleton)
#'
#' @param fis the fis object to be converted
#' @param option the convert option.'s2n': singleton to non-singleton
#' @param ... For 's2n': fuzzification.method, fuzzification.params, firing.method. 
#' See details below for more information.
#' @details
#' * fuzzification.method, fuzzification.params, firing.method - see \code{\link{addvar}}
#'
#' Usage:
#' 1. convertfis(fis, option, mf.params, fuzzification.method, fuzzification.params)
#' 1. convertfis(fis, option, mf.params, fuzzification.method, fuzzification.params, firing.method)
#' @return Membership grade(s)
#' @examples
#' fis <- tipper()
#' fis.ns.1 <- convertfis(fis, option='s2n', fuzzification.method='gauss', fuzzification.params=1)
#' fis.ns.2 <- convertfis(fis, option='s2n', fuzzification.method='gauss', fuzzification.params=1,
#'                          firing.method='tnorm.min.max')
#' @author Chao Chen
#' @export
#' @md
convertfis <- function(fis, option='s2n', ...) {
    params <- list(...)
    params.len <- length(params)

    if(option == 's2n') {
        if(params.len != 2 && params.len != 3) {
            stop("Incorrect number of parameters for non-singleton fuzzification!")
        }

        fis <- convertfis.s2n(fis, ...)
    } else {
        stop("Incorrect option for fis conversion!")
    }
}

convertfis.s2n <- function(fis, fuzzification.method, fuzzification.params, firing.method = 'tnorm.min.max') {

    input.num <- length(fis$input)
    fis.ns <- fis
    if(input.num > 0) {
        for(i in 1:input.num) {
            fis.ns$input[[i]]$firing.method <- firing.method
            fis.ns$input[[i]]$fuzzification.method <- fuzzification.method
            fis.ns$input[[i]]$fuzzification.params <- fuzzification.params
        }
    }

    fis.ns
}


#' Insert a variable
#'
#' Adds an input or output variable to a fis object.
#' @param fis A fis must be provided.
#' @param varType Should be either 'input' or 'output' which represents the type of variable to be created and added.
#' @param varName A string representing the name of the variable.
#' @param varBounds Also known as the 'range', this should be a vector giving a range for the variable, such as 1:10.
#' @param method fuzzification or defuzzification method.
#' * fuzzification: 'gauss', 'gbell', 'tri', or user-defined.
#' * defuzzification: 'centroid', 'cos', 'coh', 'csum' or user-defined.
#' @param params the required parameters for the corresponding fuzzification or defuzzification method.
#' For example, the required parameters for \code{\link{gbell.fuzzification}} are c(a,b)
#' @param firing.method the chosen method for getting the firing strength (for non-singleton fuzzification).
#' * 'tnorm.min.max' - minimum t-norm with maximum membership grade as the firing strength
#' * 'tnorm.prod.max' - product t-norm with maximum membership grade as the firing strength
#' * 'tnorm.min.defuzz.\[method\]' - the firing strength is based on minimum t-norm, and the chosen defuzzification method (e.g. tnorm.min.defuzz.centroid)
#' * 'tnorm.prod.defuzz.\[method\] - the firing strength is based on product t-norm, and the chosen defuzzification method (e.g. tnorm.prod.defuzz.bisector)
#' * 'similarity.set' - Set-theoretic similarity: the ratio between the intersection and the union of two fuzzy sets
#' @return A fis with the new variable added.
#' @examples
#' fis <- newfis('tipper')
#' fis <- addvar(fis, 'input', 'service', c(0, 10))
#' fis <- addvar(fis, 'input', 'service', c(0, 10), 'gauss', 0.5, 'tnorm.min.max')
#' @export
#' @md
addvar <- function(fis, varType, varName, varBounds, method = NULL, params = NULL, firing.method = 'tnorm.min.max') {
    if (varType == "input") {
        fis$input <- append(fis$input,
                            list(list(name = varName, range = range(varBounds), fuzzification.method = method, fuzzification.params = params, mf = NULL, firing.method = firing.method)))
    }
    else {
        fis$output <- append(fis$output,
                             list(list(name = varName, range = range(varBounds), mf = NULL, defuzzification.method = method)))
    }
    fis
}


#' Insert a membership function.
#'
#' Adds a membership function to a variable of a fis object.
#' @param fis A fis structure is to be provided.
#' @param varType Should be either 'input' or 'output', which relates to the type of variable (stored on the existing fis structure) that the membership function will be added to.
#' @param varIndex Should be an integer value representing the index value of the input or output variable that the membership function will be added to (base 1).
#' @param mfName Membership function name to be declared, for example (Poor,Good)
#' @param mfType Membership function type to be declared, for example (trimf, trapmf)
#' @param mfParams The value of membership function.
#' @return A fis structure with the new membership function added.
#' @examples
#' fis <- newfis('tipper')
#' fis <- addvar(fis, 'input', 'service', c(0, 10))
#' fis <- addmf(fis, 'input', 1, 'poor', 'gaussmf', c(1.5, 0))
#' @export
addmf <- function(fis, varType, varIndex, mfName, mfType, mfParams) {
    if (varType == "input") {
        if (varIndex <= length(fis$input)) {
            fis$input[[varIndex]]$mf <- append(fis$input[[varIndex]]$mf,
                                               list(list(name = mfName, type = mfType, params = mfParams, perturbation = NULL)))
        }
    }
    else {
        if (varIndex <= length(fis$output)) {
            fis$output[[varIndex]]$mf <- append(fis$output[[varIndex]]$mf,
                                                list(list(name = mfName, type = mfType, params = mfParams, perturbation = NULL)))
        }
    }
    fis
}


#' Inserts a rule
#'
#' Adds a rule to a fis object.
#'
#' @param fis  A fis structure is to be provided.
#' @param ruleList  A vector of length m + n + 2, where m is the number of input variables of a fis. \cr Each column in 'm' has a number which refers to the membership function of that input variable. \cr Columns under 'n' refer to an output variable of a fis, where the value refers to the membership function of that output variable. \cr Finally, the '2' remaining columns refer to the weight to be applied to the rule (m + n + 1) and the fuzzy operator for the rule's antecedent (1 = AND, 2 = OR).
#' @details
#' For example, if one has a fis with 2 input variables, and 1 output variable, each of which have 3 membership functions (the amount of membership functions need not be the same). The following rule: 1 3 2 1 2 will mean m = 2 (for 2 input variables), n = 1 (for 1 output variable), and the last 2 columns represent weight and fuzzy operator for the rule's antecedent respectively. \cr\cr
#' The first column refers to the first input variable's membership function at index 1.\cr\cr
#' The second column refers to the second input variable's membership function at index 3.\cr\cr
#' The third column refers to the first output variable's membership function at index 2.\cr\cr
#' The fourth column refers to the weight to be applied to the rule.\cr\cr
#' The fifth column refers to the fuzzy operator for the rule's antecedent (in this case it represents 'OR').
#' @return A fis structure with the new rule added.
#' @examples
#' fis <- tipper()
#' ruleList <- rbind(c(1,1,1,1,2), c(2,0,2,1,1), c(3,2,3,1,2))
#' fis <- addrule(fis, ruleList)
#' @export
addrule <- function(fis, ruleList) {
    fis$rule <- rbind(fis$rule, rbind(ruleList))
    rownames(fis$rule) <- NULL
    fis
}

#' Showing rule from fis object
#'
#' All the rule is showing from fis object
#'
#' @param fis A fis must be provided.
#' @return Show the total of rules inside fis object
#' @examples
#' fis <- tipper()
#' ruleList <- rbind(c(1,1,1,1,2), c(2,0,2,1,1), c(3,2,3,1,2))
#' fis <- addrule(fis, ruleList)
#' showrule(fis)
#' @export
showrule <- function(fis) {
    NumInputs = length(fis$input)
    NumOutputs = length(fis$output)
    NumRules = nrow(fis$rule)
    frow = 0

    if (!is.null(NumRules)) {
        for (i in 1:NumRules) {
            frow = frow + 1
            cat(frow, '. If ', sep = '');
            for (j in 1:NumInputs)
            {
                if (fis$rule[i, j] != 0)
                {
                    cat('(', fis$input[[j]]$name, ' is ', sep = '')
                    if (fis$rule[i, j] < 0) cat('not ', sep = '')
                    cat(fis$input[[j]]$mf[[abs(fis$rule[i, j])]]$name, ') ', sep = '')
                }
                if (j < NumInputs &&
                    fis$rule[i, j] != 0 &&
                    fis$rule[i, j + 1] != 0)
                {
                    if (fis$rule[i, NumInputs + NumOutputs + 2] == 1)
                        cat('and ', sep = '')
                    else
                        cat('or ', sep = '')
                }
            }
            cat('then ', sep = '')
            for (j in 1:NumOutputs)
            {
                if (fis$rule[i, NumInputs + j] != 0)
                {
                    cat('(', fis$output[[j]]$name, ' is ', sep = '')
                    if (fis$rule[i, NumInputs + j] < 0) cat('not ', sep = '')
                    cat(fis$output[[j]]$mf[[abs(fis$rule[i, NumInputs + j])]]$name, ') ', sep = '')
                }
            }
            cat('(', fis$rule[i, NumInputs + NumOutputs + 1], ')\n', sep = '')
        }
    }
}

#' Produce a graphical evaluated fuzzy inference system.
#'
#' Produces a three dimensional graphical view of a specific fis object. This function is only  works for FIS structures with 3 variables. It will only work for 2 inputs, and 1 output.
#'
#' @param fis A fis must be provided.
#' @param ix1 Optional input (1)
#' @param ix2 Optional input (2)
#' @param ox1 Optional output
#' @return  A three dimensional graphical model generated from the fis and other optional parameters.
#' @examples
#' fis <- tipper()
#' gensurf(fis)
#' @export
gensurf <- function(fis, ix1 = 1, ix2 = 2, ox1 = 1) {
    i1 = fis$input[[ix1]]
    i2 = fis$input[[ix2]]
    o1 = fis$output[[ox1]]

    x = seq(i1$range[1], i1$range[2], length = 15)
    y = seq(i2$range[1], i2$range[2], length = 15)
    m = as.matrix(expand.grid(x, y))

    o = evalfis(m, fis)
    z = matrix(o[, ox1], 15, 15, byrow = F)

    h = (z[-15, -15] +
        z[-1, -15] +
        z[-15, -1] +
        z[-1, -1]) / 4
    h = floor((h - min(h, na.rm = T)) / (max(h, na.rm = T) - min(h, na.rm = T)) * 14 + .5) + 1

    persp(x, y, z,
          xlab = i1$name, ylab = i2$name, zlab = o1$name,
          theta = -30, phi = 30, col = rainbow(15)[16 - h], ticktype = 'detailed')
}


#' Plots a 2D graph of all membership functions in a variable.
#'
#' Plots a 2D graph of all membership functions from the specified variable which must be part of a fis object.
#'
#' @param fis Requires an existing fis as an argument.
#' @param varType Can be either 'input' or 'output', representing the type of variable.
#' @param xx primary inputs for extra lines
#' @param timelimit for perturbation
#' @param varIndex A numerical integer, representing the index of the input or output variable whose membership functions shall be plotted (base 1).
#' @param xlab X axis label using font, size and color
#' @param ylab Y axis label, same font attributes as xlab
#' @param main The main title (on top)
#' @return A two dimensional graph displaying all the membership functions of a given variable.
#' @examples
#' fis <- tipper()
#' plotmf(fis, "input", 1)
#' @export
plotmf <- function(fis, varType, varIndex, xx = NULL, timelimit = 0, xlab = NULL, ylab = NULL, main = NULL) {

    if (varType == 'input') {
        var <- fis$input[[varIndex]]
    } else {
        var <- fis$output[[varIndex]]
    }

    if (is.null(xlab))xlab = var$name
    plotvar(var, xx, timelimit, xlab, ylab, main)
}


#' Plot membership functions for an ANFIS object
#'
#' Plots a 2D graph of all membership functions from the specified variable which must be part of an anfis object.
#'
#' @param anfis Requires an existing anfis as an argument.
#' @param varType Can be either 'input' or 'output', representing the type of variable.
#' @param xx primary inputs for extra lines
#' @param timelimit for perturbation
#' @param varIndex A numerical integer, representing the index of the input or output variable whose membership functions shall be plotted (base 1).
#' @param xlab X axis label using font, size and color
#' @param ylab Y axis label, same font attributes as xlab
#' @param main The main title (on top)
#' @return A two dimensional graph displaying all the membership functions of a given variable.
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
#' anfis.plotmf(anfis, 'input', 1)
#' anfis.plotmf(anfis.final, 'input', 1)
#' @export
anfis.plotmf <- function(anfis, varType, varIndex, xx = NULL, timelimit = 0, xlab = NULL, ylab = NULL, main = NULL) {

    if (varType == 'input') {
        var <- list()
        var$name <- anfis$layer[[1]]$node[[varIndex]]$name
        var$range <- anfis$layer[[1]]$node[[varIndex]]$range
        var$mf <- list()

        fan.out <- anfis$layer[[1]]$node[[varIndex]]$fan.out
        for (i in fan.out) {
            name <- anfis$layer[[2]]$node[[i]]$mf.name
            type <- anfis$layer[[2]]$node[[i]]$mf.type
            params <- anfis$layer[[2]]$node[[i]]$mf.params
            var$mf <- append(var$mf, list(list(name = name, type = type, params = params)))
        }
    } else {
        stop("plots for output membership functions not supported")
    }

    if (is.null(xlab))xlab = var$name
    plotvar(var, xx, timelimit, xlab, ylab, main)
}


plotvar <- function(var, xx = NULL, timelimit = 0, xlab = NULL, ylab = NULL, main = NULL) {
    point_n = 501

    x = seq(var$range[1], var$range[2], length = point_n)
    x <- round(x, 4)
    y = c(rep(0, point_n - 1), 1)
    plot_graph(x, y, title = main, xlab = xlab, ylab = ylab)

    if (timelimit == 0) {
        for (i in 1:length(var$mf)) {
            y = evalmf(x, var$mf[[i]]$type, var$mf[[i]]$params)
            if (length(y) >= (2 * length(x))) {
                polygon(c(x, rev(x)), c(y[1:point_n], y[(2 * point_n):(point_n + 1)]), border = T, col = gray((i + 2) / 10))
            } else {
                lines(x, y, col = i, lwd = 3)
            }

            if (length(xx) > 0) {
                for (j in 1:length(xx)) {
                    tmp <- evalmf(xx[j], var$mf[[i]]$type, var$mf[[i]]$params)
                    lines(c(xx[j], xx[j], 0), c(0, tmp[1], tmp[1]), col = "red", lty = "dashed");
                    lines(c(xx[j], xx[j], 0), c(0, tmp[2], tmp[2]), col = "red", lty = "dashed");
                }
            }
        }
        for (i in 1:length(var$mf)) {
            y = evalmf(x, var$mf[[i]]$type, var$mf[[i]]$params)
            if (length(y) == point_n * 2)y <- y[1:point_n]
            tx = which.max(y)

            if (tx >= point_n * 0.95) {
                tx = x[point_n * 0.95]
            }
            else
            {
                tx = x[tx + 20]
            }
            text(tx, 1.03 - i / 20, var$mf[[i]]$name, font = 2, cex = .8)
        }
    } else if (timelimit > 0) {
        posx <- c()
        first <- TRUE

        for (t in 1:timelimit) {
            for (i in 1:length(var$mf)) {
                repeat {
                    params <- var$mf[[i]]$params
                    e <- 0
                    perturbation <- var$mf[[i]]$perturbation
                    if (length(perturbation) > 0) {
                        tmp <- perturb(t, params, perturbation)
                        params <- tmp$params
                        e <- tmp$e
                    }
                    y = evalmf(x, var$mf[[i]]$type, params) + e
                    if (all(!is.nan(y)))break
                }
                y[y > 1] = 1
                y[y < 0] = 0
                if (length(y) >= (2 * length(x))) {
                    polygon(c(x, rev(x)), c(y[1:point_n], y[(2 * point_n):(point_n + 1)]), border = F, col = gray((i + 2) / 10))
                } else {
                    lines(x, y, col = gray((i * 2) / 10))
                }

                if (length(xx) > 0) {
                    for (j in 1:length(xx)) {
                        tmp <- evalmf(xx[j], var$mf[[i]]$type, var$mf[[i]]$params)
                        lines(c(xx[j], xx[j], 0), c(0, tmp[1], tmp[1]), col = "red", lty = "dashed");
                        lines(c(xx[j], xx[j], 0), c(0, tmp[2], tmp[2]), col = "red", lty = "dashed");
                    } #end of for
                } #end of if
                if (first)posx[i] <- x[match(TRUE, y == max(y)) + 10]
            } #end of for
            first = FALSE;
        } #end of for
        for (i in 1:length(var$mf)) {
            text(posx[i], 0.95, var$mf[[i]]$name, col = "red", cex = .8)
        }
    }

    if (length(xx) > 0)
        for (i in 1:length(xx)) {
            mtext(paste(xx[i]), side = 1, col = "red", at = c(xx[i]), cex = .8);
        }
}


# plot_graph - used in plotmf
plot_graph <- function(x, y, xlab = NULL, ylab = NULL, title = NULL) {
    plot(x, y, type = "n", col = 2, bty = "o", xaxs = "i", yaxs = "i", cex.axis = .8, cex.lab = .8, main = title, xaxt = "n", yaxt = "n", ylab = "", xlab = "", las = 1)

    point_n <- length(x)
    axis(1, labels = FALSE, tick = TRUE, pos = 0)
    mtext(x[1], side = 1, line = 1, outer = FALSE, adj = 0, cex = .8)

    if (length(xlab) == 0) {
        mtext("x", side = 1, line = 1, outer = FALSE, padj = 0, cex = .8)
    } else {
        mtext(xlab, side = 1, line = 1, outer = FALSE, padj = 0, cex = .8)
    }

    mtext(x[point_n], side = 1, line = 1, outer = FALSE, adj = 1, cex = .8)
    axis(2, labels = FALSE, tick = TRUE, pos = 0)
    mtext(y[1], side = 2, line = 1, outer = FALSE, at = c(0), cex = .8, las = 1)

    if (length(ylab) == 0) {
        mtext("m", side = 2, line = 1, outer = FALSE, padj = 0, cex = .8, las = 1, font = 5)
    } else {
        if (ylab %in% c('u', 'm')) {
            mtext(ylab, side = 2, line = 1, outer = FALSE, padj = 0, cex = .8, las = 1, font = 5)
        } else {
            mtext(ylab, side = 2, line = 1, outer = FALSE, padj = 0, cex = .8)
        }
    }

    mtext(y[point_n], side = 2, line = 1, outer = FALSE, at = c(1), cex = .8, las = 1)
}


#' Write a fis object to a .fis file.
#'
#' Write a fis object to a file with the .fis extension.
#'
#' @param fis The fuzzy inference system data structure to be saved.
#' @param fileName filename
#' @export
writefis <- function(fis, fileName = 'fuzzy.fis') {
    fileText = NULL

    fileText[1] = "% R-Fuzzy (C) J.M.Garibaldi, 1st Oct 2004 $Revision: 0.1$"
    NumInputs = length(fis$input)
    NumOutputs = length(fis$output)

    NumInputMFs = NULL
    for (i in 1:NumInputs) {
        NumInputMFs[i] = length(fis$input[[i]]$mf)
    }

    NumOutputMFs = NULL
    for (i in 1:NumOutputs) {
        NumOutputMFs[i] = length(fis$output[[i]]$mf)
    }

    NumRules = nrow(fis$rule)

    fileText[2] = "[System]"
    fileText[3] = paste("Name='", fis$name, "'", sep = "")
    fileText[4] = paste("Type='", fis$type, "'", sep = "")
    fileText[5] = paste("NumInputs=", NumInputs, sep = "")
    fileText[6] = paste("NumOutputs=", NumOutputs, sep = "")
    fileText[7] = paste("NumRules=", NumRules, sep = "")
    fileText[8] = paste("AndMethod='", fis$andMethod, "'", sep = "")
    fileText[9] = paste("OrMethod='", fis$orMethod, "'", sep = "")
    fileText[10] = paste("ImpMethod='", fis$impMethod, "'", sep = "")
    fileText[11] = paste("AggMethod='", fis$aggMethod, "'", sep = "")
    fileText[12] = paste("DefuzzMethod='", fis$defuzzMethod, "'", sep = "")
    fileText[13] = paste("mfType='", fis$mfType, "'", sep = "")
    fileText[14] = ""
    line = 15

    for (i in 1:NumInputs) {
        fileText[line] = paste("[Input", i, "]", sep = "")
        line = line + 1
        fileText[line] = paste("Name='", fis$input[[i]]$name, "'", sep = "")
        line = line + 1
        fileText[line] = paste("Range=[", fis$input[[i]]$range[1], " ", fis$input[[i]]$range[2], "]", sep = "")
        line = line + 1
        fileText[line] = paste("fuzzification.method='", fis$input[[i]]$fuzzification.method, "'", sep = "")
        line = line + 1
        fileText[line] = paste("fuzzification.params=[", paste(fis$input[[i]]$fuzzification.params, collapse = " "), "]", sep = "")
        line = line + 1
        fileText[line] = paste("firing.method='", fis$input[[i]]$firing.method, "'", sep = "")
        line = line + 1
        fileText[line] = paste("NumMFs=", NumInputMFs[i], sep = "")
        line = line + 1

        for (j in 1:NumInputMFs[i]) {
            part1 = paste("MF", j, "='", fis$input[[i]]$mf[[j]]$name, "':", sep = "")
            part2 = paste("'", fis$input[[i]]$mf[[j]]$type, "',", sep = "")
            part3 = paste("[", paste(fis$input[[i]]$mf[[j]]$params, collapse = " "), "]", sep = "")
            fileText[line] = paste(part1, part2, part3, sep = "")
            line = line + 1
        }

        fileText[line] = ""
        line = line + 1
    }

    for (i in 1:NumOutputs) {
        fileText[line] = paste("[Output", i, "]", sep = "")
        line = line + 1
        fileText[line] = paste("Name='", fis$output[[i]]$name, "'", sep = "")
        line = line + 1
        fileText[line] = paste("Range=[", fis$output[[i]]$range[1], " ", fis$output[[i]]$range[2], "]", sep = "")
        line = line + 1
        fileText[line] = paste("NumMFs=", NumOutputMFs[i], sep = "")
        line = line + 1

        for (j in 1:NumOutputMFs[i]) {
            part1 = paste("MF", j, "='", fis$output[[i]]$mf[[j]]$name, "':", sep = "")
            part2 = paste("'", fis$output[[i]]$mf[[j]]$type, "',", sep = "")
            part3 = paste("[", paste(fis$output[[i]]$mf[[j]]$params, collapse = " "), "]", sep = "")
            fileText[line] = paste(part1, part2, part3, sep = "")
            line = line + 1
        }

        fileText[line] = ""
        line = line + 1
    }

    fileText[line] = "[Rules]"
    line = line + 1
    for (i in 1:NumRules) {
        part1 = paste(fis$rule[i, 1:NumInputs], collapse = " ")
        part2 = paste(fis$rule[i, (NumInputs + 1):(NumInputs + NumOutputs)], collapse = " ")
        part3 = paste(" (", fis$rule[i, NumInputs + NumOutputs + 1], ") : ", fis$rule[i, NumInputs + NumOutputs + 2], sep = "")
        fileText[line] = paste(part1, ", ", part2, part3, sep = "")
        line = line + 1
    }

    fileText[line] = ""
    writeLines(fileText, fileName)
}


#' Read a fis object from a .fis file.
#'
#' Reads a fis object from a file with the .fis extension, and converts it into a data structure to be used within the environment.
#'
#' @param fileName Should be an absolute path given as a string to the file to be read, with escaped backslashes.
#' @return A fis structure with its values generated from that of the files.
#' @export
readfis <- function(fileName) {
    fileText <- readLines(fileName)
    if (length(fileText) == 0)
        stop('Zero length file!')
    fis <- list()
    line <- 1

    # structure parameters
    line = charmatch('[System]', fileText)
    if (is.na(line) || line == 0)
        stop(paste("No '[System]' line in file", fileName))
    line <- line + 1

    # defaults in case the user has omitted them
    Name <- 'untitled'
    Type <- 'mamdani'
    AndMethod <- 'min'
    OrMethod <- 'max'
    ImpMethod <- 'min'
    AggMethod <- 'max'
    DefuzzMethod <- 'centroid'
    mfType <- 't1'

    # evaluate the values from the file
    while (all(is.na(charmatch(c('[In', '[Out', '[Rules'), fileText[line])))) {
        eval(parse(text = fileText[line]))
        line <- line + 1
    }

    # create a FIS with the given structure
    fis <- list(name = Name, type = Type, mfType = mfType,
                andMethod = AndMethod, orMethod = OrMethod,
                impMethod = ImpMethod, aggMethod = AggMethod,
                defuzzMethod = DefuzzMethod,
                input = NULL, output = NULL, rule = NULL)

    # now begin with the inputs
    for (varIndex in 1:NumInputs) {
        while (is.na(charmatch('[Input', fileText[line])))
            line <- line + 1

        fuzzification.method <- 'singleton.fuzzification'
        fuzzification.params <- NULL
        firing.method <- 'tnorm.min.max'

        # name and range (needs processing)
        eval(parse(text = fileText[line + 1]))
        rangeStr <- fileText[line + 2]
        rangeSplit <- unlist(strsplit(rangeStr, "[][ ]"))
        rangeText <- paste(rangeSplit[1], 'c(', paste(rangeSplit[-1], collapse = ','), ')')
        eval(parse(text = rangeText))
        eval(parse(text = fileText[line + 3]))
        # eval(parse(text = fileText[line + 4]))
        paramText <- fileText[line + 4]
        paramText <- unlist(strsplit(paramText, "[=]"))
        paramText <- unlist(strsplit(paramText[2], "[][ ]"))
        paramText <- paste('fuzzification.params=', 'c(', paste(paramText[-1], collapse = ','), ')')
        eval(parse(text = paramText))
        paramText <- fileText[line + 5]
        eval(parse(text = paramText))

        # now add the variable to the FIS
        fis <- addvar(fis, 'input', Name, Range, fuzzification.method, fuzzification.params, firing.method)

        # number of membership functions
        eval(parse(text = fileText[line + 6]))
        line <- line + 7

        for (MFIndex in 1:NumMFs) {
            # MF1='very low':'trapmf',[13 13 32 35]
            mfStr <- fileText[line]
            #browser()
            mfSplit <- unlist(strsplit(mfStr, "[',]"))
            mfName <- mfSplit[2]
            mfType <- mfSplit[4]
            paramSplit <- unlist(strsplit(mfSplit[6], "[][ ]"))
            paramText <- paste('mfParams=', 'c(', paste(paramSplit[-1], collapse = ','), ')')
            eval(parse(text = paramText))

            # now add the membership function to the FIS
            fis <- addmf(fis, 'input', varIndex, mfName, mfType, mfParams)
            line <- line + 1
        }
    }

    # now for the outputs
    for (varIndex in 1:NumOutputs) {
        while (is.na(charmatch('[Output', fileText[line])))
            line <- line + 1

        # name and range (needs processing)
        eval(parse(text = fileText[line + 1]))
        rangeStr <- fileText[line + 2]
        rangeSplit <- unlist(strsplit(rangeStr, "[][ ]"))
        rangeText <- paste(rangeSplit[1], 'c(', paste(rangeSplit[-1], collapse = ','), ')')
        eval(parse(text = rangeText))

        # now add the variable to the FIS
        fis <- addvar(fis, 'output', Name, Range)

        # number of membership functions
        eval(parse(text = fileText[line + 3]))
        line <- line + 4

        for (MFIndex in 1:NumMFs) {
            mfStr <- fileText[line]
            #browser()
            mfSplit <- unlist(strsplit(mfStr, "[',]"))
            mfName <- mfSplit[2]
            mfType <- mfSplit[4]
            paramSplit <- unlist(strsplit(mfSplit[6], "[][ ]"))
            paramText <- paste('mfParams= c(', paste(paramSplit[-1], collapse = ','), ')')
            eval(parse(text = paramText))

            # now add the membership function to the FIS
            fis <- addmf(fis, 'output', varIndex, mfName, mfType, mfParams)
            line <- line + 1
        }
    }

    # now for the rules
    while (is.na(charmatch('[Rules]', fileText[line])))
        line <- line + 1
    line <- line + 1

    for (ruleIndex in 1:NumRules) {
        ruleStr <- fileText[line]
        ruleSplit <- unlist(strsplit(ruleStr, "[ ,():]"))
        ruleSplit <- ruleSplit[nchar(ruleSplit) > 0]
        ruleText <- paste('ruleList= c(', paste(ruleSplit, collapse = ','), ')')
        eval(parse(text = ruleText))

        # now add the rule to the FIS
        fis = addrule(fis, ruleList)
        line <- line + 1
    }

    fis
}


#' Show a fis object.
#'
#' Shows a fis and all its data in an ordered format on the console.
#'
#' @param fis Requires a fis structure to be displayed.
#' @return Returned the organised text regarding the fis is output to console.
#' @examples
#' fis <- tipper()
#' showfis(fis)
#' @export
showfis <- function(fis) {
    NumInputs = length(fis$input)
    NumOutputs = length(fis$output)

    NumInputMFs = NULL
    for (i in 1:NumInputs) {
        NumInputMFs[i] = length(fis$input[[i]]$mf)
    }

    NumOutputMFs = NULL
    for (i in 1:NumOutputs) {
        NumOutputMFs[i] = length(fis$output[[i]]$mf)
    }

    NumRules = nrow(fis$rule)

    cat('1.  Name             ', fis$name, '\n')
    cat('2a. Type             ', fis$type, '\n')
    cat('2b. mfType           ', fis$mfType, '\n')
    cat('3.  Inputs/Outputs   ', '[', NumInputs, NumOutputs, ']', '\n')
    cat('4.  NumInputMFs      ', '[', NumInputMFs, ']', '\n')
    cat('5.  NumOutputMFs     ', '[', NumOutputMFs, ']', '\n')
    cat('6.  NumRules         ', NumRules, '\n')
    cat('7.  AndMethod        ', fis$andMethod, '\n')
    cat('8.  OrMethod         ', fis$orMethod, '\n')
    cat('9.  ImpMethod        ', fis$impMethod, '\n')
    cat('10. AggMethod        ', fis$aggMethod, '\n')
    cat('11. DefuzzMethod     ', fis$defuzzMethod, '\n')
    frow = 11

    if (NumInputs > 0) {
        for (i in 1:NumInputs) {
            cat(frow + i, '. InLabels          ', fis$input[[i]]$name, '\n', sep = '')
        }
        frow = frow + i
    }

    if (NumOutputs > 0) {
        for (i in 1:NumOutputs) {
            cat(frow + i, '. OutLabels         ', fis$output[[i]]$name, '\n', sep = '')
        }
        frow = frow + i
    }

    if (NumInputs > 0) {
        for (i in 1:NumInputs) {
            cat(frow + i, '. ', sep = '');
            cat('InRange          ', '[', fis$input[[i]]$range, ']', '\n')
        }
        frow = frow + i
    }

    if (NumOutputs > 0) {
        for (i in 1:NumOutputs) {
            cat(frow + i, '. ', sep = '');
            cat('OutRange         ', '[', fis$output[[i]]$range, ']', '\n')
        }
        frow = frow + i
    }

    if (NumInputs > 0) {
        for (i in 1:NumInputs) {
            for (j in 1:NumInputMFs[i]) {
                frow = frow + 1
                cat(frow, '. InMFLabels        ', fis$input[[i]]$mf[[j]]$name, '\n', sep = '')
            }
        }
    }

    if (NumOutputs > 0) {
        for (i in 1:NumOutputs) {
            for (j in 1:NumOutputMFs[i]) {
                frow = frow + 1
                cat(frow, '. OutMFLabels          ', fis$output[[i]]$mf[[j]]$name, '\n', sep = '')
            }
        }
    }

    if (NumInputs > 0) {
        for (i in 1:NumInputs) {
            for (j in 1:NumInputMFs[i]) {
                frow = frow + 1
                cat(frow, '. InMFTypes         ', fis$input[[i]]$mf[[j]]$type, '\n', sep = '')
            }
        }
    }

    if (NumOutputs > 0) {
        for (i in 1:NumOutputs) {
            for (j in 1:NumOutputMFs[i]) {
                frow = frow + 1
                cat(frow, '. OutMFTypes           ', fis$output[[i]]$mf[[j]]$type, '\n', sep = '')
            }
        }
    }

    if (NumInputs > 0) {
        for (i in 1:NumInputs) {
            for (j in 1:NumInputMFs[i]) {
                frow = frow + 1
                cat(frow, '. ', sep = '');
                cat('InMFParams       ', '[', fis$input[[i]]$mf[[j]]$params, ']', '\n')
            }
        }
    }

    if (NumOutputs > 0) {
        for (i in 1:NumOutputs) {
            for (j in 1:NumOutputMFs[i]) {
                frow = frow + 1
                cat(frow, '. ', sep = '');
                cat('OutMFParams      ', '[', fis$output[[i]]$mf[[j]]$params, ']', '\n')
            }
        }
    }

    if (!is.null(NumRules)) {
        for (i in 1:NumRules) {
            frow = frow + 1
            cat(frow, '. ', sep = '');
            cat('Rule Antecedent   [', fis$rule[i, 1:NumInputs], ']', '\n')
        }

        for (i in 1:NumRules) {
            frow = frow + 1
            cat(frow, '. ', sep = '');
            cat('Rule Consequent  ', fis$rule[i, (NumInputs + 1):(NumInputs + NumOutputs)], '\n')
        }

        for (i in 1:NumRules) {
            frow = frow + 1
            cat(frow, '. ', sep = '');
            cat('Rule Weight      ', fis$rule[i, NumInputs + NumOutputs + 1], '\n')
        }

        for (i in 1:NumRules) {
            frow = frow + 1
            cat(frow, '. ', sep = '');
            cat('Rule Connection  ', fis$rule[i, NumInputs + NumOutputs + 2], '\n')
        }
    }
}


#' Evaluate a Fuzzy Inference System (fis)
#'
#' Returns an evaluated crisp value for a given fis structure.
#'
#' @param  input_stack A matrix representing the input stack, number of inputs (columns) by number of outputs (rows).
#' @param  fis A fis must be provided.
#' @param  time default 1
#' @param  point_n number of discretised points, default 101
#' @param  draw whether to draw, TRUE or FALSE
#' @return Returns a matrix of evaluated values.
#' @examples
#' Input_data <- matrix((1:2),1,2)
#' fis <- tipper()
#' evalfis(Input_data, fis)
#' @export
evalfis <- function(input_stack, fis, time = 1, point_n = 101, draw = FALSE) {

    # tsk implementation, added by CC, 20/11/2019
    if (fis$type == 'tsk') {
        if (!exists("GLOBAL_FIS") || !identical(fis, .GlobalEnv$GLOBAL_FIS)) {
            .GlobalEnv$GLOBAL_FIS <- fis
            .GlobalEnv$GLOBAL_TSK <- anfis.builder(fis)
        }

        if (is.null(nrow(input_stack))) {
            input_stack <- matrix(input_stack, nrow = 1)
        }
        out_stack <- anfis.eval(.GlobalEnv$GLOBAL_TSK, input_stack)

        return(out_stack)
    }


    first <- FALSE
    dir = getwd()

    #to make sure there wasn't already a global fis, 
    #or, if there was, it's not identical to the one being passed for initialisation
    if (!exists("GLOBAL_FIS") || !identical(fis, .GlobalEnv$GLOBAL_FIS)) {
        first <- TRUE

        # Add '.GlobalEnv$' to replace super assignment '<<-' (by Chao & Tajul)

        #Get details of fis and create matrices
        #print("initialising ...")
        .GlobalEnv$GLOBAL_FIS <- fis
        .GlobalEnv$FIS_TYPE <- fis$type
        .GlobalEnv$IN_N <- length(fis$input)
        .GlobalEnv$OUT_N <- length(fis$output)
        .GlobalEnv$IN_MF_N <- NULL
        .GlobalEnv$OUT_MF_N <- NULL
        for (i in 1:.GlobalEnv$IN_N)
            .GlobalEnv$IN_MF_N[i] <- length(fis$input[[i]]$mf)
        for (i in 1:.GlobalEnv$OUT_N)
            .GlobalEnv$OUT_MF_N[i] <- length(fis$output[[i]]$mf)

        .GlobalEnv$RULE_N <- nrow(fis$rule)
        .GlobalEnv$RULE_LIST <- fis$rule[, 1:(.GlobalEnv$IN_N + .GlobalEnv$OUT_N)]
        .GlobalEnv$RULE_ANTE <- fis$rule[, 1:.GlobalEnv$IN_N]
        .GlobalEnv$RULE_CONS <- fis$rule[, (.GlobalEnv$IN_N + 1):(.GlobalEnv$IN_N + .GlobalEnv$OUT_N)]
        .GlobalEnv$RULE_WEIGHT <- fis$rule[, .GlobalEnv$IN_N + .GlobalEnv$OUT_N + 1]
        .GlobalEnv$AND_OR <- fis$rule[, .GlobalEnv$IN_N + .GlobalEnv$OUT_N + 2]
        .GlobalEnv$OUT_MF_N_sum <- sum(.GlobalEnv$OUT_MF_N)
        .GlobalEnv$IN_MF_N_sum <- sum(.GlobalEnv$IN_MF_N)
        .GlobalEnv$AND_METHOD <- fis$andMethod
        .GlobalEnv$OR_METHOD <- fis$orMethod
        .GlobalEnv$IMP_METHOD <- fis$impMethod
        .GlobalEnv$AGG_METHOD <- fis$aggMethod
        .GlobalEnv$DEFUZZ_METHOD <- fis$defuzzMethod

        # Get input params and types into globals
        .GlobalEnv$IN_TYPE <- NULL
        .GlobalEnv$IN_PARAMS <- list()
        .GlobalEnv$IN_PERTURBATION <- list()
        .GlobalEnv$OUT_PERTURBATION <- list()

        # Check to see with type of FL is being used
        .GlobalEnv$MF_TYPE <- 1
        .GlobalEnv$OUT_COUPLE <- 0
        .GlobalEnv$COUPLE <- 0
        .GlobalEnv$IN_COUPLE <- 0
        if (fis$mfType == 'it2') {
            .GlobalEnv$MF_TYPE <- 2
            .GlobalEnv$OUT_COUPLE <- c(0, .GlobalEnv$OUT_MF_N_sum)
            .GlobalEnv$COUPLE <- c(0, .GlobalEnv$RULE_N)
            .GlobalEnv$IN_COUPLE <- c(0, .GlobalEnv$IN_MF_N_sum)
        }

        #Get the type, parameters and perturbation functions for each INPUT mf
        idx = 1
        for (i in 1:.GlobalEnv$IN_N) {
            for (j in 1:.GlobalEnv$IN_MF_N[i]) {
                .GlobalEnv$IN_TYPE[idx] <- fis$input[[i]]$mf[[j]]$type
                .GlobalEnv$IN_PARAMS[idx] <- list(fis$input[[i]]$mf[[j]]$params)
                .GlobalEnv$IN_PERTURBATION[idx] <- list(fis$input[[i]]$mf[[j]]$perturbation)
                idx = idx + 1
            }
        }

        #Get perturbation function for OUTPUT mfs
        .GlobalEnv$no_output_perturbation <- TRUE
        idx = 1
        for (i in 1:.GlobalEnv$OUT_N) {
            for (j in 1:.GlobalEnv$OUT_MF_N[i]) {
                perturbation <- fis$output[[i]]$mf[[j]]$perturbation

                #if this mf has a perturbation function set flag to FALSE
                if (length(perturbation) > 0) {
                    .GlobalEnv$no_output_perturbation <- FALSE
                }
                .GlobalEnv$OUT_PERTURBATION[idx] <- list(perturbation)
                idx = idx + 1
            }
        }

        #Get range of each OUTPUT mf, and create a sequence over that range
        .GlobalEnv$OUT_RANGE <- matrix(0, .GlobalEnv$OUT_N, 2)
        .GlobalEnv$rangex <- matrix(0, point_n, .GlobalEnv$OUT_N)
        for (i in 1:.GlobalEnv$OUT_N) {
            .GlobalEnv$OUT_RANGE[i,] <- fis$output[[i]]$range
            .GlobalEnv$rangex[, i] <- seq(.GlobalEnv$OUT_RANGE[i, 1], .GlobalEnv$OUT_RANGE[i, 2], length = point_n)
        }

        # allocate matrices for rule consequents and aggregated consequents
        .GlobalEnv$OUT_RULE_CONS <- matrix(0, .GlobalEnv$RULE_N * .GlobalEnv$MF_TYPE, point_n * .GlobalEnv$OUT_N)
        .GlobalEnv$OUT_RULE_AGG <- matrix(0, .GlobalEnv$MF_TYPE, point_n * .GlobalEnv$OUT_N)
    }

    # why first is set to be true here?
    first <- TRUE
    if (first ||
        !exists("no_output_perturbation") ||
        !.GlobalEnv$no_output_perturbation) {

        .GlobalEnv$OUT_TEMP_MF <- matrix(0, (.GlobalEnv$OUT_MF_N_sum) * .GlobalEnv$MF_TYPE + .GlobalEnv$MF_TYPE, point_n)


        # -- BEGIN, 15 Mar 2020, Chao -- 
        # -- When AGG_METHOD is 'max', the following line will produce unexpected results. Hence, it is commented out.
        # .GlobalEnv$OUT_TEMP_MF[OUT_COUPLE + 1,] <- 1
        # -- END, 15 Mar 2020, Chao -- 

        idx = 1

        for (i in 1:.GlobalEnv$OUT_N) {

            #if draw is set plot each MF in this output variable and produce a PDF
            filename <- paste(dir, substring(fis$name, length(fis$name)), "output", i, ".pdf", sep = "")
            if (draw && time == 1 && !file.exists(filename)) {
                pdf(filename)
                pushViewport(plotViewport())
                plotmf(fis, "output", i, xlab = tolower(fis$output[[i]]$name), main = fis$output[[i]]$name)
                popViewport();
                dev.off()
            }

            for (j in 1:.GlobalEnv$OUT_MF_N[i]) {

                #Get the parameters, type and perturbation function for each mf in this output var
                params <- fis$output[[i]]$mf[[j]]$params
                type <- fis$output[[i]]$mf[[j]]$type
                perturbation <- .GlobalEnv$OUT_PERTURBATION[[idx]]

                #If there are any perturbation functions, perturb the membership function
                if (length(perturbation) > 0) {
                    count <- 0
                    repeat {
                        tmp <- perturb(time, params, perturbation)
                        params <- tmp$params
                        e <- tmp$e
                        tmp <- matrix(evalmf(.GlobalEnv$rangex[, i],
                                             type, params) + e, .GlobalEnv$MF_TYPE, point_n, byrow = TRUE)
                        if (count > 1000 || all(!is.nan(tmp)))break
                        count <- count + 1
                    }
                } else {
                    tmp <- matrix(evalmf(.GlobalEnv$rangex[, i], type, params), .GlobalEnv$MF_TYPE, point_n, byrow = TRUE)
                }
                #clipping
                tmp[tmp < 0] = 0
                tmp[tmp > 1] = 1

                #Copy perturbed mf
                .GlobalEnv$OUT_TEMP_MF[OUT_COUPLE + idx + 1,] <- tmp
                idx = idx + 1
            }
        }

        # restructure to fit OUT_MF
        idx = matrix(1, .GlobalEnv$RULE_N, 1) %*% cumsum(c(1, .GlobalEnv$OUT_MF_N[0:(.GlobalEnv$OUT_N - 1)])) + abs(.GlobalEnv$RULE_CONS)
        idx[RULE_CONS == 0] = 1

        .GlobalEnv$OUT_MF <- .GlobalEnv$OUT_TEMP_MF[t(idx),]

        # -- BEGIN, 15 April 2020, Chao -- 
        # -- This is part of the solution to address the issue when the FIS only has one rule.
        if (.GlobalEnv$RULE_N == 1) {
            .GlobalEnv$OUT_MF <- t(.GlobalEnv$OUT_MF)
        }
        # -- END, 15 April 2020, Chao -- 

        #Now OUT_MF is RULE_N*OUT_N X point_n matrix
        .GlobalEnv$OUT_MF[t(.GlobalEnv$RULE_CONS) < 0,] <- 1 - .GlobalEnv$OUT_MF[t(.GlobalEnv$RULE_CONS < 0),]

        #if it's a type-2 fis,
        if (.GlobalEnv$MF_TYPE == 2) {

            idx <- idx + .GlobalEnv$OUT_MF_N_sum
            tmp <- .GlobalEnv$OUT_TEMP_MF[t(idx),]

            #Now OUT_MF is RULE_N*OUT_N X point_n matrix
            if (any(.GlobalEnv$RULE_CONS < 0)) {
                tt <- t(.GlobalEnv$RULE_CONS) < 0
                tmp[tt,] <- 1 - tmp[tt,]
                temp <- tmp[tt,]
                tmp[tt,] <- .GlobalEnv$OUT_MF[tt,]
                .GlobalEnv$OUT_MF[tt,] <- temp
            }
            .GlobalEnv$OUT_MF <- rbind(.GlobalEnv$OUT_MF, tmp)
        }
        .GlobalEnv$OUT_MF <- matrix(t(.GlobalEnv$OUT_MF), .GlobalEnv$RULE_N * .GlobalEnv$MF_TYPE, point_n * .GlobalEnv$OUT_N, byrow = TRUE)
    }
    # end of initialisation

    # check for errors in the input stack
    if (is.vector(input_stack)) {
        #Put the inputs into seperate columns
        input_stack = matrix(input_stack, ncol = .GlobalEnv$IN_N, byrow = TRUE) #Error here
    }
    data_n = nrow(input_stack)

    # create output stack
    # out_stack = matrix(0, data_n * .GlobalEnv$MF_TYPE, .GlobalEnv$OUT_N)
    # modified by CC, 19/11/2019
    out_stack = matrix(0, data_n, .GlobalEnv$OUT_N)

    # loop through each input
    for (k in 1:data_n) {
        input = input_stack[k,]
        # create in_temp_mf_value
        in_temp_mf_value = matrix(0, .GlobalEnv$IN_MF_N_sum * .GlobalEnv$MF_TYPE, 1)

        idx = 1
        pr <- c()
        for (i in 1:.GlobalEnv$IN_N) {

            #create a filename
            filename <- paste(dir, substring(fis$name, length(fis$name)), "input", i, ".pdf", sep = "")

            #if draw, create a pdf and plot this input mf
            if (draw && time == 1 && !file.exists(filename)) {
                pdf(filename)
                pushViewport(plotViewport())
                plotmf(fis, "input", i, xlab = tolower(fis$input[[i]]$name), main = fis$input[[i]]$name)
                popViewport()
                dev.off()
            }

            #Create filename
            filename <- paste(dir, substring(fis$name, length(fis$name)), "-fuzzify-input-", i, "-x=", input[i], ".pdf", sep = "")

            #if draw, create pdf and plot input mf, this time alter xx 
            if (draw && time == 1 && !file.exists(filename)) {
                pdf(filename)
                pushViewport(plotViewport())
                plotmf(fis, "input", i, xlab = tolower(fis$input[[i]]$name), main = fis$input[[i]]$name, xx = input[i])
                popViewport()
                dev.off()
            }

            for (j in 1:.GlobalEnv$IN_MF_N[i]) {

                #Get details of this mf
                perturbation <- .GlobalEnv$IN_PERTURBATION[[idx]]
                params <- .GlobalEnv$IN_PARAMS[[idx]]
                type <- .GlobalEnv$IN_TYPE[idx]

                #If there are any perturbation functions
                if (length(perturbation) > 0) {
                    count <- 0
                    #Apply perturbation function
                    repeat {
                        tmp <- perturb(time, params, perturbation)
                        params <- tmp$params
                        e <- tmp$e
                        tmp <- matrix(evalmf(input[i], type, params) + e, .GlobalEnv$MF_TYPE, 1, byrow = TRUE)
                        if (count > 1000 || all(!is.nan(tmp)))break
                        count <- count + 1
                    }
                } else {
                    fuzzification.method <- fis$input[[i]]$fuzzification.method

                    if (is.null(fuzzification.method) || fuzzification.method == 'singleton.fuzzification') {
                        tmp <- matrix(evalmf(input[i], type, params), .GlobalEnv$MF_TYPE, 1, byrow = TRUE)
                    } else {
                        tmp <- matrix(evalmf(input[i], type, params, fuzzification.method, fis$input[[i]]$fuzzification.params, fis$input[[i]]$firing.method, fis$input[[i]]$range), .GlobalEnv$MF_TYPE, 1, byrow = TRUE)
                    }
                }

                #Tidy up
                tmp[tmp > 1] = 1
                tmp[tmp < 0] = 0

                in_temp_mf_value[IN_COUPLE + idx,] <- tmp
                idx = idx + 1
            }
        }

        # restructure row of in_temp_mf_value to fit in_mf_value
        ind = matrix(1, .GlobalEnv$RULE_N, 1) %*% cumsum(c(0, .GlobalEnv$IN_MF_N[0:(.GlobalEnv$IN_N - 1)])) + abs(.GlobalEnv$RULE_ANTE)
        # -- BEGIN 13 Mar 2020, Chao -- 
        # -- Bug fix (when first input is not used)
        ind[ind == 0] <- 1
        # -- END 13 Mar 2020, Chao -- 
        in_mf_value = matrix(in_temp_mf_value[ind,], .GlobalEnv$RULE_N, .GlobalEnv$IN_N)

        # replace dont-care MFs in AND rules with 1, and OR rules with 0
        in_mf_value[which(((.GlobalEnv$AND_OR == 1) * (.GlobalEnv$RULE_ANTE == 0)) == 1)] = 1
        in_mf_value[which(((.GlobalEnv$AND_OR == 2) * (.GlobalEnv$RULE_ANTE == 0)) == 1)] = 0

        # sort out NOTs (negative rule indices)
        idx = which(.GlobalEnv$RULE_ANTE < 0)
        in_mf_value[idx] = 1 - in_mf_value[idx]

        #If it's a type-2 system 
        if (.GlobalEnv$MF_TYPE == 2) {
            #set the ind to the total number of input mfs
            ind <- ind + .GlobalEnv$IN_MF_N_sum
            tmp = matrix(in_temp_mf_value[ind,], .GlobalEnv$RULE_N, .GlobalEnv$IN_N)

            # replace dont-care MFs in AND rules with 1, and OR rules with 0
            tmp[which(((.GlobalEnv$AND_OR == 1) * (.GlobalEnv$RULE_ANTE == 0)) == 1)] = 1
            tmp[which(((.GlobalEnv$AND_OR == 2) * (.GlobalEnv$RULE_ANTE == 0)) == 1)] = 0

            # take care of NOTs (negative rule indices)
            if (length(idx) > 0) {
                temp <- in_mf_value[idx]
                in_mf_value[idx] = 1 - tmp[idx]
                tmp[idx] = temp
            }
            in_mf_value = rbind(in_mf_value, tmp)
        }

        # Get firing strengths for the rules
        f_str <- matrix(0, .GlobalEnv$RULE_N * .GlobalEnv$MF_TYPE, 1)
        if (.GlobalEnv$IN_N == 1) {
            f_str = in_mf_value
        } else {
            and_ind = which(.GlobalEnv$AND_OR == 1)
            or_ind = which(.GlobalEnv$AND_OR == 2)
            f_str[and_ind,] =
                apply(rbind(in_mf_value[and_ind,]), 1, fis$andMethod)

            #Check type
            if (.GlobalEnv$MF_TYPE == 2) {
                f_str[and_ind + .GlobalEnv$RULE_N,] =
                    apply(rbind(in_mf_value[and_ind + .GlobalEnv$RULE_N,]), 1, fis$andMethod)
            }

            f_str[or_ind,] =
                apply(rbind(in_mf_value[or_ind,]), 1, fis$orMethod)

            #Check type
            if (.GlobalEnv$MF_TYPE == 2) {
                f_str[or_ind + .GlobalEnv$RULE_N,] =
                    apply(rbind(in_mf_value[or_ind + .GlobalEnv$RULE_N,]), 1, fis$orMethod)
            }
        }

        #Calculate weighted firing strength
        f_str = f_str * .GlobalEnv$RULE_WEIGHT

        if (all(f_str == 0)) { # if no rules fired,
            for (i in 1:.GlobalEnv$OUT_N) {
                for (j in 1:.GlobalEnv$MF_TYPE) {
                    out_stack[(k - 1) * .GlobalEnv$MF_TYPE + j, i] = NaN
                }
            }
        } else {
            #Check system type and implication methods
            if (.GlobalEnv$FIS_TYPE == 'mamdani') {

                tmp = matrix(f_str, nrow(f_str), point_n * .GlobalEnv$OUT_N)
                if (fis$impMethod == 'prod') {
                    .GlobalEnv$OUT_RULE_CONS <- tmp * .GlobalEnv$OUT_MF
                } else if (fis$impMethod == 'min') {
                    .GlobalEnv$OUT_RULE_CONS <- pmin(tmp, .GlobalEnv$OUT_MF)
                } else {
                    cat('user-defined implication not implemented yet\n')
                }

                #IF draw and it's type-1, plot antecedents, consequents and aggregation and put them in PDFs
                if (draw && .GlobalEnv$MF_TYPE == 1 && time == 1) {
                    input_set <- paste(tolower(fis$input[[1]]$name), "=", input[1])
                    for (i in 2:.GlobalEnv$IN_N)
                        input_set <- paste(input_set, "and", tolower(fis$input[[i]]$name), "=", input[i])

                    for (i in 1:.GlobalEnv$RULE_N) {
                        filename <- paste(dir, substring(fis$name, length(fis$name)), "-ante-eval-", i, "-(", sep = "")
                        s <- paste(filename, input_set, ").pdf", sep = "")

                        if (time == 1 && !file.exists(s)) {
                            pdf(s)
                            pushViewport(plotViewport())
                            drawAnteEval(fis, i, input, in_mf_value[i,], point_n = point_n, title = NULL, label = TRUE)
                            popViewport()
                            dev.off()
                        }

                        filename <- paste(dir, substring(fis$name, length(fis$name)), "-cons-eval-", i, "-(", sep = "")
                        s <- paste(filename, input_set, ").pdf", sep = "")
                        if (time == 1 && !file.exists(s)) {
                            pdf(s)
                            pushViewport(plotViewport())
                            drawConsEval(fis, i, support = f_str[i,], OUT_RULE_CONS = .GlobalEnv$OUT_RULE_CONS[i,], point_n = point_n, title = NULL, label = TRUE)
                            popViewport()
                            dev.off()
                        }

                        filename <- paste(dir, substring(fis$name, length(fis$name)), "-rule-eval-", i, "-(", sep = "")
                        s <- paste(filename, input_set, ").pdf", sep = "")
                        if (time == 1 && !file.exists(s)) {
                            pdf(s)
                            pushViewport(plotViewport())
                            drawRuleEval(fis, i, input = input, in_mf_value = in_mf_value[i,], support = f_str[i,], OUT_RULE_CONS = .GlobalEnv$OUT_RULE_CONS[i,], point_n = point_n, title = NULL, label = TRUE)
                            popViewport()
                            dev.off()
                        }
                    }

                    filename <- paste(dir, substring(fis$name, length(fis$name)), "-rule-evaluation-(", sep = "")
                    s <- paste(filename, input_set, ").pdf", sep = "")
                    if (time == 1 && !file.exists(s)) {
                        pdf(s)

                        pushViewport(viewport(width = 0.9, height = 0.9))

                        drawRuleAggr(fis, input = input, in_mf_value = in_mf_value, f_str = f_str, OUT_RULE_CONS = .GlobalEnv$OUT_RULE_CONS, point_n = point_n, title = NULL)
                        popViewport()
                        dev.off()
                    }
                }


                # Check type and aggregate rule consequents
                for (ii in 1:.GlobalEnv$MF_TYPE)
                    .GlobalEnv$OUT_RULE_AGG[ii,] <- apply(.GlobalEnv$OUT_RULE_CONS[(ii - 1) * .GlobalEnv$RULE_N + (1:.GlobalEnv$RULE_N),], 2, fis$aggMethod)
                #perform type reduction if it's type-2
                if (.GlobalEnv$MF_TYPE == 2) {
                    if (fis$defuzzMethod == "centroid") {
                        rs <- trCentroid(rbind(.GlobalEnv$rangex, .GlobalEnv$rangex), .GlobalEnv$OUT_RULE_AGG)

                        #plot output mf
                        if (draw) plot_it2(.GlobalEnv$rangex, .GlobalEnv$OUT_RULE_AGG, "Y")
                    } else if (fis$defuzzMethod == "cos") { #Centre of Spread, for NSFS
                        if (first ||
                            !exists("no_output_perturbation") ||
                            !.GlobalEnv$no_output_perturbation) {
                            COS <- matrix(0, .GlobalEnv$RULE_N * 2, .GlobalEnv$OUT_N)
                            COS <- trCentroid(rbind(.GlobalEnv$rangex, .GlobalEnv$rangex), .GlobalEnv$OUT_MF)
                        }
                        rs <- trCOS(COS, f_str)
                    } else if (fis$defuzzMethod == "coh") { #Centre of Height?
                        height <- matrix(0, .GlobalEnv$RULE_N, .GlobalEnv$OUT_N)
                        B <- matrix(0, .GlobalEnv$RULE_N * 2, .GlobalEnv$OUT_N)
                        sp <- matrix(0, .GlobalEnv$RULE_N, .GlobalEnv$OUT_N)

                        for (i in 1:.GlobalEnv$RULE_N) {
                            for (j in 1:.GlobalEnv$OUT_N) {
                                maxp <- max(.GlobalEnv$OUT_RULE_CONS[i, (1:point_n) + (j - 1) * point_n])
                                height[i, j] <- mean(.GlobalEnv$rangex[which(.GlobalEnv$OUT_RULE_CONS[i, (1:point_n) + (j - 1) * point_n] == maxp), j])
                                tmp <- abs(.GlobalEnv$rangex[, j] - height[i, j])
                                pos = which(tmp == min(tmp))[1]
                                B[i, j] <- .GlobalEnv$OUT_RULE_CONS[i, pos]
                                B[i + .GlobalEnv$RULE_N, j] <- .GlobalEnv$OUT_RULE_CONS[i + .GlobalEnv$RULE_N, pos]
                                sp[i, j] <- (.GlobalEnv$OUT_RULE_CONS[i, pos] + .GlobalEnv$OUT_RULE_CONS[i + .GlobalEnv$RULE_N, pos]) / 2
                            }
                        }

                        rs <- iterative_method(rbind(height, height), B)
                    } else if (fis$defuzzMethod == "csum") { #cumulative sum?
                        .GlobalEnv$SUMS <- matrix(0, point_n * 2, .GlobalEnv$OUT_N)
                        for (i in 1:.GlobalEnv$OUT_N) {
                            for (j in 1:point_n) {
                                .GlobalEnv$SUMS[j, i] <- sum(.GlobalEnv$OUT_RULE_CONS[1:.GlobalEnv$RULE_N, (i - 1) * point_n + j])
                                .GlobalEnv$SUMS[j + point_n, i] <- sum(.GlobalEnv$OUT_RULE_CONS[(1 + .GlobalEnv$RULE_N):(.GlobalEnv$RULE_N + .GlobalEnv$RULE_N), (i - 1) * point_n + j])
                            }
                        }
                        rs <- iterative_method(rbind(.GlobalEnv$rangex, .GlobalEnv$rangex), .GlobalEnv$SUMS)
                    } else if (fis$defuzzMethod == "user") {   #User defined
                        rs <- iterative_method(matrix(.GlobalEnv$OUT_MF[, 1], 50, 1), f_str)
                    }

                    # defuzzify each output
                    for (i in 1:.GlobalEnv$OUT_N) {
                        # modified by CC, 19/11/2019
                        # out_stack[2 * k - 1, i] = mean(rs[, i])
                        # out_stack[2 * k, i] = rs[1, i] - mean(rs[, i])
                        out_stack[k, i] = mean(rs[, i])
                    }
                } else {

                    # defuzzify each output
                    for (i in 1:.GlobalEnv$OUT_N) {

                        out_stack[k, i] =
                            defuzz(.GlobalEnv$rangex[, i],
                                   .GlobalEnv$OUT_RULE_AGG[1, ((i - 1) * point_n + 1):(i * point_n)], fis$defuzzMethod)
                        # add by Tajul, and modified by Chao
                        .GlobalEnv$D_x <- .GlobalEnv$rangex[, i]
                        .GlobalEnv$D_y <- .GlobalEnv$OUT_RULE_AGG[1, ((i - 1) * point_n + 1):(i * point_n)]
                        # end
                    }

                    if (draw) {
                        filename <- paste(dir, substring(fis$name, length(fis$name)), "-all-evaluations-(", sep = "")
                        s <- paste(filename, input_set, ").pdf", sep = "")
                    }
                    if (draw && time == 1 && !file.exists(filename)) {
                        #pdf(s)
                        pushViewport(viewport(width = 0.9, height = 0.9))
                        drawAllSteps(fis, input = input, in_mf_value = in_mf_value, f_str = f_str, OUT_RULE_CONS = .GlobalEnv$OUT_RULE_CONS,
                                     OUT_RULE_AGG = .GlobalEnv$OUT_RULE_AGG, output = out_stack[k,], point_n = point_n, title = NULL)
                        popViewport()
                        #dev.off()
                    }
                }
            }
            else if (.GlobalEnv$FIS_TYPE == 'sugeno')
            {
                cat('sugeno inference not implemented yet\n')
            }
            else {
                cat('unknown inference type\n')
            }
        } #end of if
    }

    # add by Tajul
    .GlobalEnv$D_out <- round(out_stack, 2)
    #end

    out_stack
}


#' Defuzzify a set of values.
#'
#' Defuzzifies a given set of values using a specified range and defuzzification type producing a crisp value.
#'
#' @param x The range to be applied in the function (numeric vector).
#' @param mf The values to be applied in the function (numeric vector).
#' @param type The defuzzification method type, which should be either 'centroid', 'bisector', 'mom', 'som' or 'lom'.
#' @return Returns a defuzzified crisp value (double).
#' @examples
#' Crisp_value = defuzz(1:10, c(1.5, 5), "centroid")
#' @export
defuzz <- function(x, mf, type) {
    if (type == "centroid")
        if (sum(mf, na.rm = T) != 0) {
            sum(mf * x, na.rm = T) / sum(mf, na.rm = T)
        } else {
            mean(x, na.rm = T)
        }
    else if (type == "bisector") {
        cs = cumsum(mf)
        a2 = sum(mf) / 2
        xs = match(TRUE, cs > a2)
        x[xs - 1] +
            (a2 - cs[xs - 1]) / mf[xs] +
            (x[xs] - x[xs - 1]) / 2
    } else if (type == "mom") {
        mean(x[which(mf == max(mf))])
    } else if (type == "som") {
        x[min(which(mf == max(mf)))]
    } else if (type == "lom") {
        x[max(which(mf == max(mf)))]
    } else {
        NA
    }
}


#############################################
# IT2 defuzzification - need to go through these

#CENTROID
#B: 2*M X point_n * OUT_N matrix
#y: 2*point_n X  OUT_N matrix
#outs : 2*M x OUT_N
trCentroid <- function(y, B) {

    # M <- nrow(B) %/% 2
    # point_n <- nrow(y) %/% 2
    OUT_N <- ncol(y)
    # out <- matrix(0, nrow(B), OUT_N)

    ff <- rotate2(B, OUT_N)
    out <- iterative_method(y, ff)

    out
}


#COS : matrix RULE_N*2 x OUT_N
#f : RULE_N*2 x 1(OUT_N)
#out: OUT_N*2 x 1
trCOS <- function(COS, f) {
    tmp <- center_spread(COS)
    cc <- tmp$c
    s <- tmp$s

    tmp <- center_spread(f)
    h <- tmp$c
    OUT_N <- ncol(COS)
    delta <- tmp$s

    outs <- matrix(0, 2, OUT_N)
    for (i in 1:OUT_N) {
        tmp <- interval_wtdavg(cc[, i], s[, i], h, delta)
        outs[1, i] <- tmp$r_out
        outs[2, i] <- tmp$l_out
    }

    outs
}


#y : M*2  x OUT_N
#f : M*2*N  x OUT_N
#out : 2*N  x OUT_N
#Used in centroid
iterative_method <- function(y, f) {

    OUT_N <- ncol(y)
    N <- nrow(f) %/% nrow(y)
    M <- nrow(y) %/% 2

    tmp <- center_spread(f)
    h <- tmp$c #N*M x OUT_N
    delta <- tmp$s  #N*M x OUT_N

    tmp <- center_spread(y)
    c <- tmp$c
    s <- tmp$s #c,s : M x OUT_N

    out <- matrix(0, 2 * N, OUT_N)
    for (i in 1:OUT_N) {
        for (j in 0:(N - 1)) {
            tmp <- interval_wtdavg(c[, i], s[, i], h[(1:M) + j * M, i], delta[(1:M) + j * M, i])
            out[j + 1, i] <- tmp$r_out
            out[j + 1 + N, i] <- tmp$l_out
        }
    }
    out
}


#x : 2*N x M matrix
#$c,$s : N*M matrix
center_spread <- function(x) {
    N = nrow(x) %/% 2
    M = ncol(x)
    cc <- matrix(0, N, M)
    s <- matrix(0, N, M)
    for (i in 1:N) {
        cc[i,] = (x[i,] + x[i + N,]) / 2
        s[i,] = (x[i,] - cc[i,])
    }

    list(c = cc, s = s)
}


#c,s,h,delta are all vetors length M
interval_wtdavg <- function(c, s, h, delta) {
    lower = h - delta;
    upper = h + delta;

    l_out = adapt(c - s, lower, upper, -1)$outextreme;
    r_out = adapt(c + s, lower, upper, 1)$outextreme;
    list(l_out = l_out, r_out = r_out)
}


#OUT_MF : RULE_N*2 X point_n*OUT_N
#out : RULE_N*2*point_n X OUT_N
#used in Centroid
rotate2 <- function(OUT_MF, OUT_N) {
    point_n <- ncol(OUT_MF) %/% OUT_N
    out <- matrix(0, nrow(OUT_MF) * point_n, OUT_N)
    for (i in 1:OUT_N) {
        tmp <- OUT_MF[, (1:point_n) + (i - 1) * point_n]
        #tmp : RULE_N*2 X point_n
        out[, i] <- t(tmp)
    }
    out
}


# Inputs : "ypoint", "lower" and "upper" are all M-dimensional vectors.
# "ypoint" contains the "y_l"s. "lower" and "upper"
# contain, respectively, the "w_lower" and "w_upper" values for each
# weight "w_l". If "maxflag > 0" (scalar), "S" is maximized, else
# it is minimized.
# Used in interval_wtdavg
adapt <- function(ypoint, lower, upper, maxflag) {
    tmp <- sort(ypoint, index = TRUE)
    z <- tmp$x
    ix <- tmp$ix
    lower_sort <- lower[ix]
    upper_sort <- upper[ix]
    lz <- length(z)

    hl = (lower_sort + upper_sort) / 2;
    S = sum(z * hl) / sum(hl);   # starting point

    eps = 1e-5;   # small quantity to avoid floating point equality problems

    count = 0;
    theta = hl;
    S_new = S + 10 * eps;

    if ((abs(S - z[1]) < eps) | (abs(S - z[lz]) < eps)) {
        outextreme = S;
    } else {
        while (abs(S - S_new) > eps) {
            count = count + 1;

            if (count > 1)
                S = S_new;

            in1 = which(z > (S - eps));
            min1 = min(in1);

            if (min1 > 2) {
                in2 = 1:(min1 - 1);
            } else {
                in2 = 1;
            }

            if (maxflag > 0) {
                theta[in1] = upper_sort[in1];
                theta[in2] = lower_sort[in2];
            } else {
                theta[in1] = lower_sort[in1];
                theta[in2] = upper_sort[in2];


                # To avoid division by zero if all lower_sort=0
                if (abs(S - z[min1]) < eps)
                    theta[min1] = upper_sort[min1];
            } #end if maxflag


            S_new = sum(z * theta) / sum(theta);
        } #end while
        outextreme = S_new;
    } # end if

    list(outextreme = outextreme, count = count, theta = theta)
}


#' @title TSK FIS builder
#' @description
#' To build a one-output TSK FIS by automatically generating the input membership functions and the fuzzy rules
#' @param x.range a vector/matrix as the range of input(s)
#' @param input.num the number of inputs
#' @param input.mf.num a list of the number of membership functions for all inputs
#' @param input.mf.type designed for different membershp function types, however, currently, 'T1' for gbellmf, else 'it2gbellmf'
#' @param rule.num the number of rules
#' @param rule.which selected rules to be used in the full rule list, for example, c(1,2,3) specify the first three rules
#' @param defuzzMethod "default"
#' @param params.ante parameter settings for initialising antecedent membership functions
#' @param params.conse parameter settings for initialising consequent membership functions
#' @author Chao Chen
#' @export

fis.builder <- function(x.range, input.num, input.mf.num, input.mf.type, rule.num = prod(input.mf.num), rule.which = NULL, defuzzMethod = "default", params.ante, params.conse) {
    fis <- newfis('newFIS', fisType = "tsk", andMethod = "prod", orMethod = "max", impMethod = "min", aggMethod = "max", defuzzMethod = defuzzMethod)

    if (!is.null(rule.which)) {
        rule.which <- sort(unique(rule.which))
        if (rule.num != length(rule.which) || !all(rule.which %in% 1:prod(input.mf.num)))
        {
            stop("rule.which is not right!")
        }
    }

    k <- 1
    x.range <- matrix(x.range, ncol = 2)
    for (i in 1:input.num) {
        if (nrow(x.range) == 1) {
            input.range <- c(x.range)
        } else {
            input.range <- x.range[i,]
        }

        fis <- addvar(fis, 'input', paste0('X', i), input.range, 'singleton.fuzzification')
        #fis <- addvar(fis, 'input', paste0('X', i), input.range, 'it2gbell.fuzzification', init.params.it2gbell(input.range,0)[1,])

        if (missing(params.ante)) {
            params.gbell <- init.params.gbell(input.range, input.mf.num[i])
            params.it2gbell <- init.params.it2gbell(input.range, input.mf.num[i])
        } else {
            params.gbell <- as.matrix(params.ante[[i]])
            params.it2gbell <- params.gbell
        }

        for (j in 1:input.mf.num[i]) {
            #if(input.mf.type[i, j] == 1) {
            if (input.mf.type == "T1") {
                fis <- addmf(fis, 'input', i, paste0('A', k), 'gbellmf', params.gbell[j,])
            } else {
                fis <- addmf(fis, 'input', i, paste0('A', k), 'it2gbellmf', params.it2gbell[j,])
            }
            k <- k + 1
        }
    }

    fis <- addvar(fis, 'output', 'Y', NULL, 'defuzzification.method')

    if (missing(params.conse)) {
        params.linearmf <- matrix(0, nrow = rule.num, ncol = input.num + 1)
    } else {
        params.linearmf <- params.conse
    }

    for (i in 1:rule.num) {
        fis <- addmf(fis, 'output', 1, paste0('C', i), 'linearmf', params.linearmf[i,])
    }

    rule.ante <- NULL
    rule.ante <- expand.grid(lapply(input.mf.num, seq))
    if (is.null(rule.which)) {
        rule <- as.matrix(rule.ante)[sort(sample(1:nrow(rule.ante), rule.num)),]
    } else {
        rule <- as.matrix(rule.ante)[rule.which,]
    }
    rule <- cbind(rule, 1:rule.num, 0, 1)
    colnames(rule) <- NULL

    fis <- addrule(fis, as.matrix(rule))
    fis
}

#' Produces an example fis object for Waiter-Tipping.
#'
#' A function used primarily for example purposes, it creates a fis with two input (service & food), output variables (tip) and their membership functions.
#'
#' @return A fis is return
#' @examples
#' fis <- tipper()
#' @export
tipper <- function() {
    fis = newfis('tipper', andMethod='prod')
    fis = addvar(fis, 'input', 'service', c(0, 10))
    fis = addvar(fis, 'input', 'food', c(0, 10))
    fis = addvar(fis, 'output', 'tip', c(0, 30))

    fis = addmf(fis, 'input', 1, 'poor', 'gaussmf', c(1.5, 0, 1))
    fis = addmf(fis, 'input', 1, 'good', 'gaussmf', c(1.5, 5, 1))
    fis = addmf(fis, 'input', 1, 'excellent', 'gaussmf', c(1.5, 10, 1))

    fis = addmf(fis, 'input', 2, 'rancid', 'trapmf', c(0, 0, 1, 3, 1))
    fis = addmf(fis, 'input', 2, 'delicious', 'trapmf', c(7, 9, 10, 10, 1))

    fis = addmf(fis, 'output', 1, 'cheap', 'trimf', c(0, 5, 10, 1))
    fis = addmf(fis, 'output', 1, 'average', 'trimf', c(10, 15, 20, 1))
    fis = addmf(fis, 'output', 1, 'generous', 'trimf', c(20, 25, 30, 1))

    rl = rbind(c(1, 1, 1, 1, 2), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 2))
    fis = addrule(fis, rl)
    fis
}


# testing
#' Produces an example fis object which can be used for ANFIS.
#'
#' A function used primarily for example purposes, it creates a fis with two input (service & food), output variables (tip) and their membership functions.
#'
#' @return A fis is return
#' @examples
#' fis <- anfis.tipper()
#' @export
anfis.tipper <- function() {
    fis = newfis('tipper', fisType = 'tsk', andMethod = 'prod', orMethod = 'max', impMethod = 'min', aggMethod = 'max')

    fis = addvar(fis, 'input', 'service', c(0, 10), 'singleton.fuzzification')
    #fis= addvar(fis, 'input', 'service', c(0, 10), 'it2gbell.fuzzification', c(1,2,3))

    fis = addvar(fis, 'input', 'food', c(0, 15), 'singleton.fuzzification')
    #fis= addvar(fis, 'input', 'food', c(0, 15), 'it2gbell.fuzzification', c(1,2,3))

    fis = addvar(fis, 'output', 'tip', c(0, 30), 'defuzzification.method', c(1, 2, 3))
    #fis= addvar(fis, 'output', 'test', c(0, 100))

    fis = addmf(fis, 'input', 1, 'poor', 'gbellmf', c(1, 1.5, 3))

    fis = addmf(fis, 'input', 1, 'good', 'gbellmf', c(1, 1.5, 5))
    #fis= addmf(fis, 'input', 1, 'good', 'it2gbellmf', c(1, 2, 1.5, 5))

    fis = addmf(fis, 'input', 1, 'excellent', 'gbellmf', c(2, 1.5, 10))
    #fis= addmf(fis, 'input', 1, 'excellent', 'it2gbellmf', c(1, 2, 1.5, 10))

    fis = addmf(fis, 'input', 2, 'rancid', 'gbellmf', c(1, 2, 3))
    #fis= addmf(fis, 'input', 2, 'rancid', 'it2gbellmf', c(1, 2, 2, 3))

    fis = addmf(fis, 'input', 2, 'delicious', 'gbellmf', c(4, 10, 10))
    #fis= addmf(fis, 'input', 2, 'delicious', 'it2gbellmf', c(3, 4, 10, 10))

    fis = addmf(fis, 'output', 1, 'cheap', 'linearmf', c(1, 5, 10))
    fis = addmf(fis, 'output', 1, 'average', 'linearmf', c(10, 15, 20))
    fis = addmf(fis, 'output', 1, 'generous', 'linearmf', c(20, 25, 30))
    #fis= addmf(fis, 'output', 1, 'generous', 'itlinearmf', c(10, 20, 25, 30))

    #fis= addmf(fis, 'output', 2, 'l', 'trimf', c(0, 0, 50))
    #fis= addmf(fis, 'output', 2, 'm', 'trimf', c(0, 50, 100))
    #fis= addmf(fis, 'output', 2, 'h', 'trimf', c(50, 100, 100))

    rl = rbind(c(1, 1, 1, 1, 1), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 1))
    #rl = rbind(c(1,1,1,1,1,2), c(2,0,2,2,1,1), c(3,2,3,3,1,2))
    fis = addrule(fis, rl)

    fis
}


#' Produces an example fis object (TSK type), which can also be optimised by ANFIS.
#'
#' A function used primarily for example purposes, it creates a fis with two input (service & food), output variables (tip) and their membership functions.
#'
#' @return A fis is return
#' @examples
#' fis <- tipper.tsk()
#' @export
tipper.tsk <- function() {
    fis = newfis('tipper', fisType = 'tsk', andMethod = 'prod')
    # fis = newfis('tipper', fisType = 'tsk')

    fis = addvar(fis, 'input', 'service', c(0, 10), 'singleton.fuzzification')
    # fis= addvar(fis, 'input', 'service', c(0, 10), 'it2gbell.fuzzification', c(1,2,3))

    fis = addvar(fis, 'input', 'food', c(0, 15), 'singleton.fuzzification')
    #fis= addvar(fis, 'input', 'food', c(0, 15), 'it2gbell.fuzzification', c(1,2,3))

    fis = addvar(fis, 'output', 'tip', c(0, 30))
    #fis= addvar(fis, 'output', 'test', c(0, 100))

    fis = addmf(fis, 'input', 1, 'poor', 'gbellmf', c(1, 1.5, 3, 1))

    fis = addmf(fis, 'input', 1, 'good', 'gbellmf', c(1, 1.5, 5, 1))
    #fis= addmf(fis, 'input', 1, 'good', 'it2gbellmf', c(1, 2, 1.5, 5, 1))

    fis = addmf(fis, 'input', 1, 'excellent', 'gbellmf', c(2, 1.5, 10, 1))
    #fis= addmf(fis, 'input', 1, 'excellent', 'it2gbellmf', c(1, 2, 1.5, 10, 1))

    fis = addmf(fis, 'input', 2, 'rancid', 'gbellmf', c(1, 2, 3, 1))
    #fis= addmf(fis, 'input', 2, 'rancid', 'it2gbellmf', c(1, 2, 2, 3, 1))

    fis = addmf(fis, 'input', 2, 'delicious', 'gbellmf', c(4, 10, 10, 1))
    #fis= addmf(fis, 'input', 2, 'delicious', 'it2gbellmf', c(3, 4, 10, 10, 1))

    fis = addmf(fis, 'output', 1, 'cheap', 'linearmf', c(0, 0.5, 0.5))
    fis = addmf(fis, 'output', 1, 'average', 'linearmf', c(10, 0.5, 0.5))
    fis = addmf(fis, 'output', 1, 'generous', 'linearmf', c(20, 0.5, 0.5))
    #fis= addmf(fis, 'output', 1, 'generous', 'itlinearmf', c(20, 25, 0.5, 0.5))

    rl = rbind(c(1, 1, 1, 1, 1), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 1))
    #rl = rbind(c(1,1,1,1,1,2), c(2,0,2,2,1,1), c(3,2,3,3,1,2))
    fis = addrule(fis, rl)

    fis
}


#' Produces an example it2fis object for Waiter-Tipping.
#'
#' A function used primarily for example purposes, it creates a it2 fis with two input (service & food), output variables (tip) and their membership functions.
#'
#' @return A fis object
#' @examples
#' it2fis <- it2tipper()
#' @export
it2tipper <- function() {

    fis = newfis('tipper', mfType = 'it2')

    fis = addvar(fis, 'input', 'service', c(0, 10), 'singleton.fuzzification')
    fis = addvar(fis, 'input', 'food', c(0, 15), 'it2gbell.fuzzification', c(1, 2, 3))
    fis = addvar(fis, 'output', 'tip', c(0, 30), 'defuzzification.method')

    fis = addmf(fis, 'input', 1, 'poor', 'it2gbellmf', c(1, 2, 1.5, 0, 1, 1))
    fis = addmf(fis, 'input', 1, 'good', 'it2gbellmf', c(1, 2, 1.5, 5, 1, 1))
    fis = addmf(fis, 'input', 1, 'excellent', 'it2gbellmf', c(1, 2, 1.5, 10, 1, 1))

    fis = addmf(fis, 'input', 2, 'rancid', 'it2gbellmf', c(1, 2, 2, 3, 1, 1))
    fis = addmf(fis, 'input', 2, 'delicious', 'it2gbellmf', c(1, 2, 2, 12, 1, 1))

    fis = addmf(fis, 'output', 1, 'cheap', 'it2trimf', c(3, 7.5, 12, 0, 7.5, 15, 0.8, 1))
    fis = addmf(fis, 'output', 1, 'average', 'it2trimf', c(10.5, 15, 19.5, 7.5, 15, 22.5, 0.8, 1))
    fis = addmf(fis, 'output', 1, 'generous', 'it2trimf', c(18, 22.5, 27, 15, 22.5, 30, 0.8, 1))

    rl = rbind(c(1, 1, 1, 1, 2), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 2))
    fis = addrule(fis, rl)

    fis
}


#' Produces an example non-singleton fis object for Waiter-Tipping.
#'
#' A function used primarily for example purposes, it creates a nsfis with two input (service & food), output variables (tip) and their membership functions.
#' @return A non-singleton fis object
#' @examples
#' fis <- tipper.ns()
#' @author Yu Zhao
#' @export
tipper.ns <- function() {
    fis = newfis("tipper")
    fis = addvar(fis, "input", "service", c(0, 10), 'gauss', 0.5, 'tnorm.min.max')
    fis = addvar(fis, "input", "food", c(0, 10), 'gauss', 0.5, 'tnorm.min.defuzz.centroid')
    fis = addvar(fis, "output", "tip", c(0, 30))

    fis = addmf(fis, "input", 1, "poor", "gaussmf", c(1.5, 0, 1))
    fis = addmf(fis, "input", 1, "good", "gaussmf", c(1.5, 5, 1))
    fis = addmf(fis, "input", 1, "excellent", "gaussmf", c(1.5, 10, 1))

    fis = addmf(fis, "input", 2, "rancid", "trapmf", c(0, 0, 1, 3, 1))
    fis = addmf(fis, "input", 2, "delicious", "trapmf", c(7, 9, 10, 10, 1))

    fis = addmf(fis, "output", 1, "cheap", "trimf", c(0, 5, 10, 1))
    fis = addmf(fis, "output", 1, "average", "trimf", c(10, 15, 20, 1))
    fis = addmf(fis, "output", 1, "generous", "trimf", c(20, 25, 30, 1))
    rl = rbind(c(1, 1, 1, 1, 2), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 2))
    fis = addrule(fis, rl)

    fis
}


#' Graphic User Interface for Waiter-Tipping
#'
#' Graphic User Interface for Waiter-Tipping to display the membership function (input & output) and rules.
#'
#' @return Return graphic user interface for Waiter-Tipping
#' @author Tajul Razak
#' @examples
#' fis <- tipperGUI()
#' @import shiny splines
#' @export
tipperGUI = function() {

    ui = fluidPage(
        titlePanel(title = h1("Type-1 Fuzzy Logic : Waiter-Tipping", align = "center")),
        br(),
        sidebarLayout(
            sidebarPanel((h3("Parameter ")),
                         selectInput("service", "Choose a MF type for input 'Service':",
                                     list("Trapeziod", "Gaussian", "Tringular")),
                         selectInput("food", "Choose a MF type for input 'Food':",
                                     list("Trapeziod", "Gaussian", "Tringular")),
                         selectInput("tip", "Choose a MF type for output 'Tip':",
                                     list("Trapeziod", "Gaussian", "Tringular")),
                         radioButtons("out", "Display MF ",
                                      list("Service" = 1,
                                           "Food" = 2,
                                           "Tip" = 3))
            ),
            mainPanel(
                tabsetPanel(
                    type = "tab",
                    tabPanel("Membership Function", plotOutput("plot1")),
                    tabPanel("Rules", verbatimTextOutput("out"))
                    #  tabPanel("Output", uiOutput("output1"))
                )
            )
        )
    )


    server = function(input, output) {

        output$plot1 <- renderPlot({

            fis = newfis('tipper')

            fis = addvar(fis, 'input', 'service', c(0, 10))
            fis = addvar(fis, 'input', 'food', c(0, 10))
            fis = addvar(fis, 'output', 'tip', c(0, 30))

            # MF type for service
            if (input$service == "Trapeziod")
            {
                fis = addmf(fis, 'input', 1, 'poor', 'trapmf', c(-2.844, -0.6886, 0.6886, 2.844))
                fis = addmf(fis, 'input', 1, 'good', 'trapmf', c(2.156, 4.311, 5.689, 7.844))
                fis = addmf(fis, 'input', 1, 'excellent', 'trapmf', c(7.156, 9.311, 10.69, 12.84))
            }
            else if (input$service == "Gaussian")
            {
                fis = addmf(fis, 'input', 1, 'poor', 'gaussmf', c(1.5, 0))
                fis = addmf(fis, 'input', 1, 'good', 'gaussmf', c(1.5, 5))
                fis = addmf(fis, 'input', 1, 'excellent', 'gaussmf', c(1.5, 10))
            }
            else if (input$service == "Tringular")
            {
                fis = addmf(fis, 'input', 1, 'poor', 'trimf', c(-3.532, 0, 3.532))
                fis = addmf(fis, 'input', 1, 'good', 'trimf', c(1.468, 5, 8.532))
                fis = addmf(fis, 'input', 1, 'excellent', 'trimf', c(6.468, 10, 13.53))
            }

            # MF type for food
            if (input$food == "Trapeziod")
            {
                fis = addmf(fis, 'input', 2, 'rancid', 'trapmf', c(0, 0, 1, 3))
                fis = addmf(fis, 'input', 2, 'delicious', 'trapmf', c(7, 9, 10, 10))
            }
            else if (input$food == "Gaussian")
            {
                fis = addmf(fis, 'input', 2, 'rancid', 'gaussmf', c(1.189, 0.6))
                fis = addmf(fis, 'input', 2, 'delicious', 'gaussmf', c(1.172, 9.4))
            }
            else if (input$food == "Tringular")
            {
                fis = addmf(fis, 'input', 2, 'rancid', 'trimf', c(-2.2, 0.6001, 3.4))
                fis = addmf(fis, 'input', 2, 'delicious', 'trimf', c(6.639, 9.4, 12.16))
            }

            # MF type for tip
            if (input$tip == "Trapeziod")
            {
                fis = addmf(fis, 'output', 1, 'cheap', 'trapmf', c(0.5, 4.5, 5.5, 9.5))
                fis = addmf(fis, 'output', 1, 'average', 'trapmf', c(10.5, 14.5, 15.5, 19.5))
                fis = addmf(fis, 'output', 1, 'generous', 'trapmf', c(20.5, 24.5, 25.5, 29.5))
            }
            else if (input$tip == "Gaussian")
            {
                fis = addmf(fis, 'output', 1, 'cheap', 'gaussmf', c(2.123, 5))
                fis = addmf(fis, 'output', 1, 'average', 'gaussmf', c(2.123, 15))
                fis = addmf(fis, 'output', 1, 'generous', 'gaussmf', c(2.123, 25))
            }
            else if (input$tip == "Tringular")
            {
                fis = addmf(fis, 'output', 1, 'cheap', 'trimf', c(0, 5, 10))
                fis = addmf(fis, 'output', 1, 'average', 'trimf', c(10, 15, 20))
                fis = addmf(fis, 'output', 1, 'generous', 'trimf', c(20, 25, 30))
            }


            #add rules
            rl = rbind(c(1, 1, 1, 1, 2), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 2))
            fis = addrule(fis, rl)

            # showrules(fis)

            if (input$out == 1)
            {
                plotmf(fis, "input", 1, main = "Membership function plots")
            }
            else if (input$out == 2)
            {
                plotmf(fis, "input", 2, main = "Membership function plots")
            }
            else if (input$out == 3)
            {
                plotmf(fis, "output", 1, main = "Membership function plots")
            }
        })

        output$out = renderPrint({
            fis = tipper()
            showrule(fis)
            #  evalfis(c(1,2,5), fis)
        })
    }


    shinyApp(ui = ui, server = server)
}

#' Graphic User Interface for Waiter-Tipping (another style)
#'
#' Another style of Graphic User Interface for Waiter-Tipping to display the membership function (input & output) and rules.
#'
#' @return Return graphic user interface for Waiter-Tipping
#' @author Tajul Razak
#' @examples
#' fis <- tipperGUI2()
#' @import shiny splines
#' @export
tipperGUI2 = function() {


    ui = fluidPage(
        titlePanel(title = h1("Type-1 Fuzzy Logic : Waiter-Tipping", align = "center")),
        br(),
        sidebarLayout(
            sidebarPanel(
                radioButtons("out", "Please Choose Variable : ",
                             list("Service" = 1,
                                  "Food" = 2,
                                  "Tip" = 3)),
                conditionalPanel(
                    condition = "input.out == 1",
                    selectInput("service", "Choose a MF type for input 'Service':",
                                list("Trapeziod", "Gaussian", "Tringular")),
                    conditionalPanel(
                        condition = "input.service == 'Trapeziod'",
                        selectInput("out1", "Linguistic Variable':",
                                    list("Poor", "Good", "Excellent")),
                        conditionalPanel(
                            condition = "input.out1 == 'Poor'",
                            sliderInput("a1", "A", min = 0, max = 10, value = -2.844, step = 0.1),
                            sliderInput("b1", "B", min = 0, max = 10, value = -0.6886, step = 0.1),
                            sliderInput("c1", "C", min = 0, max = 10, value = 0.6886, step = 0.1),
                            sliderInput("d1", "D", min = 0, max = 10, value = 2.844, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out1 == 'Good'",
                            sliderInput("a2", "A", min = 0, max = 10, value = 2.156, step = 0.1),
                            sliderInput("b2", "B", min = 0, max = 10, value = 4.311, step = 0.1),
                            sliderInput("c2", "C", min = 0, max = 10, value = 5.689, step = 0.1),
                            sliderInput("d2", "D", min = 0, max = 10, value = 7.844, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out1 == 'Excellent'",
                            sliderInput("a3", "A", min = 0, max = 10, value = 7.156, step = 0.1),
                            sliderInput("b3", "B", min = 0, max = 10, value = 9.311, step = 0.1),
                            sliderInput("c3", "C", min = 0, max = 10, value = 10.69, step = 0.1),
                            sliderInput("d3", "D", min = 0, max = 10, value = 12.84, step = 0.1)
                        )
                    ),
                    conditionalPanel(
                        condition = "input.service == 'Gaussian'",
                        selectInput("out2", "Linguistic Variable':",
                                    list("Poor", "Good", "Excellent")),
                        conditionalPanel(
                            condition = "input.out2 == 'Poor'",
                            sliderInput("g1", "Mean", min = 0, max = 10, value = 1.5, step = 0.1),
                            sliderInput("g2", "Spread", min = 0, max = 10, value = 0, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out2 == 'Good'",
                            sliderInput("g3", "Mean", min = 0, max = 10, value = 1.5, step = 0.1),
                            sliderInput("g4", "Spread", min = 0, max = 10, value = 5, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out2 == 'Excellent'",
                            sliderInput("g5", "Mean", min = 0, max = 10, value = 1.5, step = 0.1),
                            sliderInput("g6", "Spread", min = 0, max = 10, value = 10, step = 0.1)
                        )
                    ),
                    conditionalPanel(
                        condition = "input.service == 'Tringular'",
                        selectInput("out3", "Linguistic Variable':",
                                    list("Poor", "Good", "Excellent")),
                        conditionalPanel(
                            condition = "input.out3 == 'Poor'",
                            sliderInput("t1", "Start", min = 0, max = 10, value = -3.532, step = 0.1),
                            sliderInput("t2", "Peak", min = 0, max = 10, value = 0, step = 0.1),
                            sliderInput("t3", "Stop", min = 0, max = 10, value = 3.532, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out3 == 'Good'",
                            sliderInput("t4", "Start", min = 0, max = 10, value = 1.468, step = 0.1),
                            sliderInput("t5", "Peak", min = 0, max = 10, value = 5, step = 0.1),
                            sliderInput("t6", "Stop", min = 0, max = 10, value = 8.532, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out3 == 'Excellent'",
                            sliderInput("t7", "Start", min = 0, max = 10, value = 6.468, step = 0.1),
                            sliderInput("t8", "Peak", min = 0, max = 10, value = 10, step = 0.1),
                            sliderInput("t9", "Stop", min = 0, max = 10, value = 13.53, step = 0.1)
                        )
                    )
                ),
                conditionalPanel(
                    condition = "input.out == 2",
                    selectInput("food", "Choose a MF type for input 'Food':",
                                list("Trapeziod", "Gaussian", "Tringular")),
                    conditionalPanel(
                        condition = "input.food == 'Trapeziod'",
                        selectInput("out4", "Linguistic Variable':",
                                    list("Rancid", "Delicious")),
                        conditionalPanel(
                            condition = "input.out4 == 'Rancid'",
                            sliderInput("a4", "A", min = 0, max = 10, value = 0, step = 0.1),
                            sliderInput("b4", "B", min = 0, max = 10, value = 0, step = 0.1),
                            sliderInput("c4", "C", min = 0, max = 10, value = 1, step = 0.1),
                            sliderInput("d4", "D", min = 0, max = 10, value = 3, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out4 == 'Delicious'",
                            sliderInput("a5", "A", min = 0, max = 10, value = 7, step = 0.1),
                            sliderInput("b5", "B", min = 0, max = 10, value = 9, step = 0.1),
                            sliderInput("c5", "C", min = 0, max = 10, value = 10, step = 0.1),
                            sliderInput("d5", "D", min = 0, max = 10, value = 10, step = 0.1)
                        )
                    ),
                    conditionalPanel(
                        condition = "input.food == 'Gaussian'",
                        selectInput("out5", "Linguistic Variable':",
                                    list("Rancid", "Delicious")),
                        conditionalPanel(
                            condition = "input.out5 == 'Rancid'",
                            sliderInput("g7", "Mean", min = 0, max = 10, value = 1.189, step = 0.1),
                            sliderInput("g8", "Spread", min = 0, max = 10, value = 0.6, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out5 == 'Delicious'",
                            sliderInput("g9", "Mean", min = 0, max = 10, value = 1.172, step = 0.1),
                            sliderInput("g10", "Spread", min = 0, max = 10, value = 9.4, step = 0.1)
                        )
                    ),
                    conditionalPanel(
                        condition = "input.food == 'Tringular'",
                        selectInput("out6", "Linguistic Variable':",
                                    list("Rancid", "Delicious")),
                        conditionalPanel(
                            condition = "input.out6 == 'Rancid'",
                            sliderInput("t10", "Start", min = 0, max = 10, value = -2.2, step = 0.1),
                            sliderInput("t11", "Peak", min = 0, max = 10, value = 0.6001, step = 0.1),
                            sliderInput("t12", "Stop", min = 0, max = 10, value = 3.4, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out6 == 'Delicious'",
                            sliderInput("t13", "Start", min = 0, max = 10, value = 6.639, step = 0.1),
                            sliderInput("t14", "Peak", min = 0, max = 10, value = 9.4, step = 0.1),
                            sliderInput("t15", "Stop", min = 0, max = 10, value = 12.16, step = 0.1)
                        )
                    )
                ),
                conditionalPanel(
                    condition = "input.out == 3",
                    selectInput("tip", "Choose a MF type for output 'Tip':",
                                list("Trapeziod", "Gaussian", "Tringular")),
                    conditionalPanel(
                        condition = "input.tip == 'Trapeziod'",
                        selectInput("out7", "Linguistic Variable':",
                                    list("Cheap", "Average", "Generous")),
                        conditionalPanel(
                            condition = "input.out7 == 'Cheap'",
                            sliderInput("a6", "A", min = 0, max = 30, value = 0.5, step = 0.1),
                            sliderInput("b6", "B", min = 0, max = 30, value = 4.5, step = 0.1),
                            sliderInput("c6", "C", min = 0, max = 30, value = 5.5, step = 0.1),
                            sliderInput("d6", "D", min = 0, max = 30, value = 9.5, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out7 == 'Average'",
                            sliderInput("a7", "A", min = 0, max = 30, value = 10.5, step = 0.1),
                            sliderInput("b7", "B", min = 0, max = 30, value = 14.5, step = 0.1),
                            sliderInput("c7", "C", min = 0, max = 30, value = 15.5, step = 0.1),
                            sliderInput("d7", "D", min = 0, max = 30, value = 19.5, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out7 == 'Generous'",
                            sliderInput("a8", "A", min = 0, max = 30, value = 20.5, step = 0.1),
                            sliderInput("b8", "B", min = 0, max = 30, value = 24.5, step = 0.1),
                            sliderInput("c8", "C", min = 0, max = 30, value = 25.5, step = 0.1),
                            sliderInput("d8", "D", min = 0, max = 30, value = 29.5, step = 0.1)
                        )
                    ),
                    conditionalPanel(
                        condition = "input.tip == 'Gaussian'",
                        selectInput("out8", "Linguistic Variable':",
                                    list("Cheap", "Average", "Generous")),
                        conditionalPanel(
                            condition = "input.out8 == 'Cheap'",
                            sliderInput("g11", "Mean", min = 0, max = 30, value = 2.123, step = 0.1),
                            sliderInput("g12", "Spread", min = 0, max = 30, value = 5, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out8 == 'Average'",
                            sliderInput("g13", "Mean", min = 0, max = 30, value = 2.123, step = 0.1),
                            sliderInput("g14", "Spread", min = 0, max = 30, value = 15, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out8 == 'Generous'",
                            sliderInput("g15", "Mean", min = 0, max = 30, value = 2.123, step = 0.1),
                            sliderInput("g16", "Spread", min = 0, max = 30, value = 25, step = 0.1)
                        )
                    ),
                    conditionalPanel(
                        condition = "input.tip == 'Tringular'",
                        selectInput("out9", "Linguistic Variable':",
                                    list("Cheap", "Average", "Generous")),
                        conditionalPanel(
                            condition = "input.out9 == 'Cheap'",
                            sliderInput("t16", "Start", min = 0, max = 30, value = 0, step = 0.1),
                            sliderInput("t17", "Peak", min = 0, max = 30, value = 5, step = 0.1),
                            sliderInput("t18", "Stop", min = 0, max = 30, value = 10, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out9 == 'Average'",
                            sliderInput("t19", "Start", min = 0, max = 30, value = 10, step = 0.1),
                            sliderInput("t20", "Peak", min = 0, max = 30, value = 15, step = 0.1),
                            sliderInput("t21", "Stop", min = 0, max = 30, value = 20, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.out9 == 'Generous'",
                            sliderInput("t22", "Start", min = 0, max = 30, value = 20, step = 0.1),
                            sliderInput("t23", "Peak", min = 0, max = 30, value = 25, step = 0.1),
                            sliderInput("t24", "Stop", min = 0, max = 30, value = 30, step = 0.1)
                        )
                    )
                )
            ),
            mainPanel(
                tabsetPanel(
                    type = "tab",
                    tabPanel("Membership Function", plotOutput("plot1")),
                    tabPanel("Rules", verbatimTextOutput("out"))
                    #  tabPanel("Output", uiOutput("output1"))
                )
            )
        )
    )


    server = function(input, output) {

        output$plot1 <- renderPlot({

            fis = newfis('tipper')

            fis = addvar(fis, 'input', 'service', c(0, 10))
            fis = addvar(fis, 'input', 'food', c(0, 10))
            fis = addvar(fis, 'output', 'tip', c(0, 30))

            # MF type for service
            if (input$service == "Trapeziod")
            {
                fis = addmf(fis, 'input', 1, 'poor', 'trapmf', c(input$a1, input$b1, input$c1, input$d1))
                fis = addmf(fis, 'input', 1, 'good', 'trapmf', c(input$a2, input$b2, input$c2, input$d2))
                fis = addmf(fis, 'input', 1, 'excellent', 'trapmf', c(input$a3, input$b3, input$c3, input$d3))
            }
            else if (input$service == "Gaussian")
            {
                fis = addmf(fis, 'input', 1, 'poor', 'gaussmf', c(input$g1, input$g2))
                fis = addmf(fis, 'input', 1, 'good', 'gaussmf', c(input$g3, input$g4))
                fis = addmf(fis, 'input', 1, 'excellent', 'gaussmf', c(input$g5, input$g6))
            }
            else if (input$service == "Tringular")
            {
                fis = addmf(fis, 'input', 1, 'poor', 'trimf', c(input$t1, input$t2, input$t3))
                fis = addmf(fis, 'input', 1, 'good', 'trimf', c(input$t4, input$t5, input$t6))
                fis = addmf(fis, 'input', 1, 'excellent', 'trimf', c(input$t7, input$t8, input$t9))
            }

            # MF type for food
            if (input$food == "Trapeziod")
            {
                fis = addmf(fis, 'input', 2, 'rancid', 'trapmf', c(input$a4, input$b4, input$c4, input$d4))
                fis = addmf(fis, 'input', 2, 'delicious', 'trapmf', c(input$a5, input$b5, input$c5, input$d5))
            }
            else if (input$food == "Gaussian")
            {
                fis = addmf(fis, 'input', 2, 'rancid', 'gaussmf', c(input$g7, input$g8))
                fis = addmf(fis, 'input', 2, 'delicious', 'gaussmf', c(input$g9, input$g10))
            }
            else if (input$food == "Tringular")
            {
                fis = addmf(fis, 'input', 2, 'rancid', 'trimf', c(input$t10, input$t11, input$t12))
                fis = addmf(fis, 'input', 2, 'delicious', 'trimf', c(input$t13, input$t14, input$t15))
            }

            # MF type for tip
            if (input$tip == "Trapeziod")
            {
                fis = addmf(fis, 'output', 1, 'cheap', 'trapmf', c(input$a6, input$b6, input$c6, input$d6))
                fis = addmf(fis, 'output', 1, 'average', 'trapmf', c(input$a7, input$b7, input$c7, input$d7))
                fis = addmf(fis, 'output', 1, 'generous', 'trapmf', c(input$a8, input$b8, input$c8, input$d8))
            }
            else if (input$tip == "Gaussian")
            {
                fis = addmf(fis, 'output', 1, 'cheap', 'gaussmf', c(input$g11, input$g12))
                fis = addmf(fis, 'output', 1, 'average', 'gaussmf', c(input$g13, input$g14))
                fis = addmf(fis, 'output', 1, 'generous', 'gaussmf', c(input$g15, input$g16))
            }
            else if (input$tip == "Tringular")
            {
                fis = addmf(fis, 'output', 1, 'cheap', 'trimf', c(input$t16, input$t17, input$t18))
                fis = addmf(fis, 'output', 1, 'average', 'trimf', c(input$t19, input$t20, input$t21))
                fis = addmf(fis, 'output', 1, 'generous', 'trimf', c(input$t22, input$t23, input$t24))
            }


            #add rules
            rl = rbind(c(1, 1, 1, 1, 2), c(2, 0, 2, 1, 1), c(3, 2, 3, 1, 2))
            fis = addrule(fis, rl)

            # showrules(fis)

            if (input$out == 1)
            {
                plotmf(fis, "input", 1, main = "Membership function plots")
            }
            else if (input$out == 2)
            {
                plotmf(fis, "input", 2, main = "Membership function plots")
            }
            else if (input$out == 3)
            {
                plotmf(fis, "output", 1, main = "Membership function plots")
            }
        })

        output$out = renderPrint({
            fis = tipper()
            showrule(fis)
            #  evalfis(c(1,2,5), fis)
        })
    }


    shinyApp(ui = ui, server = server)
}

#' Show a Graphic User Interface of fis object
#'
#' Show a Graphic User Interface to display membership function plots for input and output, rules and evaluate the fis.
#'
#' @param fis Requires a fis structure to display a GUI.
#' @param advancedGUI TRUE/FALSE; if TRUE, an advanced GUI with more features is provided (provided by science@sboldt.com). 
#' @details
#' This function is purposed to display all the membership plots and rules of fis object in Graphic User Interface (GUI). It also provide a function to evaluate the fis object.\cr
#'
#' showGUI(fis) will display the GUI of fis object.
#' @return Return the GUI to display membership function for input and output together with rules.
#' @examples
#' fis <- tipper()
#' fis <- showGUI(fis)
#' @author Tajul Razak
#' @import shiny splines plyr
#' @export
showGUI = function(fis, advancedGUI=FALSE){
    # restructured by science@sboldt.com
    # date: 28th April 2021 - second edition
    
    ## set variables
    # set: NumOutputs, NumInput, LabelList, LabelList2, MinMaxList
    NumOutputs = length(fis$output)
    NumInput = length(fis$input)
    
    # LabelList is list of labels
    LabelList = vector("list", 10)
    
    # MinMaxList contains lower and upper limit of input values
    MinMaxList = vector("list", 10)    
    
    # give some init values
    for (a in 1:length(MinMaxList)) {
        LabelList[[a]] = "NULL"
        MinMaxList[[a]][1] = 0
        MinMaxList[[a]][2] = 1
    }
    
    # import values
    for (a in 1:NumInput) {
        LabelList[[a]] = fis$input[[a]]$name
        MinMaxList[[a]] = fis$input[[a]]$range
    }
    
    # LabelList2 is LabelList + label of output MF
    LabelList2 = LabelList[1:NumInput]
    LabelList2[[NumInput + 1]] = fis$output[[1]]$name
    
    
    # Number of rules:
    nRules <- nrow(fis$rule)
    
    # nRules5 is used by subpage to show max 5 rules at once: quotient and remainder
    nRules5 <- c(nRules %/% 5, nRules %% 5)
    
    # list of rules
    lRules <- vector("list")
    lRules_label <- vector("character")
    
    rule_index <- 0
    
    # split up rules into smaller lists with max. five elements
    # value nRules5[[1]] contains the info how many lists with 5 elements are needed (quotient)
    while (rule_index < nRules5[[1]]) {
        lRules[[rule_index + 1]] <- c((rule_index * 5 + 1):(rule_index * 5 + 5))
        lRules_label[[rule_index + 1]] <- paste(min(lRules[[rule_index + 1]]), max(lRules[[rule_index + 1]]), sep = '-')
        rule_index <- rule_index + 1
    }
    
    # value nRules5[[2]] contains the info how long the list with less than 5 elements is (remainder)
    if (nRules5[[2]]) {
        lRules[[rule_index + 1]] <- c((rule_index * 5 + 1):(rule_index * 5 + nRules5[[2]]))
        lRules_label[[rule_index + 1]] <- paste(min(lRules[[rule_index + 1]]), max(lRules[[rule_index + 1]]), sep = '-')
    }

    ui = fluidPage(
        title = paste("Type-1 Fuzzy Logic : ", fis$name),
        titlePanel(title = h1(
            paste("Type-1 Fuzzy Logic : ", fis$name), align = "center"
        )),
        br(),
        sidebarLayout(
            sidebarPanel(
                # allow user to see input information
                selectInput("NumInput", "Number of input':", list(NumInput)),
                selectInput("out1", "Number of output':", list(NumOutputs)),
                selectInput("out2", "Select input / output variable", LabelList2),
                
                # 
                radioButtons("eva", "Evaluate FIS ", list( "Reset" = 1, "evalfis()" = 2)),
                
                # add another slider if statement is true
                # maybe this can also be done as list/array
                conditionalPanel(condition = "input.eva == 2", sliderInput("slider1", LabelList[[1]], min = MinMaxList[[1]][1], max = MinMaxList[[1]][2], value = MinMaxList[[1]][1], step = 0.1),
                                 conditionalPanel(condition = "input.NumInput >= 2", sliderInput("slider2", LabelList[[2]], min = MinMaxList[[2]][1], max = MinMaxList[[2]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 3", sliderInput("slider3", LabelList[[3]], min = MinMaxList[[3]][1], max = MinMaxList[[3]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 4", sliderInput("slider4", LabelList[[4]], min = MinMaxList[[4]][1], max = MinMaxList[[4]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 5", sliderInput("slider5", LabelList[[5]], min = MinMaxList[[5]][1], max = MinMaxList[[5]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 6", sliderInput("slider6", LabelList[[6]], min = MinMaxList[[6]][1], max = MinMaxList[[6]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 7", sliderInput("slider7", LabelList[[7]], min = MinMaxList[[7]][1], max = MinMaxList[[7]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 8", sliderInput("slider8", LabelList[[8]], min = MinMaxList[[8]][1], max = MinMaxList[[8]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 9", sliderInput("slider9", LabelList[[9]], min = MinMaxList[[9]][1], max = MinMaxList[[9]][2], value = MinMaxList[[1]][1], step = 0.1)),
                                 conditionalPanel(condition = "input.NumInput >= 10", sliderInput("slider10", LabelList[[10]], min = MinMaxList[[10]][1], max = MinMaxList[[10]][2], value = MinMaxList[[1]][1], step = 0.1))
                ),
                
                # exit button to quit page
                actionButton("do", "Exit")
                
            ),
            if (advancedGUI){
                mainPanel(
                    tabsetPanel(
                        type = "tab",
                        
                        # Show uses membership function
                        tabPanel("Membership Function",plotOutput("plot1"), downloadButton(outputId = "plot1down", label = "Download the plot")),
                        
                        # Rules in words
                        tabPanel("Rules", verbatimTextOutput("RulesInWords")),
                        
                        # Show how much a single rule applies and contributes
                        tabPanel("Plot rules (single)",selectInput("rule", "Choose a rule:", c(1:nrow(fis$rule))), plotOutput("plotRules")),
                        
                        # Show how much each rule applies and contributes (all rules at once)
                        tabPanel("Plot rules (all)", plotOutput("plotRulesAll")),
                        
                        # Show how much each rule applies and contributes (max. rules at once)
                        tabPanel("Plot rules (chosen)",selectInput("rule2", "Choose a set rules:", lRules_label), plotOutput("plotRulesAll2")),
                        
                        # added rules as graph (like defuzzifier but single contributes by rule are visible)
                        tabPanel("Output-MF", plotOutput("plotRules2")),
                        
                        # Defuzzifier
                        tabPanel("Defuzzifier", plotOutput("defu"), downloadButton(outputId = "defudown", label = "Download the plot"))
                        #  tabPanel("Output", uiOutput("output1"))
                    )
                )
            } else { # if not advancedGUI
                mainPanel(
                    tabsetPanel(
                        type = "tab",
                        
                        # Show uses membership function
                        tabPanel("Membership Function",plotOutput("plot1"), downloadButton(outputId = "plot1down", label = "Download the plot")),
                        
                        # Rules in words
                        tabPanel("Rules", verbatimTextOutput("RulesInWords")),
                        
                        # Defuzzifier
                        tabPanel("Defuzzifier", plotOutput("defu"), downloadButton(outputId = "defudown", label = "Download the plot"))
                        #  tabPanel("Output", uiOutput("output1"))
                    )
                )
            
            } 
        )
    )

    server = function(input, output){
        
        
        # import input values:
        get_inputvalues <- function(input) {
            inputvalues <- c(input$slider1)
            # add another input if statement is true
            if (NumInput >= 2) {inputvalues <- c(inputvalues, input$slider2)}
            if (NumInput >= 3) {inputvalues <- c(inputvalues, input$slider3)}
            if (NumInput >= 4) {inputvalues <- c(inputvalues, input$slider4)}
            if (NumInput >= 5) {inputvalues <- c(inputvalues, input$slider5)}
            if (NumInput >= 6) {inputvalues <- c(inputvalues, input$slider6)}
            if (NumInput >= 7) {inputvalues <- c(inputvalues, input$slider7)}
            if (NumInput >= 8) {inputvalues <- c(inputvalues, input$slider8)}
            if (NumInput >= 9) {inputvalues <- c(inputvalues, input$slider9)}
            if (NumInput >= 10) {inputvalues <- c(inputvalues, input$slider10)}
            return(inputvalues)
        }
        
        ## get_inputvalues <- function(input){
        ##     inputvalues <- vector("list", 10)
        ##     return(inputvalues)}

        # Plot membership functions
        output$plot1 <- renderPlot({
            # importing values
            # name 'xx' is taken from plotmf: primary inputs for extra lines
            
            # if input value was chosen: plot input MF (with index 'a')
            for (a in 1:NumInput) {
                if ((advancedGUI) && (input$eva == 2)) {
                    xx =  get_inputvalues(input)[[a]]
                } else {
                    xx = NULL
                }
                
                if (input$out2 == LabelList2[[a]]) {plotmf(fis, "input", a, main = "Membership function plots", xx = xx)}
            }

            # if output is chosen as 'out2', just show membership functions
            if (input$out2 == fis$output[[1]]$name) {plotmf(fis, "output", 1, main = "Membership function plots")
            }
        })        

        # Download membership functions
        output$plot1down <- downloadHandler(
            filename = function() {
                paste("MFs", "pdf", sep = ".")
            },
            # content is a function with argument file. content writes the plot to the device
            content = function(file) {
                pdf(file)
                
                # if input value was chosen: plot input MF (with index 'a')
                for (a in 1:NumInput) {
                    if ((advancedGUI) && (input$eva == 2)) {
                        xx =  get_inputvalues(input)[[a]]
                    } else {
                        xx = NULL
                    }
                    
                    if (input$out2 == LabelList2[[a]]) {plotmf(fis, "input", a, main = "Membership function plots", xx = xx)}
                }
                
                # if output value was chosen: plot output MF (with index '1')       
                if (input$out2 == fis$output[[1]]$name) {plotmf(fis, "output", 1, main = "Membership function plots")}
                
                dev.off()
            }
        )

        # Defuzzification plot
        output$defu <- renderPlot({
            
            inputvalues <- get_inputvalues(input)
            out_name = fis$output[[1]]$name
            evalfis(inputvalues, fis)
            
            # blue line shows output MF
            plot(D_x, D_y, type = "l", main = c(out_name, D_out), col = "blue",lwd = 2, ylim = c(0, 1))
            
            if (advancedGUI){
                # red line show defuzzified value
                lines(c(D_out, D_out), c(-1, 2), col = "red", lty = "dashed")
            }
        })

        # Defuzzification plot (download)
        output$defudown <- downloadHandler(
            filename =  function() {paste("MFs", "pdf", sep = ".")},
            # content is a function with argument file. content writes the plot to the device
            content = function(file) {
                
                inputvalues <- get_inputvalues(input)
                out_name = fis$output[[1]]$name
                evalfis(inputvalues, fis)
                
                pdf(file)
                plot(D_x, D_y, type = "l", main = c(out_name, D_out), col = "blue", lwd = 2)
                if (advancedGUI){
                    # red line show defuzzified value
                    lines(c(D_out, D_out), c(-1, 2), col = "red", lty = "dashed")
                }
                dev.off()
            }
        )
        observeEvent(input$do, {
            stopApp()
        })

        # Rules in words
        output$RulesInWords = renderPrint({
            showrule(fis)
        })

        # Show how much a single rule applies and contributes
        output$plotRules <- renderPlot({
            
            inputvalues <- get_inputvalues(input)
            r <- as.numeric(input$rule)
            evalfis(inputvalues, fis)
            plotmf(fis, "output", 1, main = "Membership function plots")
            
            if (input$eva == 2) {
                ## lines(D_x,OUT_RULE_CONS[r,], ylim=c(0,1))
                polygon(c(D_x[1], D_x, D_x[length(D_x)]), c(0, OUT_RULE_CONS[r, ], 0), border = NA,col = 'blue')
            }
        })

        # Show how much each rule applies and contributes (all rules at once)
        output$plotRulesAll <- renderPlot({
            
            inputvalues <- get_inputvalues(input)
            
            ## out_name = fis$output[[1]]$name
            evalfis(inputvalues, fis)
            
            # print(OUT_RULE_CONS)
            ## plot(D_x,D_y, type="l", main=c(out_name,D_out), col="blue", lwd=2, ylim=c(0,1))
            ## lines(c(D_out,D_out),c(-1,2),col="red",lty="dashed");
            
            par(mfrow = c(nrow(OUT_RULE_CONS), 1), mar = c(1, 1, 1, 1))
            
            for (r in 1:nrow(OUT_RULE_CONS)) {
                # TODO: if (r>1){} # add extra space between plots
                plotmf(fis, "output", 1, main = paste("rule",as.character(r), sep=" "))#, ylim = c(0, 1))
                if (input$eva == 2) {
                    ## lines(D_x,OUT_RULE_CONS[r,], ylim=c(0,1))
                    polygon( c(D_x[1], D_x, D_x[length(D_x)]), c(0, OUT_RULE_CONS[r, ], 0), border = NA, col = 'blue')}
            }
        })

        # Show how much each rule applies and contributes (max. five rules at once)
        output$plotRulesAll2 <- renderPlot({
            
            inputvalues <- get_inputvalues(input)
            ## out_name = fis$output[[1]]$name
            evalfis(inputvalues, fis)
            
            rule_index = match(input$rule2, lRules_label)            
            
            # print(OUT_RULE_CONS)
            ## plot(D_x,D_y, type="l", main=c(out_name,D_out), col="blue", lwd=2, ylim=c(0,1))
            ## lines(c(D_out,D_out),c(-1,2),col="red",lty="dashed");
            
            par(mfrow = c(length(lRules[[rule_index]]), 1), mar = c(1, 1, 1, 1))
            
            for (r in lRules[[rule_index]]) {
                # TODO: if (r>1){} # add extra space between plots
                plotmf(fis,"output", 1, main = paste("rule",as.character(r), sep=" "))#, ylim = c(0, 1))
                if (input$eva == 2) {
                    ## lines(D_x,OUT_RULE_CONS[r,], ylim=c(0,1))
                    polygon( c(D_x[1], D_x, D_x[length(D_x)]), c(0, OUT_RULE_CONS[r, ], 0), border = NA, col = 'blue')}
            }
        })
        
        # added rules as graph (like defuzzifier but single contributes by rule are visible)
        output$plotRules2 <- renderPlot({
            
            inputvalues <- get_inputvalues(input)
            evalfis(inputvalues, fis)
            
            for (r in 1:nrow(OUT_RULE_CONS)) {
                if (r == 1) {
                    plot(D_x,OUT_RULE_CONS[r, ], ylim = c(0, 1), type = 'l', xlab=fis$output[[1]]$name ,ylab="Output MF") }
                else {
                    lines(D_x, OUT_RULE_CONS[r, ])}
            }
        })
    }
  
    shinyApp(ui = ui, server = server)
}

