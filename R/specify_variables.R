#' Examine only therapeutic groups from Evil Insurers Marketscan sample
#'
#' This function allows you to screen out the therapeutic class variables when running
#' \code{SuperLearner()}. It relies on the therapeutic class variables all beginning with
#' the prefix \code{tcls}. I should probably make this a more general function that allows
#' the user to specify whatever variable prefix they want.
#' @param X data frame columns (should be specified when calling \code{SuperLearner()})
#' @keywords SuperLearner EvilInsurers
#' @export
#' @examples Cannot think of a good example here -- pretty much the one setting for using it
#'

tgrp.fun <- function(X, ...){
  whichvars <- c(rep.int(TRUE, ncol(X)))
  names(whichvars) <- colnames(X)
  tclsvars <- grep("tcls", names(X), value=T)
  whichvars[tclsvars] <- FALSE
  whichvars <- unname(whichvars)
  return(whichvars)
}


#' Lasso screener selects pre-specified variables
#'
#' Lasso screener for \code{SuperLearner()} that always retains specified variables and passes
#' approximately \code{nVar} variables to \code{SuperLearner()}. When the number of non-zero
#' coefficients exceeds \code{nVar}, a larger value of the regularization parameter \code{lambda}
#' is chosen to select a smaller set of variables that excludes the ties.
#'
#' @param X data frame
#' @param Y outcome variable (specified in \code{SuperLearner()})
#' @param var.index indices of variables to always be included by the screener
#' @param nVar number of non-zero variables to be selected
#'
#' @seealso See \code{\link[glmnet]{glmnet}} for additional details on implementing lasso
#'
#' @keywords SuperLearner EvilInsurers
#' @export
#' @examples If you do not know the indices of the variables you always want to include,
#'   you can get them from the variable name, where newdat is the dataframe:
#'   var.index <- c(which(colnames(newdat)=="tcls14"),
#'                  which(colnames(newdat)=="tcls251"))


screen.glmnet.fix <- function(Y, X, family, alpha = 1, minscreen = 2, nVar = 10, nfolds = 10, nlambda = 100,fixed.var.index=var.index,...) {
  if(!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = 'deviance',
                             nfolds = nfolds, family = family$family, alpha = alpha,
                             nlambda = nlambda, pmax= nVar, parallel=T)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  # the [-1] removes the intercept; taking the coefs from the fit w/ lambda that gives minimum cvm
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen,
            increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
  }
  whichVariable[c(var.index)] <- TRUE
  return(whichVariable)
}


#' Random forest screener selects pre-specified variables
#'
#' Random forest screener for \code{\link[SuperLearner]{SuperLearner()}} that selects
#' user specified variables in addition to variables chosen data-adaptively
#'
#' If you do not care about the exact number of variables the screener chooses, use this function
#' rather than \code{screen.rf.fix.exact}. This function is faster, but will not necessarily return
#' exactly \code{nVar} variables to \code{SuperLearner()}. \code{screen.rf.fuzzy} selects the top
#' \code{nVar} variables, and then also makes sure the user specified variables are also passed
#' to \code{SuperLearner()}. If the user specified variables are in the top \code{nVar} variables,
#' then \code{nVar} variables will be passed to \code{SuperLearner()}. If any of the user
#' specified variables are outside the top \code{nVar} variables, then more than \code{nVar}
#' variables will be passed to \code{SuperLearner()}.
#'
#' @section Super Learner:
#' See \code{SuperLearner()} documentation for information on additional arguments and
#' instructions on implementing \code{SuperLearner()}.
#'
#' @seealso \code{\link{screen.glmnet.fix}} for lasso screener, \code{\link{screen.rf.exact}}
#' for exact random forest screener.
#'
#'
#' @param X data frame
#' @param Y outcome variable (specified in \code{SuperLearner()})
#' @param var.index indices of variables to always be included by the screener
#' @param nVar number of variables for the screener to select
#' @keywords SuperLearner EvilInsurers
#' @export
#' @examples If you do not know the indices of the variables you always want to include,
#'  you can get them from the variable name, where newdat is the dataframe name:
#'
#'  var.index <- c(which(colnames(newdat)=="sex"), which(colnames(newdat)=="age"))

screen.rf.fuzzy <- function (Y, X, family, nVar = 10, ntree = 500, mtry = ifelse(family$family=="gaussian",
                                                                           floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), nodesize = ifelse(family$family=="gaussian", 5, 1), ...) {
  rank.rf.fit <- randomForest::randomForest(Y ~ ., data = X, ntree = ntree, mtry = mtry, nodesize = nodesize, keep.forest = FALSE)
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  whichVariable[c(var.index)] <- TRUE
  return(whichVariable)
}


#' Random forest screener limits selected variables and selects pre-specified variables
#'
#' Random forest screener for \code{\link[SuperLearner]{SuperLearner()}} that selects specified individual variables
#' and specified overall number of variables.
#'
#' This function can be pretty slow, because currently it
#' operates by searching the rankings for the user selected ("fixed") variables. If the fixed variables
#' are included in the top \code{nVar} then it does not change anything. If the fixed variables are
#' not included in the top \code{nVar}, then it selects a subset of top \code{nVar}; e.g., the overall number
#' of variables to select is 10 and 2 of the fixed variables are outside the top 10, it will select the
#' top 8, and convert the 2 fixed variables outside the top 10 to be \code{TRUE}.
#'
#' @section Super Learner:
#' See \code{SuperLearner()} documentation for information on additional arguments and
#' instructions on implementing \code{SuperLearner()}.
#'
#' @family other SLscreeners
#' @seealso \code{\link{screen.glmnet.fix}} for lasso screener
#'
#'
#' @param X data frame
#' @param Y outcome variable (specified in SuperLearner())
#' @param var.index indices of variables to always be included by the screener
#' @param nVar number of variables for the screener to select
#' @param nFix number of individual variables that are alaways passed to SuperLearner()
#' @keywords SuperLearner EvilInsurers
#' @export
#' @examples If you do not know the indices of the variables you always want to include,
#'  you can get them from the variable name, where newdat is the dataframe name:
#'
#'  var.index <- c(which(colnames(newdat)=="sex"), which(colnames(newdat)=="age"),
#'                which(colnames(newdat)=="emp_active"))
#'


screen.rf.exact <- function (Y, X, family, nVar = 10, nFix=3, fixed.var.index=var.index, ntree = 500, mtry = ifelse(family$family=="gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), nodesize = ifelse(family$family=="gaussian", 5, 1), ...)
{
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = X, ntree = ntree, mtry = mtry, nodesize = nodesize, keep.forest = FALSE)
  }

  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ ., data=X, ntree=ntree, mtry = mtry, nodesize = nodesize, keep.forest = FALSE)
  }
  ranks <- rank(-rank.rf.fit$importance)
  varInt <- (rank(-rank.rf.fit$importance) <= (nVar))
  varInt[c(fixed.var.index)] <- TRUE
  var_replace <- function(i){
    if(sum(varInt)==(nVar+i)){
      varFinal <- (rank(-rank.rf.fit$importance) <= (nVar-i))
      varFinal[c(fixed.var.index)] <- TRUE
      return(varFinal)
    }
  }
  results <- lapply(seq(0,nFix,1), var_replace) # applying the function over the specified number of fixed variables
  whichVariable <- unlist(Filter(Negate(is.null), results)) # getting rid of the null results and converting from a list to a vector
  return(whichVariable)
}

