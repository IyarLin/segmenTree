#' @title Import lift_method for lift modeling with rpart
#'
#' @description \code{import_lift_method} is a function that
#' imports a list of functions to serve as a user defined method
#' with the \code{rpart} function. see example below for more details on 
#' it's usage.
#' 
#' @param f_n a function that takes as input the split child node sample
#' size and returns a scalar/vector of the same length. this
#' number is used to weight the competing sub populations when
#' making the split. see example below for a use case. if NULL
#' the sample size is ignored when comparing splits.
#'
#' @return a list containing eval, split and init functions.
#' @example examples/segmenTree_example.R
#' @details the \code{rpart} function accepts in the method argument
#' a user defined list. This function imports the list that
#' implements a segment tree. The \code{y} variable in the input 
#' data.frame must be a 2 column matrix who's first column 
#' contains the response variable values (either binary or numeric, 
#' categorical isn't supported) and the second column is a
#' binary only treatment indicator.

#' @seealso \code{\link{extract_segments}}, \code{\link{tune_cp}}
#'
#' @export

import_lift_method <- function (f_n = NULL) 
{
  lift_method = list(eval = function(y, wt, parms) {
    mean_y_treatment <- mean(y[y[, 2] == 1, 1])
    mean_y_treatment <- replace(mean_y_treatment, 
                                       is.nan(mean_y_treatment), 0)
    mean_y_control <- mean(y[y[, 2] == 0, 1])
    mean_y_control <- replace(mean_y_control, 
                                     is.nan(mean_y_control), 0)

    lift <- mean_y_treatment - mean_y_control
    
    deviance <- nrow(y)^4
    
    list(label = lift, deviance = deviance)
  }, split = function(y, wt, x, parms, continuous) {
    n <- nrow(y)
    if(is.null(f_n)) f_n <- function(x) 1
    if (continuous) {
      positive_cases <- sum(y[, 1])
      treatment_cases <- sum(y[, 2])
      cases_left <- 1:(n - 1)
      cases_right <- (n - 1):1
      positive_treatment_left <- cumsum(y[, 1] * y[, 2])[-n]
      positive_treatment_right <- sum(y[, 1] * y[, 2]) - 
        positive_treatment_left
      treatment_n_left <- cumsum(y[, 2])[-n]
      treatment_n_right <- treatment_cases - treatment_n_left
      positive_control_left <- positive_cases - positive_treatment_left
      positive_control_right <- positive_cases - positive_treatment_right
      control_n_left <- n - treatment_n_left
      control_n_right <- n - treatment_n_right
      mean_y_treatment_left <- positive_treatment_left/treatment_n_left
      mean_y_control_left <- positive_control_left/control_n_left
      mean_y_treatment_left[is.nan(mean_y_treatment_left)] <- 0
      mean_y_control_left[is.nan(mean_y_control_left)] <- 0
      lift_left <- mean_y_treatment_left - mean_y_control_left
      lift_left[is.nan(lift_left)] <- 0
      lift_left_abs <- abs(lift_left)
      
      mean_y_treatment_right <- positive_treatment_right/treatment_n_right
      mean_y_control_right <- positive_control_right/control_n_right
      mean_y_treatment_right[is.nan(mean_y_treatment_right)] <- 0
      mean_y_control_right[is.nan(mean_y_control_right)] <- 0
      lift_right <- mean_y_treatment_right - mean_y_control_right
      lift_right[is.nan(lift_right)] <- 0
      lift_right_abs <- abs(lift_right)
      
      if(parms$preferred_lift == "both"){
        goodness <- pmax(lift_left_abs*f_n(cases_left), 
                         lift_right_abs*f_n(cases_right))
      } else if (parms$preferred_lift == "higher"){
        goodness <- pmax(lift_left*f_n(cases_left), 
                         lift_right*f_n(cases_right))
      } else {
        goodness <- pmax(-lift_left*f_n(cases_left), 
                         -lift_right*f_n(cases_right))
      }
      list(goodness = goodness, direction = sign(lift_left - lift_right))
    } else {
      cases_x <- tapply(y[, 1], x, length)
      positive_cases_x <- tapply(y[, 1], x, sum)
      treatment_cases_x <- tapply(y[, 2], x, sum)
      control_cases_x <- cases_x - treatment_cases_x
      
      positive_treatment_x <- tapply(y[, 1] * y[, 2], x, sum)
      mean_y_treatment_x <- positive_treatment_x/treatment_cases_x
      mean_y_treatment_x[is.nan(mean_y_treatment_x)] <- 0
      
      positive_control_x <- positive_cases_x - positive_treatment_x
      mean_y_control_x <- positive_control_x/control_cases_x
      mean_y_control_x[is.nan(mean_y_control_x)] <- 0
      
      lifts <- mean_y_treatment_x - mean_y_control_x
      ux <- sort(unique(x))
      ord <- order(abs(lifts))
      nx <- length(ux)
      
      cases_x_left <- cumsum(cases_x[ord])[-nx]
      cases_x_right <- (n - cases_x_left)[-nx]
      positive_cases <- sum(y[, 1])
      treatment_cases <- sum(y[, 2])
      positive_treatment_left <- cumsum(positive_treatment_x[ord])[-nx]
      positive_treatment_right <- sum(positive_treatment_x) - 
        positive_treatment_left
      treatment_n_left <- cumsum(treatment_cases_x[ord])[-nx]
      treatment_n_right <- sum(treatment_cases_x) - treatment_n_left
      mean_y_treatment_left <- positive_treatment_left/treatment_n_left
      mean_y_treatment_left[is.nan(mean_y_treatment_left)] <- 0
      mean_y_treatment_right <- positive_treatment_right/treatment_n_right
      mean_y_treatment_right[is.nan(mean_y_treatment_right)] <- 0
      positive_control_left <- cumsum(positive_control_x[ord])[-nx]
      positive_control_right <- sum(positive_control_x) - 
        positive_control_left
      control_n_left <- cases_x_left - treatment_n_left
      control_n_right <- cases_x_right - treatment_n_right
      mean_y_control_left <- positive_control_left/control_n_left
      mean_y_control_left[is.nan(mean_y_control_left)] <- 0
      mean_y_control_right <- positive_control_right/control_n_right
      mean_y_control_right[is.nan(mean_y_control_right)] <- 0
      lift_left <- mean_y_treatment_left - mean_y_control_left
      lift_left[is.nan(lift_left)] <- 0
      lift_left_abs <- abs(lift_left)
      lift_right <- mean_y_treatment_right - mean_y_control_right
      lift_right[is.nan(lift_right)] <- 0
      lift_right_abs <- abs(lift_right)
      
      if(parms$preferred_lift == "both"){
        goodness <- pmax(lift_left_abs*f_n(cases_x_left), 
                         lift_right_abs*f_n(cases_x_right))
      } else if (parms$preferred_lift == "higher"){
        goodness <- pmax(lift_left*f_n(cases_x_left), 
                         lift_right*f_n(cases_x_right))
      } else {
        goodness <- pmax(-lift_left*f_n(cases_x_left), 
                         -lift_right*f_n(cases_x_right))
      }

      list(goodness = goodness, direction = ux[ord])
    }
  }, init = function(y, offset, parms = NULL, wt) {
    if(is.null(parms$preferred_lift)) parms$preferred_lift <- "both"
    if(!parms$preferred_lift %in% c("both", "higher", "lower")) stop("parms$preferred_lift must be one of 'both', 'higher' or 'lower'")
    if (!is.matrix(y) | ncol(y) != 2) {
      stop("y has to be a 2 column matrix")
    }
    if (!all(y[, 2] %in% c(0, 1))) {
      stop("y[, 2] treatment column must be binary")
    }
    if (!missing(offset) && length(offset) > 0) {
      warning("offset argument ignored")
    }
    if(!identical(wt, rep(1, nrow(y)))) stop("segmenTree does not support weights argument wt")
    sfun <- function(yval, dev, wt, ylevel, digits) {
      paste(" lift=", format(signif(yval, digits)), ", deviance=", 
            format(signif(dev, digits)), sep = "")
    }
    environment(sfun) <- .GlobalEnv
    list(y = y, parms = parms, numresp = 1, numy = 2, summary = sfun)
  })
  return(lift_method)
}