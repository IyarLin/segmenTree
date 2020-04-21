#' @title Import lift_method for lift modeling with rpart
#'
#' @description \code{import_lift_method} is a function that
#' imports a list of functions to serve as a user defined method
#' with the rpart function.
#' See example below for more details on it's usage.
#'
#' @return A list containing eval, split and init functions.
#' @example examples/segmenTree_example.R
#' @details The rpart function accepts in the method argument
#' a user defined list. This function imports the list that
#' implements a causal tree. In addition to the method a user
#' also needs to input the parms argument with the baseline lift
#' in the population and the significance level for lift confidence
#' intervals returned in the yval2 object. See example below for more
#' details on how to use.
#' @seealso \code{\link{extract_segments}}
#'
#' @export

import_lift_method <- function(){
  lift_method = list(
    eval = function(y, wt, parms) {
      # y - matrix, first column is binary result, second is binary treatment
      positive_rate_treatment <- mean(y[y[, 2] == 1, 1])
      positive_rate_treatment <- replace(positive_rate_treatment, is.nan(positive_rate_treatment), 0)
      positive_rate_control <- mean(y[y[, 2] == 0, 1])
      positive_rate_control <- replace(positive_rate_control, is.nan(positive_rate_control), 0)
      treatment_cases <- sum(y[, 2] == 1)
      control_cases <- sum(y[, 2] == 0)

      lift <- positive_rate_treatment - positive_rate_control

      lift_to_p <- (lift + 1)/2
      deviance <- lift_to_p*(1-lift_to_p)*nrow(y)/parms$n



      list(label = lift, deviance = deviance)
    },
    split = function(y, wt, x, parms, continuous)
    {
      n <- nrow(y)
      if (continuous) {
        # continuous x variable
        positive_cases <- sum(y[, 1])
        treatment_cases <- sum(y[, 2])

        cases_left <- 1:(n-1)
        cases_right <- (n-1):1
        positive_treatment_left <- cumsum(y[, 1] * y[, 2])[-n]
        positive_treatment_right <- sum(y[, 1] * y[, 2]) - positive_treatment_left
        treatment_n_left <- cumsum(y[, 2])[-n]
        treatment_n_right <- treatment_cases - treatment_n_left

        positive_control_left <- positive_cases - positive_treatment_left
        positive_control_right <- positive_cases - positive_treatment_right
        control_n_left <- n - treatment_n_left
        control_n_right <- n - treatment_n_right

        positive_rate_treatment_left <- positive_treatment_left/treatment_n_left
        positive_rate_control_left <- positive_control_left/control_n_left

        positive_rate_treatment_left[is.nan(positive_rate_treatment_left)] <- 0
        positive_rate_control_left[is.nan(positive_rate_control_left)] <- 0

        lift_left <- positive_rate_treatment_left - positive_rate_control_left

        positive_rate_treatment_right <- positive_treatment_right/treatment_n_right
        positive_rate_control_right <- positive_control_right/control_n_right

        positive_rate_treatment_right[is.nan(positive_rate_treatment_right)] <- 0
        positive_rate_control_right[is.nan(positive_rate_control_right)] <- 0

        lift_right <- positive_rate_treatment_right - positive_rate_control_right

        lift_to_p_left <- (lift_left + 1)/2
        lift_to_p_right <- (lift_right + 1)/2
        goodness <- cases_left * (0.25 - lift_to_p_left*(1-lift_to_p_left)) +
          cases_right * (0.25 - lift_to_p_right*(1-lift_to_p_right))

        list(goodness = goodness, direction = sign(lift_left - lift_right))
      } else {
        # Categorical X variable
        cases_x <- tapply(y[, 1], x, length)
        positive_cases_x <- tapply(y[, 1], x, sum)
        treatment_cases_x <- tapply(y[, 2], x, sum)
        control_cases_x <- cases_x - treatment_cases_x

        positive_treatment_x <- tapply(y[, 1] * y[, 2], x, sum)
        positive_rate_treatment_x <- positive_treatment_x/treatment_cases_x
        positive_rate_treatment_x[is.nan(positive_rate_treatment_x)] <- 0

        positive_control_x <- positive_cases_x - positive_treatment_x
        positive_rate_control_x <- positive_control_x/control_cases_x
        positive_rate_control_x[is.nan(positive_rate_control_x)] <- 0

        lifts <- positive_rate_treatment_x - positive_rate_control_x

        ux <- sort(unique(x))
        ord <- order(abs(lifts))
        nx <- length(ux)

        cases_x_left <- cumsum(cases_x[ord])[-nx]
        cases_x_right <- (n - cases_x_left)[-nx]

        positive_cases <- sum(y[, 1])
        treatment_cases <- sum(y[, 2])

        positive_treatment_left <- cumsum(positive_treatment_x[ord])[-nx]
        positive_treatment_right <- sum(positive_treatment_x) - positive_treatment_left
        treatment_n_left <- cumsum(treatment_cases_x[ord])[-nx]
        treatment_n_right <- sum(treatment_cases_x) - treatment_n_left
        positive_rate_treatment_left <- positive_treatment_left/treatment_n_left
        positive_rate_treatment_left[is.nan(positive_rate_treatment_left)] <- 0
        positive_rate_treatment_right <- positive_treatment_right/treatment_n_right
        positive_rate_treatment_right[is.nan(positive_rate_treatment_right)] <- 0

        positive_control_left <- cumsum(positive_control_x[ord])[-nx]
        positive_control_right <- sum(positive_control_x) - positive_control_left
        control_n_left <- cases_x_left - treatment_n_left
        control_n_right <- cases_x_right - treatment_n_right
        positive_rate_control_left <- positive_control_left/control_n_left
        positive_rate_control_left[is.nan(positive_rate_control_left)] <- 0
        positive_rate_control_right <- positive_control_right/control_n_right
        positive_rate_control_right[is.nan(positive_rate_control_right)] <- 0

        lift_left <- positive_rate_treatment_left - positive_rate_control_left
        lift_right <- positive_rate_treatment_right - positive_rate_control_right

        lift_to_p_left <- (lift_left + 1)/2
        lift_to_p_right <- (lift_right + 1)/2
        goodness <- cases_x_left * (0.25 - lift_to_p_left*(1-lift_to_p_left)) +
          cases_x_right * (0.25 - lift_to_p_right*(1-lift_to_p_right))

        list(goodness = goodness, direction = ux[ord])
      }
    },
    init = function(y, offset, parms, wt) {
      if (!is.matrix(y) | ncol(y) != 2){
        stop("y has to be a 2 column matrix")
      }
      if (any(!y[, 1] %in% c(0, 1) | !y[, 2] %in% c(0, 1))){
        stop("y input columns must be binary")
      }
      if (!missing(offset) && length(offset) > 0){
        warning("offset argument ignored")
      }
      if(is.null(parms$alpha)) parms$alpha <- 0.05

      parms$baseline_lift <- mean(y[y[, 2] == 1, 1]) - mean(y[y[, 2] == 0, 1])
      parms$n <- nrow(y)
      sfun <- function(yval, dev, wt, ylevel, digits ) {
        paste(" lift=", format(signif(yval, digits)),
              ", deviance=" , format(signif(dev, digits)),
              sep = '')
      }
      environment(sfun) <- .GlobalEnv
      list(y = y, parms = parms, numresp = 1, numy = 2, summary = sfun)
    }
  )
  return(lift_method)
}
