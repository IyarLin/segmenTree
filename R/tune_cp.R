#' @title Tune cp hyper parameter using cross validation
#' @description This function uses cross validation to find the optimal cp.
#' The metric optimsed is the out of sample modified gini
#'
#' @details The baseline is measured as the actual lift in every test set
#' node. Nodes with too few observations (such that lift is nan) are regardedd as maximum error (0.25)

#' @return Mean out of sample error for every input cp value
#' @example examples/segmenTree_example.R
#' @param rpart_fit An object of class rpart
#' @param cp_num How many cp values to evaluate
#' @param train_frac fraction of observations to be used to train the model in each cross validation sample
#' @param M number of cross validation samples

#' @seealso \code{\link{import_lift_method}}
#'
#' @export

tune_cp <- function(rpart_fit, cp_num = 10, train_frac = 0.8, M = 10){
  if(is.null(rpart_fit$x)) stop("must run rpart with x = T")
  if(is.null(rpart_fit$y)) stop("must run rpart with y = T")
  dat <- data.frame(y = rpart_fit$y, rpart_fit$x)
  cp_vec <- rpart_fit$cptable[1:cp_num]
  ans <- matrix(nrow = M, ncol = length(cp_vec), dimnames = list(1:M, cp_vec))
  for(i in 1:nrow(ans)){
    train_ind <- sample.int(nrow(dat), nrow(dat) * train_frac)
    train_set <- dat[train_ind, ]
    test_set <- dat[-train_ind, ]
    for(j in 1:ncol(ans)){
      train_model <- rpart(y ~ ., data = train_set,
                           control = rpart.control(cp = cp_vec[j]),
                           method = lift_method,
                           parms = list())
      preds <- predict(train_model, test_set)
      unique_preds <- unique(preds)
      pred_ind <- match(preds, unique_preds)
      actuals <- tapply(test_set$y, c(pred_ind, pred_ind), function(y_mat){
        y_mat <- matrix(y_mat, ncol = 2)
        lift <- mean(y_mat[y_mat[, 2] == 1, 1]) - mean(y_mat[y_mat[, 2] == 0, 1])
        return(list(lift, nrow(y_mat)))
      }, simplify = T)
      agg_lift <- data.frame(n = sapply(actuals, function(x) x[[2]]),
                              lift_actual = sapply(actuals, function(x) x[[1]]),
                              pred = unique_preds)
      agg_lift$lift_actual[is.nan(agg_lift$lift_actual)] <- 0

      ans[i, j] <- weighted.mean(agg_lift$lift_actual * sign(agg_lift$pred), agg_lift$n)
    }
  }
  return(ans)
}
