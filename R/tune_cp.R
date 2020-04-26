#' @title Tune cp hyper parameter using cross validation
#' @description This function uses cross validation to find the optimal cp.
#'
#' @details For every cp value the tree is pruned accordignly. Next the data is
#' repeateddly split to train and test, where's the train is used to estimate lift
#' in each node. The metric calculated at each noded is the multiplication of
#' the train and test lifts. If their signs agree a positive outcome is obtained. 
#' Otherwise, a negative one. 

#' @return A matrix containin in every row the sampling result for every cp in it's column
#' @example examples/segmenTree_example.R
#' @param rpart_fit An object of class rpart
#' @param cp_num How many cp values to evaluate
#' @param train_frac fraction of observations to be used to train the model in each cross validation sample
#' @param M number of cross validation rounds

#' @seealso \code{\link{cp_elbow}}
#'
#' @export

tune_cp <- function(rpart_fit, cp_num = 10, train_frac = 0.8, M = 10){
  if(is.null(rpart_fit$y)) stop("must run rpart with y = T")
  if(cp_num > nrow(rpart_fit$cptable)) cp_num <- nrow(rpart_fit$cptable)
  cp_vec <- rpart_fit$cptable[1:cp_num, 1]
  ans <- matrix(nrow = M, ncol = length(cp_vec), dimnames = list(1:M, cp_vec))
  for(j in 1:ncol(ans)){
    pruned_tree <- prune(rpart_fit, cp = cp_vec[j])
    leaves <- which(pruned_tree$frame$var == "<leaf>")
    for(i in 1:nrow(ans)){
      train_logic <- sample(c(rep(T, nrow(rpart_fit$y) * train_frac), 
                       rep(F, nrow(rpart_fit$y) - nrow(rpart_fit$y) * train_frac)))
      
      res <- sapply(leaves, function(leaf){
        train_est <- mean(pruned_tree$y[train_logic & pruned_tree$where == leaf & pruned_tree$y[, 2] == 1, 1]) - 
          mean(pruned_tree$y[train_logic & pruned_tree$where == leaf & pruned_tree$y[, 2] == 0, 1])
        test_est <- mean(pruned_tree$y[!train_logic & pruned_tree$where == leaf & pruned_tree$y[, 2] == 1, 1]) - 
          mean(pruned_tree$y[!train_logic & pruned_tree$where == leaf & pruned_tree$y[, 2] == 0, 1])
        return(c(train_est * test_est, sum(pruned_tree$where == leaf)))
      })
      res[is.nan(res)] <- 0
      ans[i, j] <- weighted.mean(res[1, ], res[2, ])
    }
  }
  return(ans)
}
