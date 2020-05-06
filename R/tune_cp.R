#' @title Tune cp hyper parameter using cross validation
#' @description This function uses cross validation to find the optimal cp.
#'
#' @details For every cp value the tree is pruned accordingly. Next the data is
#' repeatedly split to train and test, where's the train is used to estimate lift
#' in each node. The metric calculated at each node is the multiplication of
#' the train and test lifts. If their signs agree a positive outcome is obtained. 
#' Otherwise, a negative one. 

#' @return a list `measures` with 2 elements: 
#' 
#' 1. cp_values X 2 matrix 
#' array\[,,1\] within node variance: 
#' array\[,,2\] contains between node variance
#' 
#' 1. Optimal cp found to maximize the ratio between the 2 measures
#' 
#' @example examples/segmenTree_example.R
#' @param segment_tree an object of class \code{rpart} fitted with the method 
#' from \code{\link{import_lift_method}} and \code{y = T}.
#' @param supress_plots should plots be suppressed?
#'
#' @export

tune_cp <- function(segment_tree, supress_plots = F){
  if(is.null(segment_tree$y)) stop("must run rpart with y = T")
  cp_vec <- segment_tree$cptable[, 1]
  ans <- matrix(nrow = length(cp_vec), ncol = 2, 
               dimnames = list(cp_vec, c("between node variance", "within node variance")))
  for(i in 1:nrow(ans)){
    pruned_tree <- prune(segment_tree, cp = cp_vec[i])
    leaves <- which(pruned_tree$frame$var == "<leaf>")
    res <- sapply(leaves, function(leaf){
      positive_rate_treatment <- mean(pruned_tree$y[pruned_tree$where == leaf & pruned_tree$y[, 2] == 1, 1])
      positive_rate_treatment <- replace(positive_rate_treatment, 
                                         is.nan(positive_rate_treatment), 0)
      treatment_cases <- sum(pruned_tree$where == leaf & pruned_tree$y[, 2] == 1)
      
      positive_rate_control <- mean(pruned_tree$y[pruned_tree$where == leaf & pruned_tree$y[, 2] == 0, 1])
      positive_rate_control <- replace(positive_rate_control, 
                                       is.nan(positive_rate_control), 0)
      control_cases <- sum(pruned_tree$where == leaf & pruned_tree$y[, 2] == 0)
      
      node_lift <- positive_rate_treatment - positive_rate_control
      
      node_variance <- positive_rate_treatment*(1-positive_rate_treatment)/treatment_cases + 
        positive_rate_control*(1-positive_rate_control)/control_cases
      
      node_size <- sum(pruned_tree$where == leaf)
      return(c(node_lift, node_variance, node_size))
    })
    ans[i, 1] <- var(res[1, ])*length(leaves)
    ans[i, 2] <- weighted.mean(res[2, ], res[3, ])
  }
  
  ans[is.na(ans)] <- 0
  
  if(!supress_plots){
    par(mfrow = c(1,3))
    plot(cp_vec, ans[, 1], main = "between node variance", ylab = "")
    plot(cp_vec, ans[, 2], main = "within node variance", ylab = "")
    plot(cp_vec, ans[, 1]/ans[, 2], main = "ratio", ylab = "")
  }
  
  return(list(measures = ans, optimal_cp = cp_vec[which.max(ans[, 1]/ans[, 2])]))
}
