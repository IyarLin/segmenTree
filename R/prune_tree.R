#' @title Tune cp hyper parameter using cross validation
#' @description This function calculates within and between node 
#' variance in order to obtain optimal \code{cp} value.
#'
#' @details For every cp value a pruned tree is obtained. In every
#' the node lift and variance are calculated. next between and within
#' node variance are calculated. The optimal \code{cp} is the value
#' which maximizes the ratio between/within variance.

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

prune_tree <- function(segment_tree, supress_plots = F){
  if(is.null(segment_tree$y)) stop("must run rpart with y = T")
  cp_vec <- segment_tree$cptable[, 1]
  binary_y <- all(segment_tree$y[, 1] %in% c(0, 1))
  ans <- matrix(nrow = length(cp_vec), ncol = 2, 
               dimnames = list(cp_vec, c("between node variance", "within node variance")))
  for(i in 1:nrow(ans)){
    pruned_tree <- prune(segment_tree, cp = cp_vec[i])
    leaves <- which(pruned_tree$frame$var == "<leaf>")
    res <- sapply(leaves, function(leaf){
      mean_y_treatment <- mean(pruned_tree$y[pruned_tree$where == leaf & pruned_tree$y[, 2] == 1, 1])
      mean_y_treatment <- replace(mean_y_treatment, 
                                         is.nan(mean_y_treatment), 0)
      treatment_cases <- sum(pruned_tree$where == leaf & pruned_tree$y[, 2] == 1)
      
      mean_y_control <- mean(pruned_tree$y[pruned_tree$where == leaf & pruned_tree$y[, 2] == 0, 1])
      mean_y_control <- replace(mean_y_control, 
                                       is.nan(mean_y_control), 0)
      control_cases <- sum(pruned_tree$where == leaf & pruned_tree$y[, 2] == 0)
      
      node_lift <- mean_y_treatment - mean_y_control
      
      if(binary_y){
        node_variance <- mean_y_treatment*(1-mean_y_treatment)/treatment_cases + 
          mean_y_control*(1-mean_y_control)/control_cases
      } else { # continuous
        node_variance <- var(pruned_tree$y[pruned_tree$where == leaf & 
                                             pruned_tree$y[, 2] == 1, 1]) + 
          pruned_tree$y[pruned_tree$where == leaf & pruned_tree$y[, 2] == 0, 1]
      }
      
      
      
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
    par(mfrow=c(1,1))
  }
  
  return(list(measures = ans, optimal_cp = cp_vec[which.max(ans[, 1]/ans[, 2])]))
}
