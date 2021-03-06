#' @title Extract lift segments from an rpart object in a data.frame format
#'
#' @description \code{extract_segments} takes as
#' input a fitted \code{rpart} object using the lift_method
#' and returns the resulting segments in table form. See example
#' below for more details on it's usage.
#'
#' @param segment_tree an object of class \code{rpart} fitted with the method 
#' from \code{\link{import_lift_method}} and with \code{x = T} and \code{y = T}.
#' @param alpha optional alpha value for confidence intervals
#' @return a data.frame containing the resulting segments.
#' @example examples/segmenTree_example.R
#' @seealso \code{\link{import_lift_method}}
#' @export


extract_segments <- function(segment_tree, alpha = NULL){
  if(is.null(segment_tree$y)) stop("rpart-fit must be the result of an rpart call with y = T")
  if(is.null(segment_tree$x)) stop("rpart-fit must be the result of an rpart call with x = T")
  value_reduction <- function(df){
    if(df$eqn[1]=="="){
      if(nrow(df) == 1){
        ans <- paste0(df$val, collapse = ",")
      } else {
        ans <- paste0(Reduce(intersect, sapply(df$val, function(x) unlist(strsplit(x, ",")))), collapse = ",")
      }
      return(list(paste0("=", paste0("{", ans, "}"))))
    } else {
      df$val <- as.numeric(df$val)
      ans_low <- max(df$val[df$eqn == ">"])
      if(is.infinite(ans_low)) ans_low <- min(segment_tree$x[, colnames(segment_tree$x) == df$var[1]], na.rm = T)
      ans_high <- min(df$val[df$eqn == "<"])
      if(is.infinite(ans_high)) ans_high <- max(segment_tree$x[, colnames(segment_tree$x) == df$var[1]], na.rm = T)
      return(list(c(paste0(">=", ans_low), paste0("<", ans_high))))
    }
  }
  
  leaves <- segment_tree$frame
  leaves$node <- as.integer(row.names(leaves))
  row.names(leaves) <- NULL
  leaves <- leaves[leaves$var == "<leaf>", ]
  
  paths <- path.rpart(segment_tree, nodes = leaves$node, print.it = F)
  segments <- sapply(paths, function(x){
    x <- x[-1] # remove root
    var <- gsub("[<|>|\\=].*", "", x)
    val <- gsub(".*[<|>|\\=]", "", x)
    val[val == ""] <- "empty_string"
    eqn_ind <- regexpr(pattern = "[<|>|\\=]", text = x)
    eqn <- substring(x, first = eqn_ind, last = eqn_ind)
    ans <- data.frame(var, val, eqn, stringsAsFactors = F)
    ans2 <- lapply(split(x = ans, f = ans$var), value_reduction)
    ans3 <- unname(mapply(function(var_vec, ans_vec){
      paste0(paste(var_vec, unlist(ans_vec), sep = ""), collapse = ",")
    }, names(ans2), ans2))
    return(paste0(ans3, collapse = ", "))
  }
  )
  
  segments <- data.frame(segment = segments,
                         n = leaves$n,
                         lift = leaves$yval,
                         row.names = NULL)
  segments <- segments[order(segments$lift, decreasing = T), ]
  
  if(!is.null(alpha)){
    binary_y <- all(segment_tree$y[, 1] %in% c(0, 1))
    preds <- predict(segment_tree)
    unique_preds <- segments$lift
    pred_ind <- match(preds, unique_preds)
    CI <- tapply(segment_tree$y, c(pred_ind, pred_ind), function(y_mat){
      y_mat <- matrix(y_mat, ncol = 2)
      mean_y_treatment <- mean(y_mat[y_mat[, 2] == 1, 1])
      mean_y_control <- mean(y_mat[y_mat[, 2] == 0, 1])
      lift <- mean_y_treatment - mean_y_control
      treatment_cases <- sum(y_mat[, 2] == 1)
      control_cases <- sum(y_mat[, 2] == 0)
      if(binary_y){
        lift_sd <- sqrt(mean_y_treatment*(1 - mean_y_treatment)/treatment_cases +
                          mean_y_control*(1 - mean_y_control)/control_cases)
      } else {
        lift_sd <- sqrt(var(y_mat[y_mat[, 2] == 1, 1]) + 
                          var(y_mat[y_mat[, 2] == 0, 1]))
      }
      
      lift_lower = lift - qnorm(alpha/2, lower.tail = F)*lift_sd
      lift_upper = lift + qnorm(alpha/2, lower.tail = F)*lift_sd
      return(list(lift, lift_lower, lift_upper))
    }, simplify = F)
    
    segments$lift_lower <- sapply(CI, function(x) x[[2]])
    segments$lift_upper <- sapply(CI, function(x) x[[3]])
  }
  return(segments)
}