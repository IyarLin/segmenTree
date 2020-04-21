#' @title Extract lift segments from an rpart object in a table form
#'
#' @description \code{extract_segments} takes as
#' input a fitted rpart object using the lift_method
#' and returns the resulting segments in table form. See example
#' below for more details on it's usage.
#'
#' @param rpart_fit An object of class \code{rpart} fitted with the method from \code{\link{import_lift_method}} and with \code{y = T}.
#' @param alpha Optional alpha value for confiddence intervals
#' @return A data.frame containing the resulting segments.
#' It contains confidednce intervals with the alpha parameter specified in the \code{parms} argument to the \code{rpart} function.
#' @example examples/segmenTree_example.R
#' @seealso \code{\link{import_lift_method}}
#' @export


extract_segments <- function(rpart_fit, alpha = NULL){
  if(is.null(rpart_fit$y)) stop("rpart-fit must be the result of an rpart call with y = T")
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
      ans_high <- min(df$val[df$eqn == "<"])
      return(list(c(paste0(">=", ans_low), paste0("<", ans_high))))
    }
  }

  leaves <- rpart_fit$frame
  leaves$node <- as.integer(row.names(leaves))
  row.names(leaves) <- NULL
  leaves <- leaves[leaves$var == "<leaf>", ]

  paths <- path.rpart(rpart_fit, nodes = leaves$node, print.it = F)
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
    preds <- predict(rpart_fit)
    unique_preds <- segments$lift
    pred_ind <- match(preds, unique_preds)
    CI <- tapply(rpart_fit$y, c(pred_ind, pred_ind), function(y_mat){
      y_mat <- matrix(y_mat, ncol = 2)
      positive_rate_treatment <- mean(y_mat[y_mat[, 2] == 1, 1])
      positive_rate_control <- mean(y_mat[y_mat[, 2] == 0, 1])
      lift <- positive_rate_treatment - positive_rate_control
      treatment_cases <- sum(y_mat[, 2] == 1)
      control_cases <- sum(y_mat[, 2] == 0)
      lift_sd <- sqrt(positive_rate_treatment*(1 - positive_rate_treatment)/treatment_cases +
                        positive_rate_control*(1 - positive_rate_control)/control_cases)
      lift_lower = lift - qnorm(alpha/2, lower.tail = F)*lift_sd
      lift_upper = lift + qnorm(alpha/2, lower.tail = F)*lift_sd
      return(list(lift, lift_lower, lift_upper))
    }, simplify = F)

    segments$lift_lower <- sapply(CI, function(x) x[[2]])
    segments$lift_upper <- sapply(CI, function(x) x[[3]])
    if(!identical(unname(sapply(CI, function(x) x[[1]])), segments$lift)) stop("something is wrong with lift lower/upper. Do deep debugging")
  }
  return(segments)
}
