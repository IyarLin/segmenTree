#' @title Find optimal cp parameter using elbow method
#' @description The function uses elbow strength to find optimal cp
#'
#' @details This function uses the method found 
#' \href{https://www.datasciencecentral.com/profiles/blogs/how-to-automatically-determine-the-number-of-clusters-in-your-dat}{here}
#' to determine the optimal cp. It also plots cp and elbow strength to allow the user to a manual determination.

#' @return optimal cp. It also plots cp and elbow strengh
#' @example examples/segmenTree_example.R
#' @param rpart_fit An object of class rpart
#' @param cp_num How many cp values to evaluate

#' @seealso \code{\link{tune_cp}}
#'
#' @export

cp_elbow <- function(rpart_fit, cp_num = 10){
  cp_vec <- rpart_fit$cptable[1:cp_num, 1]
  
  dd <- data.frame(Index = 1:length(cp_vec), cp = cp_vec)
  dd$delta1 <- c(NA, dd$cp[1:9] - dd$cp[2:10])
  dd$delta2 <- c(NA, NA, dd$delta1[1:8] - dd$delta1[2:9])
  dd$elbow_strengh <- c((dd$delta2 - dd$delta1)[-1], NA)
  
  par(mfrow = c(1, 2))
  plot(cp_vec, ylab = "cp")
  plot(dd$Index, dd$elbow_strengh, xlab = "Index", ylab = "Elbow strength")
  return(cp_vec[which.max(dd$elbow_strengh)])
}