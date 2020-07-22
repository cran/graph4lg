#' Scaling function
#'
#' @description Scales values between 0 and 1
#' @param x Numeric or integer vector
#' @keywords internal
#' @export
#' @examples
#' x <- runif(min = 3, max = 15, n = 20)
#' x01 <- sc01(x)

sc01 <- function(x){
  x01 <- (x-min(x))/(max(x)-min(x))

  return(x01)
}

#' Vector of custom colors
#'
#' @description Vector of custom colors
#' @keywords internal
#' @export
#' @examples
#' mypalette[1]


mypalette <- c("#F2B950", "#1D72F5","#DF0101",
             "#FF9326","#A945FF", "#96A725",
             "#0089B2","#FDF060","#FFA6B2","#BFF217",
             "#60D5FD","#CC1577","#F2B950",
             "#7FB21D","#EC496F","#326397",
             "#B26314","#027368","#A4A4A4",
             "#610B5E", "red", "green", "darkblue",
             "orange", "blue", "grey", "black")
