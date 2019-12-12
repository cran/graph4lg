#' Compute dispersal kernel parameters
#'
#' @description The function computes the constant parameters of a dispersal
#' kernel with a negative exponential distribution
#'
#' @param p A numeric value indicating the dispersal probability at a distance
#' equal to 'd_disp' under a negative exponential distribution.
#' @param d_disp A numeric value indicating the distance to which dispersal
#' probability is equal to 'p' under a negative exponential distribution.
#' @param mode A character string indicating the value to return:
#' \itemize{
#' \item{If 'mode = 'A'' (default), the returned value 'alpha' is such that
#' exp(-alpha * d_disp) = p}
#' \item{If 'mode = 'B'', the returned value 'alpha' is such that
#' 10(-alpha * d_disp) = p}
#' }
#' @return A numeric value
#' @details If the resulting parameter when mode = "A" is a and the resulting
#' parameter when mode = "B" is b, then we have:
#' p = exp(-a.d_disp) = 10^(-b.d_disp) and a = b.ln(10)
#' @export
#' @author P. Savary
#' @examples
#' p <- 0.5
#' d_disp <- 3000
#' alpha <- kernel_param(p, d_disp, mode = "A")


kernel_param <- function(p, d_disp, mode = "A"){

  if(!inherits(p, c("numeric", "integer"))){
    stop("'p' must be of class 'numeric' or 'integer'")
  } else if(!inherits(d_disp, c("numeric", "integer"))){
    stop("'d_disp' must be of class 'numeric' or 'integer'")
  } else if(!(mode %in% c("A", "B"))){
    stop("'mode' must be equal to 'A' or 'B'")
  }

  A <- -log(p)/d_disp
  B <- A/log(10)

  #print("p(d_disp) = p = exp(-A * d_disp) = 10^(-B * d_disp) ssi A = B ln(10)")

  if( mode == "A"){
    print(paste("Returned value is A : ", A, sep = ""))
    return(A)
  } else if (mode == "B"){
    print(paste("Returned value is B : ", B, sep = ""))
    return(B)
  }

}



