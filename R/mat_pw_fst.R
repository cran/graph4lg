#' Compute a pairwise FST matrix between populations
#'
#' @description The function computes the pairwise FST matrix between
#' populations from an object of class \code{genind}
#'
#' @param x An object of class \code{genind}
#' @return A pairwise \code{matrix} of FST with as many rows and columns as
#' there are populations in the input data.
#' @details The formula used is that of Weir et Cockerham (1984).
#' This functions uses directly the function \code{pairwise.WCfst}
#' from \pkg{hierfstat}.
#' @section Warnings:
#' Negative values are converted into 0
#' @keywords internal
#' @export
#' @author P. Savary
#' @references \insertRef{weir1984estimating}{graph4lg}
#' @examples
#' \dontrun{
#' data("data_ex_genind")
#' mat_fst <- mat_pw_fst(data_ex_genind)
#' }

##################################

mat_pw_fst <- function(x){

  # If 'x' is a 'genind' object
  if (!inherits(x, "genind")){
    stop("Input value x must be a 'genind' object.")
  }

  # Compute the distance matrix between populations using hierfstat
  mat_fst <- hierfstat::pairwise.WCfst(hierfstat::genind2hierfstat(x))
  # Diagonal elements are 0
  diag(mat_fst) <- rep(0, length(diag(mat_fst)))

  return(mat_fst)
}




