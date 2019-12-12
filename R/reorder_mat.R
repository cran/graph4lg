#' Reorder the rows and columns of a symmetric matrix
#'
#' @description The function reorders the rows and columns of a symmetric
#' matrix according to a specified order.
#'
#' @param mat An object of class \code{matrix}
#' @param order A character vector with the rows and columns names of the matrix
#' in the order in which they will be ordered by the function. All its elements
#' must be rows and columns names of the matrix \code{mat}.
#' @return A reordered symmetric matrix
#' @export
#' @author P. Savary
#' @details The matrix \code{mat} must be symmetric and have rows and columns
#' names. Its values are not modified.
#' @examples
#' mat <- matrix(rnorm(36), 6)
#' mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
#' row.names(mat) <- colnames(mat) <- c("A", "C", "E", "B", "D", "F")
#' order <- c("A", "B", "C", "D", "E", "F")
#' mat <- reorder_mat(mat = mat, order = order)


reorder_mat <- function(mat, order){

  # Number of elements in the vector 'order'
  n <- length(order)

  # Check whether 'mat' is a 'matrix'
  if(!inherits(mat, "matrix")){
    stop("'mat' must be a matrix")
  # Check whether 'order' is of class 'character'
  } else if (!inherits(order, "character")){
    stop("'order' must be a character vector")
  # Check whether 'mat' is a symmetric matrix
  } else if(!(isSymmetric(mat))){
    stop("The matrix 'mat' must be symmetric")
  # Check whether 'order' has as many elements as there are rows
  # and columns in 'mat'
  } else if (n != length(colnames(mat))){
    stop("'order' must have as many elements as there are rows and
         columns in 'mat'")
  # Check whether the column names are in the 'order' vector
  } else if(length(which(colnames(mat) %in% order)) != n){
    stop("The column names of the matrix you want to reorder must
         be present in the vector 'order'")
  # Check whether the row names are in the 'order' vector
  } else if (length(which(row.names(mat) %in% order)) != n){
    print("The row names of the matrix you want to reorder must
          be present in the vector 'order'")
  } else {

    # Reorder 'mat' according to 'order'
    mat2 <- mat[order, order]

    return(mat2)
  }
}



