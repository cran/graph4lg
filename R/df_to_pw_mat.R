#' Convert an edge-list data.frame into a pairwise matrix
#'
#' @description The function converts an edge-list data.frame
#' into a symmetric pairwise matrix
#'
#' @param data An object of class \code{data.frame}
#' @param from A character string indicating the name of the column with the ID
#' of the origins
#' @param to A character string indicating the name of the column with the ID
#' of the arrivals
#' @param value A character string indicating the name of the column with the
#' values corresponding to each pair
#' @details The matrix is a symmetric matrix. Be careful, you shall not provide
#' a data.frame with different values corresponding to the pair 1-2 and 2-1 as
#' an example. Ideally, for a complete matrix, data should have n(n-1)/2 rows
#' if values are computed between n objects.
#' @return A pairwise matrix
#' @export
#' @author P. Savary
#' @examples
#' data(pts_pop_simul)
#' suppressWarnings(mat_geo <- mat_geo_dist(pts_pop_simul,
#'                  ID = "ID",
#'                  x = "x",
#'                 y = "y"))
#' g <- gen_graph_topo(mat_w = mat_geo,
#'                     mat_topo = mat_geo,
#'                     topo = "comp")
#' df <- data.frame(igraph::as_edgelist(g))
#' df$w <- igraph::E(g)$weight
#' df_to_pw_mat(df, from = "X1", to = "X2", value = "w")

df_to_pw_mat <- function(data, from, to, value){

  if(!inherits(data, "data.frame")){
    stop("'data' must be of class 'data.frame'")
  } else if(!inherits(from, "character")){
    stop("'from' must be of class 'character'")
  } else if(!inherits(to, "character")){
    stop("'to' must be of class 'character'")
  } else if(!inherits(value, "character")){
    stop("'value' must be of class 'character'")
  } else if(!(from %in% colnames(data))){
    stop("There is no column named 'from' in 'data'")
  } else if(!(to %in% colnames(data))){
    stop("There is no column named 'to' in 'data'")
  } else if(!(value %in% colnames(data))){
    stop("There is no column named 'value' in 'data'")
  }

  # Get all the names
  names <- unique(c(as.character(data[, from]), as.character(data[, to])))
  # Number of different pops
  nb <- length(names)

  # Create a symmetric matrix
  mat <- matrix(NA, ncol = nb, nrow = nb)
  colnames(mat) <- row.names(mat) <- names

  # Fill the matrix with values from data
  for(row in 1:nrow(data)){

    i <- as.character(data[row, from])
    j <- as.character(data[row, to])
    val <- data[row, value]
    mat[i, j] <- mat[j, i] <- val

  }
  diag(mat) <- 0

  return(mat)
}
