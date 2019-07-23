#' Create a graph of genetic differentiation
#' using a link's weight threshold
#'
#' @description The function allows to construct a genetic graph whose
#' links' weights are larger or lower than a specific threshold
#'
#' @param mat_w A symmetric (pairwise) \code{matrix} whose elements
#' will be the links' weights
#' @param mat_thr (optional) A symmetric (pairwise) distance \code{matrix} whose
#' values will be used for the pruning based on the threshold value.
#' @param thr The threshold value (logically between min(mat_thr) and max(mat_thr))
#' @param mode
#' \itemize{
#' \item{If 'mode = 'larger'' (default), all the links whose weight is larger
#' than 'thr' are removed.}
#' \item{If 'mode = 'lower'', all the links whose weight is lower
#' than 'thr' are removed.}
#' }
#' @return A graph object of class \code{igraph}
#' @export
#' @author P. Savary
#' @details If 'mat_thr' is not defined, 'mat_w' is used for the pruning.
#' Matrices 'mat_w' and 'mat_thr' must have the same dimensions and the
#' same rows' and columns' names.
#' Values in 'mat_thr' matrix must be positive. Negative values from
#' 'mat_w' are transformed into zeros.
#' The function works only for undirected graphs.
#' @examples
#' mat_w <- mat_gen_dist(x = data_simul_genind, dist = 'DPS')
#' suppressWarnings(mat_thr <- mat_geo_dist(pts_pop_simul,
#'                  ID = "ID",
#'                  x = "x",
#'                 y = "y"))
#' mat_thr <- mat_thr[row.names(mat_w), colnames(mat_w)]
#' graph <- gen_graph_thr(mat_w, mat_thr, thr = 6000, mode = "larger")

gen_graph_thr <- function(mat_w, mat_thr = NULL, thr, mode = "larger"){

  # Check whether mat_thr is specified, else mat_w is used as mat_thr
  if(is.null(mat_thr)){
    mat_thr <- mat_w
  # Also check whether mat_w and mat_thr have same dimensions and
  # same rows' and columns' names
  } else if(!all(dim(mat_w) == dim(mat_thr))){
    stop("Matrices 'mat_w' and 'mat_thr' must have the same dimensions.")
  } else {
    if(!all(row.names(mat_w) == row.names(mat_thr))){
      stop("Matrices 'mat_w' and 'mat_thr' must have the same rows' names.")
    }
    if(!all(colnames(mat_w) == colnames(mat_thr))){
      stop("Matrices 'mat_w' and 'mat_thr' must have the same columns' names.")
    }
  }

  # Check whether there are negative values in mat_thr
  if(min(mat_thr) < 0){
    stop("There are negative values in 'mat_thr' matrix.")
  }

  # Transform negative values from mat_w into zeros
  if(min(mat_w) < 0){
    message("Negative values in 'mat_w' matrix were transformed into zeros.")
    mat_w[mat_w < 0] <- 0
  }

  # Check whether thr is a valid numerci value
  if(!is.numeric(thr)){
    stop("'thr' must be a numeric value")
  }

  # Depending on 'mode', remove links between some population pairs
  # according to the mat_thr values and to thr
  if(mode == "larger"){
    mat_thr[mat_thr > thr] <- 0
  } else if (mode == "lower"){
    mat_thr[mat_thr < thr] <- 0
  } else {
    stop("'mode' must be 'lower' or 'larger'.")
  }

  # Create a graph from the modified version of mat_thr
  graph1 <- igraph::graph.adjacency(as.matrix(mat_thr),
                                                mode = "undirected",
                                                weighted = TRUE, diag = FALSE)
  # get the adjacency matrix of graph1
  M_adj <- igraph::as_adj(graph1, sparse = FALSE)
  # The elements which are 0 in the adjacency matrix of graph1
  # become 0 in mat_w
  mat_w[M_adj == 0] <- 0
  # Create the final from the modified version of mat_w
  graph <- igraph::graph.adjacency(as.matrix(mat_w), mode = "undirected",
                           weighted = TRUE, diag = FALSE)

  return(graph)
}








