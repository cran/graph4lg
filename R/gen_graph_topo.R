#' Create a graph of genetic differentiation with
#' a specific topology
#'
#' @description The function allows to construct a genetic graph with
#' a specific topology from genetic and/or geographical distance matrices
#'
#' @param mat_w A symmetric (pairwise) \code{matrix} whose elements
#' will be the links' weights
#' @param mat_topo (optional) A symmetric (pairwise) distance \code{matrix}
#' whose values will be used for the pruning method.
#' @param topo Which topology does the created graph have?
#' \itemize{
#' \item{If 'topo = 'gabriel'' (default), the resulting graph will be a Gabriel graph
#' (Gabriel et al., 1969). It means that there is a link between nodes x and y
#' if and only if \eqn{d_{xy}^{2} \leq \min(\sqrt{d_{xz}^{2}+d_{yz}^{2}}) },
#' with z any other node of the graph.}
#' \item{If 'topo = 'mst'', the resulting graph will have the topology
#' of a minimum spanning tree. It means that the graph will not include
#' any cycle (tree) and it will be the subgraph with a tree topology with
#' the minimum total links' weight (based on 'mat_topo' values).}
#' \item{If 'topo = 'percol'', if the link of the resulting graph with the
#' minimum weight is removed, then the graph breaks into two components.}
#' \item{If 'topo = 'comp'', a complete graph whose links are weighted with
#' values from 'mat_w' is created.}
#' }
#' @return A graph object of class \code{igraph}
#' @export
#' @author P. Savary
#' @references \insertRef{gabriel1969new}{graph4lg}
#' @details If 'mat_topo' is not defined, 'mat_w' is used for the pruning.
#' Matrices 'mat_w' and 'mat_topo' must have the same dimensions and the
#' same rows' and columns' names.
#' Values in 'mat_topo' matrix must be positive. Negative values from
#' 'mat_w' are transformed into zeros.
#' The function works only for undirected graphs.
#' @examples
#' mat_w <- mat_gen_dist(x = data_simul_genind, dist = 'DPS')
#' suppressWarnings(mat_topo <- mat_geo_dist(pts_pop_simul,
#'                  ID = "ID",
#'                  x = "x",
#'                 y = "y"))
#' mat_topo <- mat_topo[row.names(mat_w), colnames(mat_w)]
#' graph <- gen_graph_topo(mat_w, mat_topo, topo = "mst")


gen_graph_topo <- function(mat_w, mat_topo = NULL, topo = "gabriel"){

  # Check whether mat_topo is specified, else mat_w is used as mat_topo
  if(is.null(mat_topo)){
    mat_topo <- mat_w
  # Also check whether mat_w and mat_topo have same dimensions and
  # same rows' and columns' names
  } else if(!all(dim(mat_w) == dim(mat_topo))){
    stop("Matrices 'mat_w' and 'mat_topo' must have the same dimensions.")
  } else {
    if(!all(row.names(mat_w) == row.names(mat_topo))){
      stop("Matrices 'mat_w' and 'mat_topo' must have the same rows' names.")
    }

    if(!all(colnames(mat_w) == colnames(mat_topo))){
      stop("Matrices 'mat_w' and 'mat_topo' must have the same columns' names.")
    }
  }


  # Check whether there are negative values in mat_topo
  if(min(mat_topo) < 0){
    stop("There are negative values in 'mat_topo' matrix.")
  }

  # Transform negative values from mat_w into zeros
  if(min(mat_w) < 0){
    message("Negative values in 'mat_w' matrix were transformed into zeros.")
    mat_w[mat_w < 0] <- 0
  }


  ### Minimum spanning tree #####

  if (topo == "mst"){
    # Create a minimum spanning tree from mat_topo
    graph1 <- igraph::mst(igraph::graph.adjacency(as.matrix(mat_topo),
                                                  mode = "undirected",
                             weighted = TRUE, diag = FALSE))
    # Get the adjacency mtrix of graph1
    M_adj <- igraph::as_adj(graph1, sparse = FALSE)


  ### Gabriel graph ####

  } else if (topo == "gabriel"){

    # Link between x and y
    # ssi d_xy <= min(sqrt(d_xz^2+d_yz^2) | z \in S)
    M_adj <- matrix(0, nrow = nrow(mat_topo), ncol = ncol(mat_topo))
    row.names(M_adj) <- colnames(M_adj) <- row.names(mat_topo)

    # Apply the algorithm of the Gabriel graph to create an adjacency matrix
    # The loop runs along all the rows and columns of mat_topo
    for (i in 1:nrow(M_adj)){
      for (j in 1:nrow(M_adj)){
        # Only for the lower triangle of the matrix
        if(i < j){
          # Create an empty vector
          vect <- c()
          # For each row k, compute the sqrt of the sum of the squared distance
          # to i and to j
          for (k in 1:nrow(M_adj)){
            vect[k] <- sqrt( mat_topo[i, k]^2 + mat_topo[j, k]^2 )
          }
          # Remove the elements corresponding to i or j
          vect <- vect[-c(i, j)]
          # If the distance between i and j is lower or equal than the sum of
          # squared distances between i and j and any other point,
          # then there is a link between i and j
          if (mat_topo[i, j] <= min(vect)){
            M_adj[i, j] <- M_adj[j, i] <- 1
          } else {
            M_adj[i, j] <- M_adj[j, i] <- 0
          }

        }
      }
    }

  ### Percolation threshold #####

  } else if (topo == "percol"){

    # Use the function 'g_percol' to compute the graph graph1
    graph1 <- g_percol(mat_topo)

    # Get the adjacency matrix of graph1
    M_adj <- igraph::as_adj(graph1, sparse = FALSE)

  ### Complete graph #####

  } else if (topo == "comp"){

    # Create the full adjacency matrix with only 1
    n_pop <- nrow(mat_w)
    M_adj <- matrix(rep(1, n_pop * n_pop),
                    nrow = n_pop,
                    ncol = n_pop)

  } else {
    stop("You must specify a correct 'topo' option.")
  }

  # If there is no link between two points, the corresponding value
  # in the weighting matrix becomes also a zero
  mat_w[M_adj == 0] <- 0

  # Create the final graph from the modified mat_w matrix
  graph <- igraph::graph.adjacency(as.matrix(mat_w), mode = "undirected",
                           weighted = TRUE, diag = FALSE)

  return(graph)
}
