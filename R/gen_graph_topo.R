#' Create a graph of genetic differentiation with
#' a specific topology
#'
#' @description The function constructs a genetic graph with
#' a specific topology from genetic and/or geographical distance matrices
#'
#' @param mat_w A symmetric (pairwise) \code{matrix} or a \code{dist} object
#' whose elements will be the links' weights
#' @param mat_topo (optional) A symmetric (pairwise) distance \code{matrix}
#' or a \code{dist} object whose values will be used for the pruning method.
#' @param topo Which topology does the created graph have?
#' \itemize{
#' \item{If 'topo = 'gabriel'' (default), the resulting graph will be a
#' Gabriel graph (Gabriel et al., 1969). It means that there is a link
#' between nodes x and y if and only if
#' \eqn{d_{xy}^{2} \leq \min(\sqrt{d_{xz}^{2}+d_{yz}^{2}}) },
#' with z any other node of the graph.}
#' \item{If 'topo = 'mst'', the resulting graph will have the topology
#' of a minimum spanning tree. It means that the graph will not include
#' any cycle (tree) and it will be the subgraph with a tree topology with
#' the minimum total links' weight (based on 'mat_topo' values).}
#' \item{If 'topo = 'percol'', if the link of the resulting graph with the
#' minimum weight is removed, then the graph breaks into two components.}
#' \item{If 'topo = 'comp'', a complete graph whose links are weighted with
#' values from 'mat_w' is created.}
#' \item{If 'topo = 'knn'', a k-nearest neighbor graph  whose links are
#' weighted with values from 'mat_w' is created. If the distance between node i
#' and node j is among the k-th smallest distances between node i and the other
#' nodes according to distances in matrix 'mat_topo', then there is a link
#' between i and j in the resulting graph. Therefore, a node can be connected
#' to more than two nodes because the nearest node to node j is not necessarily
#' among the k nearest neighbors to node i. Let d1 be the smallest distance
#' from node i to other nodes, if there are k nodes or more at this distance
#' from node i, they are all connected to i. If there are less than k nodes
#' connected to i at a distance d1, then we consider nodes at a distance d2
#' from i. In the latter case, all the nodes at a distance d2 from i are
#' connected to i.}
#' }
#' @param k (if 'topo = 'knn'') An integer which indicates the number of
#' nearest neighbors considered to create the K-nearest neighbor graph. k must
#' be lower than the total number of nodes minus 1.
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
#' Note that the topology 'knn' works best when 'mat_topo' contains distance
#' values from a continuous value range, thereby avoiding equal distances
#' between a node and the others.  are more than k nodes located
#' at distances in the k-th smallest distances
#' If dist objects are specified, it is assumed that colnames and
#' row.names of mat_w and mat_topo refer to the same populations/locations.
#' @examples
#' mat_w <- mat_gen_dist(x = data_ex_genind, dist = 'DPS')
#' suppressWarnings(mat_topo <- mat_geo_dist(pts_pop_ex,
#'                  ID = "ID",
#'                  x = "x",
#'                 y = "y"))
#' mat_topo <- mat_topo[row.names(mat_w), colnames(mat_w)]
#' graph <- gen_graph_topo(mat_w, mat_topo, topo = "mst")



gen_graph_topo <- function(mat_w, mat_topo = NULL, topo = "gabriel", k = NULL){

  # Check whether mat_topo is specified, else mat_w is used as mat_topo
  if(is.null(mat_topo)){
    mat_topo <- mat_w
  }

  # Check whether mat_w and mat_topo are symmetric matrices or dist objects
  if(!inherits(mat_w, c("matrix", "dist"))){
    stop("'mat_w' must be an object of class 'matrix' or 'dist'.")
  } else if (!inherits(mat_topo, c("matrix", "dist"))){
    stop("'mat_lc' must be an object of class 'matrix' or 'dist'.")
  } else if (inherits(mat_w, "matrix")){
    if(!Matrix::isSymmetric(mat_w)){
      stop("'mat_w' must be a symmetric pairwise matrix.")
    }
  } else if (inherits(mat_topo, "matrix")){
    if (!Matrix::isSymmetric(mat_topo)){
      stop("'mat_topo' must be a symmetric pairwise matrix.")
    }
  } else if (inherits(mat_w, "dist")){
    mat_w <- as.matrix(mat_w)
  } else if (inherits(mat_topo, "dist")){
    mat_topo <- as.matrix(mat_topo)
  }

  # Check whether mat_w and mat_topo have same dimensions and
  # same rows' and columns' names
  if(!all(dim(mat_w) == dim(mat_topo))){
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

    ### k-nearest neighbor graph #####

  } else if (topo == "knn"){


    # If k is null or not an integer and topo == "knn", returns an error
    if(is.null(mat_topo)){
      stop("If 'topo == 'knn', then k must be an integer")
    } else if(!inherits(k, c("integer", "numeric"))){
      stop("If 'topo == 'knn', then k must be an integer")
    }

    # Take an inteher value for k
    k <- as.integer(k)

    # Number of nodes
    n_pop <- nrow(mat_topo)

    if(k > (n_pop - 1)){
      stop("k is too high given the number of nodes. It must be lower than
           the number of nodes minus 1. Use 'topo == 'comp'' instead.")
    }

    M_adj <- mat_topo

    for(i in 1:nrow(M_adj)){

      rowi <- M_adj[i, ]

      dist_k <- rowi[order(rowi)][1:(k + 1)]

      rowi[-which(rowi %in% dist_k)] <- 0

      #rowi[order(rowi)][(k + 2):length(rowi)] <- 0

      M_adj[i, ] <- rowi
    }


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
