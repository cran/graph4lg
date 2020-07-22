#' Compare the partition into modules of two graphs
#'
#' @description The function computes the Adjusted Rand Index (ARI) to
#' compare two graphs' partitions into modules or clusters more generally.
#' Both graphs must have the same number of nodes, but not necessarily the same
#' number of links. They must also have the same node names and in the
#' same order.
#'
#' @details This index takes values between -1 and 1. It measures how often
#' pairs of nodes pertaining to the same module in one graph also pertain to
#' the same module in the other graph.
#' Therefore, large values indicate that both partitions are similar.
#' The Rand Index can be defined as the frequency of agreement between two
#' classifications into discrete classes. It is the number of times a pair of
#' elements are classified into the same class or in two different classes
#' in both compared classifications, divided by the total number of possible
#' pairs of elements. The Rand Index is between 0 and 1 but its maximum value
#' depends on the number of elements. Thus, another 'adjusted' index was
#' created, the Adjusted Rand Index. According to the Hubert et
#' Arabie's formula, the ARI is computed as follows:
#' \eqn{ARI=\frac{Index - Expected index}{Maximum index - Expected index}}
#' where the values of Index, Expected index and Maximum index are computed
#' from a contingency table.
#' This function uses \code{adjustedRandIndex} from package \pkg{mclust} which
#' applies the Hubert and Arabie's formula for the ARI.
#' This function works for undirected graphs only.
#' @param x The first graph object
#' \itemize{
#' \item{If \code{mode = 'graph'} (default), \code{x} is a graph object of
#' class \code{igraph}.
#' Then, its nodes must have the same names as in graph \code{y}.}
#' \item{If \code{mode = 'data.frame'}, \code{x} refers to a column of
#' the \code{data.frame} 'data'.
#' Then \code{x} must be a character string indicating the name of the
#' column of 'data' with the modules' labels of the nodes in the first graph.
#' In that case, the column can be of class \code{numeric}, \code{character}
#' or \code{factor} but will be converted into a \code{numeric} vector
#' in any case.}
#' \item{If \code{mode = 'vector'}, \code{x} is a vector of
#' class \code{character}, \code{factor} or \code{numeric}.
#' In that case, it must have the same length as vector \code{y} and
#' will be converted into a \code{numeric} vector.}
#' }
#' @param y The second graph object
#' Same classes possible as for \code{x}. Must be of the same format as \code{x}
#' @param mode A character string indicating whether x and y are igraph objects,
#' vectors or columns from a data.frame. \code{mode} can be 'graph',
#' 'data.frame' or 'vector'.
#' @param algo (if x and y are igraph objects) A character string indicating the
#' algorithm used to create the modules with \pkg{igraph}.
#' \itemize{
#' \item{If \code{algo = 'fast_greedy'} (default),
#' function \code{cluster_fast_greedy}
#' from \pkg{igraph} is used (Clauset et al., 2004).}
#' \item{If \code{algo = 'walktrap'} (default), function \code{cluster_walktrap}
#' from \pkg{igraph} is used (Pons et Latapy, 2006) with
#' 4 steps (default options).}
#' \item{If \code{algo = 'louvain'}, function \code{cluster_louvain}
#' from \pkg{igraph} is used (Blondel et al., 2008).
#' In that case, the number of modules created in each graph is imposed.}
#' \item{If \code{algo = 'optimal'}, function \code{cluster_optimal}
#' from \pkg{igraph} is used (Brandes et al., 2008) (can be very long).
#' In that case, the number of modules created in each graph is imposed.}
#' }
#' @param node_inter (optional, if x and y are igraph objects,
#' default is 'none') A character string indicating whether the links of the
#' graph are weighted by distances or by similarity indices. It is only used
#' to compute the modularity index. It can be: \itemize{
#' \item{'distance': Link weights correspond to distances. Nodes that are close
#' to each other will more likely be in the same module.}
#' \item{'similarity': Link weights correspond to similarity indices. Nodes that
#' are similar to each other will more likely be in the same module. Inverse
#' link weights are then used to compute the modularity index.}
#' \item{'none': Links are not weighted for the computation, which is only
#' based on graph topology.}
#' }
#' Two different weightings can be used to create the modules of the two graphs.
#' \itemize{
#' \item{If \code{node_inter} is a character string, then the same link
#' weighting is used for both graphs.}
#' \item{If \code{node_inter} is a character vector of length 2, then
#' the link weighting used by the algorithm to create the modules of
#' graphs \code{x} and \code{y} is determined by the first and second elements
#' of \code{node_inter}, respectively.}
#' }
#' @param nb_modul (if x and y are igraph objects) A numeric or integer value
#' or a numeric vector with 2 elements indicating the number of modules to
#' create in both graphs.
#' \itemize{
#' \item{If \code{nb_modul} is a numeric value, then the same number of modules
#' are created in both graphs.}
#' \item{If \code{nb_modul} is a numeric vector of length 2, then the
#' numbers of modules created in graphs \code{x} and \code{y} are the
#' first and second elements of \code{nb_modul}, respectively.}
#' }
#' @param data (if x and y are columns from a data.frame) An object of class
#' data.frame with at least two columns and as many rows as there are nodes
#' in the graphs compared. The columns indicate the modules of each node in
#' 2 different classifications.
#' @return The value of the ARI
#' @export
#' @author P. Savary
#' @references \insertRef{dyer2004population}{graph4lg}
#' \insertRef{hubert1985comparing}{graph4lg}
#' \insertRef{clauset2004finding}{graph4lg}
#' \insertRef{blondel2008fast}{graph4lg}
#' \insertRef{brandes2008modularity}{graph4lg}
#' \insertRef{pons2006computing}{graph4lg}
#' @examples
#' data(data_ex_genind)
#' data(pts_pop_ex)
#' mat_dist <- suppressWarnings(graph4lg::mat_geo_dist(data=pts_pop_ex,
#'       ID = "ID",
#'       x = "x",
#'       y = "y"))
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                       order(as.character(colnames(mat_dist)))]
#' graph_obs <- gen_graph_thr(mat_w = mat_dist, mat_thr = mat_dist,
#'                             thr = 24000, mode = "larger")
#' mat_gen <- mat_gen_dist(x = data_ex_genind, dist = "DPS")
#' graph_pred <- gen_graph_topo(mat_w = mat_gen, mat_topo = mat_dist,
#'                             topo = "gabriel")
#' ARI <- graph_modul_compar(x = graph_obs, y = graph_pred)

graph_modul_compar <- function(x,
                               y,
                               mode = "graph",
                               nb_modul = NULL,
                               algo = "fast_greedy",
                               node_inter = "distance",
                               data = NULL){

  if(mode == "graph"){

    # Check whether obs_graph and pred_graph are graphs
    if(!inherits(x, "igraph")){
      stop("'x' must be a graph object of class 'igraph'.")
    } else if (!inherits(y, "igraph")){
      stop("'y' must be a graph object of class 'igraph'.")
    }

    # Check whether x and y are undirected.
    if(any(c(igraph::is.directed(x), igraph::is.directed(y)))){
      stop("This function works for undirected graphs only")
    }

    # Check whether they have the same nodes' number
    if(length(igraph::V(x)) != length(igraph::V(y))){
      stop("Both graphs must have the same nodes' number.")
    }

    n_nodes <- length(igraph::V(x))

    # Check whether the graphs' nodes have names
    if(is.null(igraph::V(x)$name)){
      stop("The nodes of 'x' must have names.")
    } else if(is.null(igraph::V(y)$name)){
      stop("The nodes of 'y' must have names.")
    }

    # Check whether the graphs have the same nodes' names and in the same order
    if(!all(igraph::V(x)$name == igraph::V(y)$name)){
      stop("Both graphs must have the same nodes' names and the nodes
           ranked in the same order.")
    }

    if(!inherits(algo, "character")){
      stop("'algo' must be a character string.")
    } else if(!(algo %in% c("fast_greedy", "walktrap",
                            "louvain", "optimal"))){
      stop("'algo' must be either 'fast_greedy', 'walktrap',
            'louvain' or 'optimal'.")
    }


    ### Node interaction mode

    if(!inherits(node_inter, "character")){
      stop("'node_inter' must be a character string.")
    }

    # node_inter length 1

    if(length(node_inter) == 1){

      if(!(node_inter %in% c("distance", "similarity", "none"))){
        stop("'node_inter' must be 'distance', 'similarity' or 'none'.")
      }

      node_inter_x <- node_inter_y <- node_inter

      # node_inter length 2

    } else if (length(node_inter) == 2){

      if(!(node_inter[1] %in% c("distance", "similarity", "none"))){
        stop("Element 1 of 'node_inter' must be 'distance', 'similarity'
             or 'none'.")
      } else {
        node_inter_x <- node_inter
      }

      if(!(node_inter[2] %in% c("distance", "similarity", "none"))){
        stop("Element 2 of 'node_inter' must be 'distance', 'similarity'
             or 'none'.")
      } else {
        node_inter_y <- node_inter
      }

    } else {
      stop("'node_inter' must be of length 1 or 2.")
    }

    # Replace none by NULL to fit with compute_graph_modul
    # function argument

    if(node_inter_x == "none"){
      node_inter_x <- NULL
    }

    if(node_inter_y == "none"){
      node_inter_y <- NULL
    }

    ##### Nb_modul arguments

    if(is.null(nb_modul)){
      nb_modul_x <- nb_modul_y <- NULL
    } else {
      if(!inherits(nb_modul, c("numeric", "integer"))){
        stop("'nb_modul' must be a numeric or integer value or vector.")
      }
      if(length(nb_modul) == 1){
        nb_modul_x <- nb_modul_y <- nb_modul
      } else if(length(nb_modul) == 2){
        nb_modul_x <- nb_modul[1]
        nb_modul_y <- nb_modul[2]
      }
    }


    mx <- compute_graph_modul(graph = x, algo = algo,
                              node_inter = node_inter_x,
                              nb_modul = nb_modul_x)

    my <- compute_graph_modul(graph = y, algo = algo,
                              node_inter = node_inter_y,
                              nb_modul = nb_modul_y)


    x1 <- as.numeric(mx[, 'module'])
    y1 <- as.numeric(my[, 'module'])


  } else if (mode == "data.frame"){

    # Check whether data is a data.frame
    if(!inherits(data, "data.frame")){
      stop("'data' must be an object of class 'data.frame'.")
    }

    # Check whether x and y are character strings
    if(!inherits(x, "character")){
      stop("x must be of class 'character' when 'data' is a data.frame.")
    }
    if(!inherits(y, "character")){
      stop("y must be of class 'character' when 'data' is a data.frame.")
    }

    # Check whether x and y are valid names of columns from data
    # and convert the columns into numeric coding of the modules.
    if(all(c(x %in% names(data), y %in% names(data)))){
      x1 <- as.numeric(as.factor(data[, x]))
      y1 <- as.numeric(as.factor(data[, y]))
    } else {
      stop("x and y must be valid columns' names from 'data'.")
    }

  } else if (mode == "vector"){

    # Check whether x and y are vectors.
    if(!is.vector(x)){
      stop("x must be a vector when 'mode = vector'.")
    }
    if(!is.vector(y)){
      stop("y must be a vector when 'mode = vector'.")
    }
    # Check whether x and y are of same length
    if(length(x) != length(y)){
      stop("x and y must be of same length when 'mode = vector'.")
    }

    x1 <- as.numeric(as.factor(x))
    y1 <- as.numeric(as.factor(y))

  } else {
    stop("You must specify a correct 'mode' option")
  }

  x1 <- as.vector(x1)
  y1 <- as.vector(y1)
  tab <- table(x1,y1)

  # Number of modules in each graph.
  n_m1 <- dim(tab)[1]
  n_m2 <- dim(tab)[2]

  # Calculations from package 'mclust'

  if(all(dim(tab) == c(1,1))){
    return(1)

  } else {

    # \sum {n_{ij} \choose 2}
    a <- sum(choose(tab, 2))
    # \sum {a_i \choose 2} - a
    b <- sum(choose(rowSums(tab), 2)) - a
    # \sum {b_j \choose 2} - a
    c <- sum(choose(colSums(tab), 2)) - a
    # \sum {n \choose 2} - a - b - c
    d <- choose(sum(tab), 2) - a - b - c

    ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
      ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))

    print(paste(n_m1,"modules in graph 1 and", n_m2, "modules in graph 2",
                sep = " "))
    print(paste("Adjusted Rand Index: ", ARI, sep = ""))

    return(ARI)
  }
}





