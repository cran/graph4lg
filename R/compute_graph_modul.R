#' Compute modules from a graph by maximising modularity
#'
#' @description The function computes  modules from a graph by maximising
#' modularity.
#' @param graph An object of class \code{igraph}. Its nodes must have names.
#' @param algo A character string indicating the algorithm used to create
#' the modules with \pkg{igraph}. \itemize{
#' \item{If \code{algo = 'fast_greedy'} (default),
#' function \code{cluster_fast_greedy} from \pkg{igraph}
#' is used (Clauset et al., 2004).}
#' \item{If \code{algo = 'walktrap'}, function \code{cluster_walktrap}
#' from \pkg{igraph} is used (Pons et Latapy, 2006) with 4 steps
#' (default options).}
#' \item{If \code{algo = 'louvain'}, function \code{cluster_louvain}
#' from \pkg{igraph} is used (Blondel et al., 2008). In that case, the number
#' of modules created in each graph is imposed.}
#' \item{If \code{algo = 'optimal'}, function \code{cluster_optimal}
#' from \pkg{igraph} is used (Brandes et al., 2008) (can be very long).
#' In that case, the number of modules created in each graph is imposed.}
#' }
#' @param node_inter (optional, default = NULL) A character string indicating
#' whether the links of the graph are weighted by distances or by similarity
#' indices. It is only used to compute the modularity index. It can be:
#' \itemize{
#' \item{'distance': Link weights correspond to distances. Nodes that are close
#' to each other will more likely be in the same module.}
#' \item{'similarity': Link weights correspond to similarity indices. Nodes that
#' are similar to each other will more likely be in the same module. Inverse
#' link weights are then used to compute the modularity index.}
#' }
#' @param nb_modul (optional , default = NULL) A numeric or integer value
#' indicating the number of modules in the graph. When this number is not
#' specified, the optimal value is retained.
#' @return A \code{data.frame} with the node names and the corresponding
#' module ID.
#' @export
#' @author P. Savary
#' @examples
#' data("data_tuto")
#' mat_gen <- data_tuto[[1]]
#' graph <- gen_graph_thr(mat_w = mat_gen, mat_thr = mat_gen,
#'                             thr = 0.8)
#' res_mod <- compute_graph_modul(graph = graph,
#'                                 algo = "fast_greedy",
#'                                 node_inter = "distance")


compute_graph_modul <- function(graph,
                                algo = "fast_greedy",
                                node_inter = NULL,
                                nb_modul = NULL){

  # Check whether graph is a graph
  if(!inherits(graph, "igraph")){
    stop("'graph' must be a graph object of class 'igraph'.")
  }

  # Check whether the graph nodes have names
  if(is.null(igraph::V(graph)$name)){
    stop("The nodes of 'graph' must have names.")
  }

  # Check algo argument
  if(!inherits(algo, "character")){
    stop("'algo' must be a character string.")
  } else if(!(algo %in% c("fast_greedy", "walktrap", "louvain", "optimal"))){
    stop("'algo' must be either 'fast_greedy', 'walktrap',
         'louvain' or 'optimal.")
  }

  # Check node_inter argument
  if(!is.null(node_inter)){
    if(!inherits(node_inter, "character")){
      stop("When specified, 'node_inter' must be a character string.")
    } else if(!(node_inter %in% c("distance", "similarity"))){
      stop("'node_inter' must be either 'distance' or 'similarity'.")
    }
  }

  # Check nb_modul
  if(!is.null(nb_modul)){
    if(!inherits(nb_modul, c("character", "integer"))){
      stop("When specified, 'nb_modul' must a numeric or integer value.")
    }
  }

  # Check whether the graphs' links are weighted to compute some metrics
  if(!is.null(node_inter)){
    if(is.null(igraph::E(graph)$weight)){
      stop(paste0("'graph' must have weighted links in order to compute the ",
                  "specified metric with 'node_inter = ", node_inter, "' option."))
    }
  }


  ############## Compute the modules ################################

  df_mod <- data.frame(ID = igraph::V(graph)$name)

  # Link weights
  if (is.null(node_inter)){

    w <- rep(1, length(igraph::E(graph)))

  } else if(node_inter == "distance"){

    w <- 1/igraph::E(graph)$weight

  } else if (node_inter == "similarity"){

    w <- igraph::E(graph)$weight

  }

  # Creation of the modules
  if (algo == "fast_greedy"){

    # If the number of modules is not defined yet,
    # it is the number of modules
    # created by default by the algorithm.
    if(is.null(nb_modul)){
      x <- igraph::cluster_fast_greedy(graph, weights = w)$membership
    } else {
      # We specify the number of modules
      x <- igraph::cut_at(igraph::cluster_fast_greedy(graph, weights = w),
                          no = nb_modul)
    }
  } else if (algo == "louvain"){
    x <- igraph::cluster_louvain(graph, weights = w)$membership

    if(!is.null(nb_modul)){
      message("With this algorithm, 'nb_modul' parameter was not used
            and the number of modules is the default number as computed
            by the algorithm.")
    }

  } else if (algo == "optimal"){
    x <- igraph::cluster_optimal(graph, weights = w)$membership

    if(!is.null(nb_modul)){
    message("With this algorithm, 'nb_modul' parameter was not used
            and the number of modules is the default number as computed
            by the algorithm.")
    }

  } else if (algo == "walktrap"){
    if(is.null(nb_modul)){
      x <- igraph::cluster_walktrap(graph, weights = w)$membership
    } else {
      x <- igraph::cut_at(igraph::cluster_walktrap(graph, weights = w),
                          no = nb_modul)
    }
  }

  # Create a data.frame with the node partition into modules and their ID
  df_mod$module <- as.character(as.factor(as.vector(x)))

  return(df_mod)
}









