#' Compute graph-theoretic metrics from a graph at the node level
#'
#' @description The function computes graph-theoretic metric values at the
#' node level.
#' @param graph An object of class \code{igraph}. Its nodes must have names.
#' @param metrics Character vector specifying the graph-theoretic
#' metrics computed at the node-level in the graphs
#' Graph-theoretic metrics can be:
#' \itemize{
#' \item{Degree (\code{metrics = c("deg", ...)})}
#' \item{Closeness centrality index (\code{metrics = c("close",...)})}
#' \item{Betweenness centrality index (\code{metrics = c("btw",...)})}
#' \item{Strength (sum of the weights of the links connected to a node)
#' (\code{metrics = c("str",...)})}
#' \item{Sum of the inverse weights of the links connected to a
#' node (\code{metrics = c("siw", ...)}, default)}
#' \item{Mean of the inverse weights of the links connected to a
#' node (\code{metrics = c("miw", ...)})}
#' }
#' By default, the vector \code{metrics} includes all these metrics.
#' @param weight Logical which indicates whether the links are weighted during
#' the calculation of the centrality indices betweenness and closeness.
#' (default: \code{weight = TRUE}). Link weights are interpreted as distances
#' when computing the shortest paths. They should then be inversely proportional
#' to the strength of the relationship between nodes (e.g. to fluxes).
#' @return A \code{data.frame} with the node names and the metrics computed.
#' @export
#' @author P. Savary
#' @examples
#' data(data_ex_genind)
#' mat_gen <- mat_gen_dist(x = data_ex_genind, dist = "DPS")
#' graph <- gen_graph_thr(mat_w = mat_gen, mat_thr = mat_gen,
#'                             thr = 0.8)
#' res_met <- compute_node_metric(graph)


compute_node_metric <- function(graph,
                                metrics = c("deg", "close", "btw",
                                            "str", "siw", "miw"),
                                weight = TRUE){


  # Check whether graph is a graph
  if(!inherits(graph, "igraph")){
    stop("'graph' must be a graph object of class 'igraph'.")
  }

  # Check whether the graph nodes have names
  if(is.null(igraph::V(graph)$name)){
    stop("The nodes of 'graph' must have names.")
  }

  # Check whether weight is logical
  if(!is.logical(weight)){
    stop("'weight' must be TRUE or FALSE.")
  }

  # Check whether the vector 'metrics' is a character vector.
  if(!inherits(metrics, "character")){
    stop("'metrics' vector must be a character vector")
  }

  # Vector of metrics options available
  metrics_options <- c("deg", "close", "btw", "siw", "miw", "str")

  # Check whether 'metrics' options are valid
  if(!all(metrics %in% metrics_options)){
    stop(paste0("You must specify valid metric names among the values:",
         paste(metrics_options, collapse = ", ")))
  }

  # Check whether the graphs' links are weighted to compute some metrics
  if(weight){
    if(is.null(igraph::E(graph)$weight)){
      stop("'graph' must have weighted links in order to compute the
           specified metric with 'weight = TRUE' option.")
    }
  }


  ############## Compute the values ################################

  df_values <- data.frame(ID = igraph::V(graph)$name)

  if("deg" %in% metrics){
      # Degree
    df_values$deg <- igraph::degree(graph)
    }

  if("close" %in% metrics){
      # Closeness with or without weights
      if(weight){
        df_values$close <- igraph::closeness(graph)
      } else if (weight == FALSE) {
        df_values$close <- igraph::closeness(graph,
                                             weights = rep(1,
                                                        length(igraph::E(graph))))
      }
    }

  if("btw" %in% metrics){
      # Betweenness with or without weights
      if(weight){
        df_values$btw <- igraph::betweenness(graph)
      } else if (weight == FALSE) {
        df_values$btw <- igraph::betweenness(graph,
                                             weights = rep(1,
                                                          length(igraph::E(graph))))
      }
    }

  if("siw" %in% metrics){
    if(weight){
      df_values$siw <- igraph::strength(graph,
                                        weights = 1/igraph::E(graph)$weight)
    } else {
      df_values$siw <- igraph::strength(graph,
                                        weights = rep(1,
                                                     length(igraph::E(graph))))
    }
  }


  if("miw" %in% metrics){
    if(weight){
      df_values$miw <- igraph::strength(graph,
                                    weights = 1/igraph::E(graph)$weight)/igraph::degree(graph)
    } else {
      df_values$miw <- igraph::strength(graph,
                                        weights = rep(1,
                                length(igraph::E(graph))))/igraph::degree(graph)
    }
  }


  if("str" %in% metrics){
    if(weight){
      df_values$str <- igraph::strength(graph)
    } else {
      df_values$str <- igraph::strength(graph, weights = rep(1,
                                              length(igraph::E(graph))))
    }
  }

  ### Create a data.frame with the values

  df_values <- df_values[, c("ID", metrics)]

  return(df_values)
}





