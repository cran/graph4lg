#' Compare the local properties of the nodes from two graphs
#'
#' @description The function computes a correlation coefficient between the
#' graph-theoretic metric values computed at the node-level in two graphs
#' sharing the same nodes. It allows to assess whether the connectivity
#' properties of the nodes in one graph are similar to that of the same nodes
#' in the other graph. Alternatively, the correlation is computed between
#' a graph-theoretic metric values and the values of an attribute associated
#' to the nodes of a graph.
#'
#' @details The correlation coefficients between the metrics can be computed
#' in different ways, as initial assumptions (e.g. linear relationship) are
#' rarely verified. Pearson's r, Spearman's rho and Kendall's tau can be
#' computed (from function \code{\link[stats]{cor}}).
#' When \code{x} is similar to \code{y}, then the correlation is computed
#' between two metrics characterizing the nodes of the same graph.
#' @param x An object of class \code{igraph}.
#' Its nodes must have the same names as in graph \code{y}.
#' @param y An object of class \code{igraph}.
#' Its nodes must have the same names as in graph \code{x}.
#' @param metrics Two-elements character vector specifying the graph-theoretic
#' metrics computed at the node-level in the graphs or the nodes' attribute
#' values to be correlated to these metrics.
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
#' Nodes' attributes must have the same names as in the \code{igraph} object,
#' and must refer to an attribute with numerical values.
#' The vector \code{metrics} is composed of two character values.
#' When a nodes' attribute has the same name as a metric computable from the
#' graph, nodes' attributes are given priority.
#' @param method A character string indicating which correlation coefficient
#' is to be computed (\code{"pearson"}, \code{"kendall"} or
#' \code{"spearman"} (default)).
#' @param weight Logical which indicates whether the links are weighted during
#' the calculation of the centrality indices betweenness and closeness.
#' (default: \code{weight = TRUE}). Links' weights are interpreted as distances
#' when computing the shortest paths. They should then be inversely proportional
#' to the strength of the relationship between nodes (e.g. to fluxes).
#' @param test Logical. Should significance testing be performed?
#' (default = TRUE)
#' @return A \code{list} summarizing the correlation analysis.
#' @export
#' @author P. Savary
#' @examples
#' data(data_ex_genind)
#' data(pts_pop_ex)
#' mat_dist <- suppressWarnings(graph4lg::mat_geo_dist(data = pts_pop_ex,
#'       ID = "ID",
#'       x = "x",
#'       y = "y"))
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                       order(as.character(colnames(mat_dist)))]
#' graph_obs <- gen_graph_thr(mat_w = mat_dist, mat_thr = mat_dist,
#'                            thr = 9500, mode = "larger")
#' mat_gen <- mat_gen_dist(x = data_ex_genind, dist = "DPS")
#' graph_pred <- gen_graph_topo(mat_w = mat_gen, mat_topo = mat_dist,
#'                             topo = "gabriel")
#' res_cor <- graph_node_compar(x = graph_obs, y = graph_pred,
#'                              metrics = c("siw", "siw"), method = "spearman",
#'                              test = TRUE, weight = TRUE)


graph_node_compar <- function(x,
                              y,
                              metrics = c("siw", "siw"),
                              method = "spearman",
                              weight = TRUE,
                              test = TRUE){

  # Check whether x and y are graphs
  if(!inherits(x, "igraph")){
    stop("'x' must be a graph object of class 'igraph'.")
  } else if (!inherits(y, "igraph")){
    stop("'y' must be a graph object of class 'igraph'.")
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
    stop("Both graphs must have the same nodes' names and the nodes ranked
         in the same order.")
  }

  # Check whether the vector 'metrics' is a two-element character vector.
  if(length(metrics)!= 2){
    stop("'metrics' vector must be of a length 2")
  } else if(!inherits(metrics, "character")) {
    stop("'metrics' vector must be a character vector")
  }

  # Vector of metrics options available
  metrics_options <- c("deg", "close", "btw", "siw", "miw", "str")

  # Names of the nodes' attributes of x and y
  attrib_x <- names(igraph::vertex.attributes(x))
  attrib_y <- names(igraph::vertex.attributes(y))

  # Check whether 'metrics' options are valid
  if(!all(metrics %in% c(metrics_options, attrib_x, attrib_y))){
    stop("You must specify valid metrics' names, either among the values
           'deg', 'close', 'btw', 'siw', 'miw' and 'str', or from the nodes'
           attributes names.")
  }

  # Check whether the graphs' links are weighted to compute some metrics
  if(metrics[1] %in% c("close", "btw", "siw", "miw", "str")){
    if(is.null(igraph::E(x)$weight)){
      stop("x must have weighted links in order to compute the
           specified metric.")
    }
  }

  if(metrics[2] %in% c("close", "btw", "siw", "miw", "str")){
    if(is.null(igraph::E(y)$weight)){
      stop("y must have weighted links in order to compute the
           specified metric.")
    }
  }

  # Check whether the metric considered from graph x is a nodes' attribute
  # or is to be computed.
  # Nodes attribute are given priority
  if(metrics[1] %in% attrib_x){
    met_from_x <- "attrib"
    num_met_x <- which(attrib_x == metrics[1])

  } else {
    met_from_x <- "comput"
  }

  # Check whether the metric considered from graph y is a nodes' attribute
  # or is to be computed.
  # Nodes attribute are given priority
  if(metrics[2] %in% attrib_y){
    met_from_y <- "attrib"
    num_met_y <- which(attrib_y == metrics[2])

  } else {
    met_from_y <- "comput"
  }

  # Check whether 'method' option is a valid one
  if(!(method %in% c("pearson", "spearman", "kendall"))){
    stop("You must specify a valid method to compute the
         correlation coefficient")
  }

  ############## Get the values to correlate ################################

  ### Graph x

  if(met_from_x == "attrib"){
    # Get the nodes' attribute from the attributes of the igraph object
    met_val_x <- igraph::vertex.attributes(x)[[num_met_x]]
  # Check whether met_val_x is of class numeric or integer
    if(!inherits(met_val_x, c("numeric", "integer"))){
      stop("Nodes' attributes must be of class 'numeric' or 'integer'")
    }

  # Compute the nodes' attributes
  } else if (met_from_x == "comput"){
    if (metrics[1] == "deg"){
      # Degree
      met_val_x <- igraph::degree(x)
    } else if (metrics[1] == "close"){
      # Closeness with or without weights
      if(weight){
        met_val_x <- igraph::closeness(x)
      } else if (weight == FALSE) {
        met_val_x <- igraph::closeness(x, weights = rep(1,
                                                        length(igraph::E(x))))
      } else {
        stop("'weight' must be TRUE or FALSE.")
      }
    } else if (metrics[1] == "btw"){
      # Betweenness with or without weights
      if(weight){
        met_val_x <- igraph::betweenness(x)
      } else if (weight == FALSE) {
        met_val_x <- igraph::betweenness(x, weights = rep(1,
                                                          length(igraph::E(x))))
      } else {
        stop("'weight' must be TRUE or FALSE.")
      }
    } else if (metrics[1] == "siw"){
      met_val_x <- igraph::strength(x,
                                    weights = 1/igraph::E(x)$weight)
    } else if (metrics[1] == "miw"){
      met_val_x <- igraph::strength(x,
                            weights = 1/igraph::E(x)$weight)/igraph::degree(x)
    } else if (metrics[1] == "str"){
      met_val_x <- igraph::strength(x)
    }
  }


  ### Graph y

  # Same steps for the graph y

  if(met_from_y == "attrib"){
    met_val_y <- igraph::vertex.attributes(y)[[num_met_y]]

    if(!inherits(met_val_y, c("numeric", "integer"))){
      stop("Nodes' attributes must be of class 'numeric' or 'integer'")
    }

  } else if (met_from_y == "comput"){
    if (metrics[2] == "deg"){
      met_val_y <- igraph::degree(y)
    } else if (metrics[2] == "close"){
      if(weight){
        met_val_y <- igraph::closeness(y)
      } else if (weight == FALSE) {
        met_val_y <- igraph::closeness(y, weights = rep(1,
                                                        length(igraph::E(y))))
      } else {
        stop("'weight' must be TRUE or FALSE.")
      }
    } else if (metrics[2] == "btw"){
      if(weight){
        met_val_y <- igraph::betweenness(y)
      } else if (weight == FALSE) {
        met_val_y <- igraph::betweenness(y, weights = rep(1,
                                                          length(igraph::E(y))))
      } else {
        stop("'weight' must be TRUE or FALSE.")
      }
    } else if (metrics[2] == "siw"){
      met_val_y <- igraph::strength(y,
                                    weights = 1/igraph::E(y)$weight)
    } else if (metrics[2] == "miw"){
      met_val_y <- igraph::strength(y,
                            weights = 1/igraph::E(y)$weight)/igraph::degree(y)
    } else if (metrics[2] == "str"){
      met_val_y <- igraph::strength(y)
    }
  }

  ### Create a data.frame with the values

  df_values <- data.frame(name = igraph::V(x)$name,
                          met_x = met_val_x,
                          met_y = met_val_y)

  ### Compute a correlation coefficient

  if (test){
    res_r <- stats::cor.test(df_values$met_x, df_values$met_y,
                             method = method, exact = FALSE)

    pval <- res_r$p.value
    r <- res_r$estimate

    res_list <- list(paste0("Metric from graph x: ", metrics[1]),
                     paste0("Metric from graph y: ", metrics[2]),
                     paste0("Method used: ", method,
                            "'s correlation coefficient"),
                     paste0("Sample size: ", n_nodes),
                     paste0("Correlation coefficient: ", r),
                     paste0("p-value of the significance test: ", pval))


  } else if (test == FALSE){
    res_r <- stats::cor(df_values$met_x, df_values$met_y, method = method)
    res_list <- list(paste0("Metric from graph x: ", metrics[1]),
                     paste0("Metric from graph y: ", metrics[2]),
                     paste0("Method used: ", method,
                            "'s correlation coefficient"),
                     paste0("Sample size: ", n_nodes),
                     paste0("Correlation coefficient: ", res_r))

  }

  return(res_list)
}





