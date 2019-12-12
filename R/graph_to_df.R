#' Convert a graph into a edge list data.frame
#'
#' @description The function converts a graph into a edge list data.frame
#'
#' @param graph A graph object of class \code{igraph}
#' @param weight Logical. If TRUE (default), then the column 'link' of the
#' output data.frame contains the weights of the links. If FALSE,
#' it contains only 0 and 1.
#' @return An object of class \code{data.frame} with a link ID, the origin nodes
#' ('from') and arrival nodes ('to') and the link
#' value ('link')(weighted or binary)
#' @details The 'graph' nodes must have names. Links must have weights if
#' 'weight = TRUE'.
#' @export
#' @author P. Savary
#' @examples
#' data(pts_pop_simul)
#' suppressWarnings(mat_geo <- mat_geo_dist(pts_pop_simul,
#'                  ID = "ID",
#'                  x = "x",
#'                 y = "y"))
#' g1 <- gen_graph_thr(mat_w = mat_geo,
#'                     mat_thr = mat_geo,
#'                     thr = 20000)
#' g1_df <- graph_to_df(g1,
#'                      weight = TRUE)

graph_to_df <- function(graph, weight = TRUE){

  # Check whether graph is a graph
  if(!inherits(graph, "igraph")){
    stop("'graph' must be a graph object of class 'igraph'.")
  }

  # Check whether the graph's nodes have names
  if(is.null(igraph::V(graph)$name)){
    stop("The nodes of 'graph' must have names.")
  }

  # Check whether the graph's links have weights if weight = TRUE
  if(weight){
    if(is.null(igraph::E(graph)$weight)){
      stop("The links of 'graph' must have weights if weight = TRUE.")
    }
  }

  # Get the weighted adjacency matrix of the graph
  mat <- igraph::as_adjacency_matrix(graph, type = "both", sparse = FALSE,
                                     attr = "weight", names = TRUE)
  # Assign a 1 value to every link in the pw matrix
  mat[mat >= 0] <- 1

  # Create a complete
  g <- gen_graph_topo(mat_w = mat,
                      mat_topo = mat,
                      topo = "comp")

  # Get the edge-list of the complete graph
  df_graph <- data.frame(igraph::as_edgelist(g))
  # Create a unique ID for every link
  df_graph$ID <- paste0(df_graph$X1, "_", df_graph$X2)

  # Get the edge-list of the input graph
  df_2 <- data.frame(igraph::as_edgelist(graph))
  # Assign the weight to every link
  if(weight){
    df_2$link <- igraph::E(graph)$weight
  } else {
    df_2$link <- 1
  }
  # Give a name to every link
  df_2$ID <- paste0(df_2$X1, "_", df_2$X2)

  # Merge both edge-lists
  df_graph <- merge(df_graph, df_2[, 3:4], by = 'ID', all.x = TRUE)
  #
  df_graph$link <- ifelse(is.na(df_graph$link), 0, df_graph$link)

  names(df_graph) <- c("ID", "from", "to", "link")

  return(df_graph)
}



