#' Visualize the topological differences between two spatial graphs on a map
#'
#' @description The function enables to compare two spatial graphs by
#' plotting them highlighting the topological similarities and differences
#' between them. Both graphs should share the same nodes and cannot
#' be directed graphs.
#'
#' @param x A graph object of class \code{igraph}.
#' Its nodes must have the same names as in graph \code{y}.
#' @param y A graph object of class \code{igraph}.
#' Its nodes must have the same names as in graph \code{x}.
#' @param crds A \code{data.frame} with the spatial
#' coordinates of the graph nodes (both \code{x} and \code{y}).
#' It must have three columns:
#' \itemize{
#' \item{ID: Name of the graph nodes (character string).
#' The names must be the same as the node names of the graphs of
#' class \code{igraph} (\code{igraph::V(graph)$name})}
#' \item{x: Longitude of the graph nodes (numeric or integer).}
#' \item{y: Latitude of the graph nodes (numeric or integer).}
#' }
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @details The graphs \code{x} and \code{y} of class \code{igraph} must have
#' node names (not necessarily in the same order as IDs in crds,
#' given a merging is done).
#' @author P. Savary
#' @examples
#' data(pts_pop_ex)
#' data(data_ex_genind)
#' mat_w <- mat_gen_dist(data_ex_genind, dist = "DPS")
#' mat_dist <- mat_geo_dist(data = pts_pop_ex,
#'                          ID = "ID",
#'                          x = "x",
#'                          y = "y")
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                    order(as.character(colnames(mat_dist)))]
#' g1 <- gen_graph_topo(mat_w = mat_w, topo = "mst")
#' g2 <- gen_graph_topo(mat_w = mat_w, mat_topo = mat_dist, topo = "gabriel")
#' g <- graph_plot_compar(x = g1, y = g2,
#'                        crds = pts_pop_ex)



graph_plot_compar <- function(x, y,
                              crds){


  #######################################################
  ####################################################

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
    stop("Both graphs must have the same nodes' names and the
         nodes ranked in the same order.")
  }

  if(!exists("crds")){
    stop("You must provide the spatial coordinates of the graph nodes.")
  } else if( nrow(crds) != length(igraph::V(x) ) ) {
    stop("'crds' must have as many rows as there are nodes in 'x' and 'y'.")
  } else if( !all( colnames(crds) == c("ID","x","y") ) ){
    stop("Column names of crds must be 'ID', 'x' and 'y'.")
  } else if( !any( as.character(crds$ID) %in% igraph::V(x)$name ) ){
    stop("The IDs of 'crds' elements are not the same as the names of
         the nodes in 'x' and 'y'.")
    # Check whether spatial coordinates are numeric or integer
  } else if(!inherits(crds$x, c("integer", "numeric"))){
    stop("'x' must be of class 'numeric' or 'integer'")
  } else if(!inherits(crds$y, c("integer", "numeric"))){
    stop("'y' must be of class 'numeric' or 'integer'")
  } else {
    crds$ID <- as.character(crds$ID)
  }

  # Get the adjacency matrix of the graph x
  adj_x <- igraph::as_adjacency_matrix(x, type = "both", attr = NULL,
                                       sparse = FALSE)

  # Get the adjacency matrix of the graph y
  adj_y <- igraph::as_adjacency_matrix(y, type = "both", attr = NULL,
                                       sparse = FALSE)

  # Add both adjacency matrix
  simil <- adj_x + adj_y
  # Remove the adjacency of one to that of the other
  diff <- adj_x - adj_y

  # Create an empty matrix
  mat_xy <- matrix(data = rep(0, n_nodes*n_nodes),
                   nrow = n_nodes,
                   ncol = n_nodes)
  row.names(mat_xy) <- colnames(mat_xy) <- row.names(adj_x)

  mat_xy[which(simil == 2)] <- 2 # If simil = 2, the link is both in x and y
  mat_xy[which(diff == 1)] <- 3 # If diff = 1, the link is in x only
  mat_xy[which(diff == -1)] <- 4 # If diff = -1, the link is in y only
  mat_xy[which(mat_xy == 0)] <- 1 # Else, neither in x nor in y

  # Create a graph g_xy from mat_xy
  g_xy <- igraph::graph.adjacency(mat_xy,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)

  # Create a data.frame describing the links of g_xy
  graph_df <- data.frame(igraph::as_edgelist(g_xy))

  # Add the weights
  graph_df$link <- igraph::E(g_xy)$weight
  names(graph_df) <- c("from", "to", "link")

  # Merge graph_df and crds to get spatial information about the nodes
  graph_df <- merge(graph_df, crds, by.x = 1, by.y = 1)
  graph_df <- merge(graph_df, crds, by.x = 2, by.y = 1)

  names(graph_df) <- c("to", "from", "link", "x", "y", "xend", "yend")

  graph_df$width <- graph_df$link

  # Adjust the width of the links to their weight
  graph_df[which(graph_df$link == 2), 'width'] <- 1
  graph_df[which(graph_df$link == 3), 'width'] <- 0.5
  graph_df[which(graph_df$link == 4), 'width'] <- 0.5

  graph_df$link <- as.factor(graph_df$link)

  # Colors
  palette_topo <- c("black", "#8C2A1C","#808080")

  # Create the plot to compare the graphs
  g <- ggplot() +
    geom_segment(data = graph_df[which(graph_df$link != "1"), ],
                 aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend,
                     color = .data$link, size = .data$width))+
    scale_size_identity()+
    geom_point(data = crds, aes(x = .data$x, y = .data$y),
               size = 6, color = "#999999") +
    geom_text(data = crds, aes(x = .data$x, y = .data$y, label = .data$ID),
              size = 3, color = "black", fontface = "bold")+
    theme_bw()+
    labs(x="Longitude", y="Latitude", color = "Topological comparison")+
    scale_color_manual(values = palette_topo,
                       labels = c("In both graphs", "In x only", "In y only"))+
    theme(legend.position = "bottom")
  #theme(panel.background = element_rect(fill = 'black', colour = 'white'))

  return(g)

}




