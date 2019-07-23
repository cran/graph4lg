#' Plot graphs
#'
#' @description The function enables to plot graphs, whether spatial or not.
#'
#' @param graph A graph object of class \code{igraph}
#' @param crds (if 'mode = 'spatial'') A \code{data.frame} with the spatial
#' coordinates of the graph's nodes. It must have three columns :
#' \itemize{
#' \item{ID: A character string indicating the name of the graph's nodes.
#' The names must be the same as the nodes' names of the graph of
#' class \code{igraph} (\code{igraph::V(graph)$name})}
#' \item{x: A character string indicating the longitude of the graphs' nodes.}
#' \item{y: A character string indicating the latitude of the graphs' nodes.}
#' }
#' @param mode A character string indicating whether the graph is
#' spatial ('mode = 'spatial'' (default)) or not ('mode = 'aspatial'')
#' @param weight A character string indicating whether the links of the graph have
#' weights (T)(default) or not (F)
#' @param width A character string indicating whether the width of the link
#' should be proportional to links' weights ("w", default) or to the inverse
#' of links' weights ("inv", convenient with distances)
#' @param pts_col (optional) A character string indicating the color
#' used to plot the nodes (default: "#F2B950"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @details When the graph is not spatial ('mode = 'aspatial''), the nodes coordinates
#' are calculated with Fruchterman et Reingold algorithm.
#' The graph object \code{x} of class \code{igraph} must have nodes' names
#' (not necessarily in the same order as IDs in crds, given a merging is done).
#' @author P. Savary
#' @references \insertRef{fruchterman1991graph}{graph4lg}
#' @examples
#' data(pts_pop_simul)
#' data(data_simul_genind)
#' mat_w <- mat_gen_dist(data_simul_genind, dist = "DPS")
#' gp <- gen_graph_topo(mat_w = mat_w, topo = "mst")
#' g <- plot_graph_lg(graph = gp,
#'                              crds = pts_pop_simul,
#'                              mode = "spatial",
#'                              weight = TRUE)



plot_graph_lg <- function(graph,
                          crds,
                          mode = "spatial",
                          weight = TRUE,
                          width = "w",
                          pts_col = "#F2B950"){


  # Check whether the graph has nodes' names
  if(is.null(igraph::V(graph)$name)){
    stop("Your graph must have nodes' names.")
  }

  # If the graph's nodes have spatial coordinates
  if(mode == "spatial"){

    # If 'crds' is given, check whether it has as many elements as there are nodes
    if(!exists("crds")){
      stop("You must provide the spatial coordinates of the graph's nodes.")
    } else if( nrow(crds) != length(igraph::V(graph) ) ) {
      stop("'crds' must have the same number of rows as there are nodes in 'graph'")
    # Check whether 'crds' has valid columns' names
    } else if( !all( colnames(crds) == c("ID","x","y") ) ){
      stop("Column names of crds must be 'ID', 'x' and 'y'.")
    # Check whether IDs from 'crds' match with the graph's nodes names
    } else if( !any( as.character(crds$ID) %in% igraph::V(graph)$name ) ){
      stop("The IDs of 'crds' elements are not the same as the names of the nodes in 'graph'")
    # Check whether spatial coordinates are numeric
    } else if(class(crds$x) != "numeric"){
      stop("'x' must be of class 'numeric'")
    } else if(class(crds$y) != "numeric"){
      stop("'y' must be of class 'numeric'")
    } else {
      crds$ID <- as.character(crds$ID)
    }

    # Create a data.frame from the links of the graph
    graph_df <- data.frame(igraph::as_edgelist(graph))

    # Add weights or not to the links
    if (weight == TRUE){
      graph_df$w <- igraph::E(graph)$weight
      names(graph_df) <- c("from", "to", "w")
    } else {
      names(graph_df) <- c("from", "to")
    }

    # Merge 'graph_df' and 'crds' in order to get spatial coordinates of the nodes
    graph_df <- merge(graph_df, crds, by.x = 1, by.y = 1)
    graph_df <- merge(graph_df, crds, by.x = 2, by.y = 1)

    # If links are weighted
    if (weight == TRUE){
      # Give colnames to 'graph_df'
      names(graph_df) <- c("to", "from", "weight", "x", "y", "xend", "yend")
      # Create a scaling function
      sc01 <- function(x){(x-min(x))/(max(x)-min(x))}

      if (width == "w"){
      # Compute the scaled weight
      graph_df$w_sc <- sc01(graph_df$weight)
      } else if (width == "inv"){
      # Compute the inverse weight and scale it
      graph_df$w_sc <- 1/sc01(graph_df$weight)
      # Replace 'Inf' values by the largest values
      graph_df[which(graph_df$w_sc == "Inf" ), 'w_sc'] <- max(graph_df[which(graph_df$w_sc < Inf), "w_sc"]) + 2
      # Rescale again
      graph_df$w_sc <- sc01(graph_df$w_sc)
      } else {
        stop("You must specify a correct 'width' option.")
      }

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                          xend = 'xend', yend = 'yend',
                                          size = 'w_sc'), color = "black")+
        scale_size_identity()+
        geom_point(data = crds, aes_string(x = 'x', y = 'y'),
                 size = 6, color = pts_col) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="Longitude", y="Latitude")
        #theme(panel.background = element_rect(fill = 'black', colour = 'white'))

    # If the links are not weighted
    } else {
      # Give columns' names to 'graph_df'
      names(graph_df) <- c("to", "from", "x", "y", "xend", "yend")

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                          xend = 'xend', yend = 'yend'),
                                          color = "black")+
        scale_size_identity()+
        geom_point(data = crds, aes_string(x = 'x', y = 'y'),
                   size = 6, color = pts_col) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="Longitude", y="Latitude")
        #theme(panel.background = element_rect(fill = 'black', colour = 'white'))
    }


  # If the graph's nodes do not have spatial coordinates
  } else if (mode == "aspatial"){

    # If links are weighted
    if (weight == TRUE){

    # Use the Fruchterman and Reingold algorithm to compute the nodes' coordinates
    # giving a weight to the links equal to the inverse of the distance they
    # are weighted with in the graph
    crds <- igraph::layout_with_fr(graph, weights = 1/(igraph::E(graph)$weight))
    } else {
    # Use the Fruchterman and Reingold algorithm to compute the nodes' coordinates
    # giving a weight of 1 to the links
    crds <- igraph::layout_with_fr(graph, weights = rep(1, length(igraph::E(graph))))
    }
    # Create the data.frame 'crds'
    crds <- data.frame(crds)
    names(crds) <- c("x", "y")

    # IDs of 'crds' elements are graph's nodes names
    crds$ID <- igraph::V(graph)$name
    crds <- crds[, c('ID', 'x', 'y')]

    # Create a data.frame from the links of the graph
    graph_df <- data.frame(igraph::as_edgelist(graph))

    # Add weights or not to the links
    if (weight == TRUE){
      graph_df$w <- igraph::E(graph)$weight
      names(graph_df) <- c("from", "to", "w")
    } else {
      names(graph_df) <- c("from", "to")
    }

    # Merge 'graph_df' and 'crds' in order to get spatial coordinates of the nodes
    graph_df <- merge(graph_df, crds, by.x = 1, by.y = 1)
    graph_df <- merge(graph_df, crds, by.x = 2, by.y = 1)

    # If links are weighted
    if (weight == TRUE){
      # Give colnames to 'graph_df'
      names(graph_df) <- c("to", "from", "weight", "x", "y", "xend", "yend")
      # Create a scaling function
      sc01 <- function(x){(x-min(x))/(max(x)-min(x))}


      if (width == "w"){
        # Compute the scaled weight
        graph_df$w_sc <- sc01(graph_df$weight)
      # Compute the inverse weight and scale it
      } else if (width == "inv"){
        graph_df$w_sc <- 1/sc01(graph_df$weight)
        # Replace 'Inf' values by the largest values
        graph_df[which(graph_df$w_sc == "Inf" ), 'w_sc'] <- max(graph_df[which(graph_df$w_sc < Inf), "w_sc"]) + 2
        # Rescale again
        graph_df$w_sc <- sc01(graph_df$w_sc)
      } else {
        stop("You must specify a correct 'width' option.")
      }

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                          xend = 'xend', yend = 'yend',
                                          size = 'w_sc'), color = "black")+
        scale_size_identity()+
        geom_point(data = crds, aes_string(x = 'x', y = 'y'),
                   size = 6, color = pts_col) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="x", y="y")
        #theme(panel.background = element_rect(fill = 'black', colour = 'white'))

    # If the links are not weighted
    } else {
      # Give columns' names to 'graph_df'
      names(graph_df) <- c("to", "from", "x", "y", "xend", "yend")

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                          xend = 'xend', yend = 'yend'),
                     color = "black")+
        scale_size_identity()+
        geom_point(data = crds, aes_string(x = 'x', y = 'y'),
                   size = 6, color = pts_col) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="x", y="y")
        #theme(panel.background = element_rect(fill = 'black', colour = 'white'))
    }

  } else {
    stop("You must specify a correct 'mode' option ('spatial' or 'aspatial')")
  }
 return(g)

}




