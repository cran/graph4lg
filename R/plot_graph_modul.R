#' Plot the graphs making visible their partition into modules
#'
#' @description The function computes a partition of the graph into modules
#' and plots the graph colouring the nodes with colors corresponding to
#' their respective modules.
#'
#' @param graph A graph object of class \code{igraph} to partition into modules
#' and to plot. The graph must be undirected. If 'crds' is not NULL, then the
#' graph nodes must have names corresponding to the ID column of 'crds'.
#' @param algo (if x and y and graph objects) A character string indicating the
#' algorithm used to create the modules with \pkg{igraph}:
#' \itemize{
#' \item{If \code{algo = 'fast_greedy'} (default), function
#' \code{cluster_fast_greedy} from \pkg{igraph} is
#' used (Clauset et al., 2004).}
#' \item{If \code{algo = 'walktrap'}, function \code{cluster_walktrap}
#' from \pkg{igraph} is used (Pons et Latapy, 2006) with
#' 4 steps (default options).}
#' \item{If \code{algo = 'louvain'}, function \code{cluster_louvain}
#' from \pkg{igraph} is used (Blondel et al., 2008).
#' In that case, the number of modules created in each graph is imposed.}
#' \item{If \code{algo = 'optimal'}, function \code{cluster_optimal}
#' from \pkg{igraph} is used (Brandes et al., 2008) (can be very long).
#' In that case, the number of modules created in each graph is imposed.}
#' }
#' @param nb_modul Numeric value indicating the number of modules to create.
#' @param weight (optional) A character string or character vector
#' indicating how to weight graph links during the calculation
#' of the modularity.
#' \itemize{
#' \item{If \code{weight = 'inv'} (default), then links are weighted with the
#' inverse values of their initial weights.}
#' \item{If \code{weight = 'w'}, then links are weighted with their initial
#' weights values.}
#' \item{If \code{weight = 'none'}, then links are not weighted during the
#' calculation.}
#' }
#' If the graph links are not weighted, then this argument is ignored.
#' Links with large weights are considered as stronger connections
#' in the modularity calculation.
#' @param crds (if 'mode = 'spatial'') A \code{data.frame} with the spatial
#' coordinates of the graph nodes. It must have three columns :
#' \itemize{
#' \item{ID: A character string indicating the name of the graph nodes.
#' The names must be the same as the node names of the graph of
#' class \code{igraph} (\code{igraph::V(graph)$name})}
#' \item{x: A numeric or integer indicating the longitude of the graph nodes.}
#' \item{y: A numeric or integer indicating the latitude of the graph nodes.}
#' }
#' @param mode A character string indicating whether the graph is
#' spatial ('mode = 'spatial'' (default)) or not ('mode = 'aspatial'')
#' @param weight_plot Logical indicating whether the links of the graph have
#' to be displayed on the plot
#' @param width Logical indicating whether the width of the links
#' on the plot should be proportional to link weights ("w", default) or
#' to the inverse of link weights ("inv", convenient with distances)
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @details When the graph is not spatial ('mode = 'aspatial''),
#' the nodes coordinates are calculated with Fruchterman et Reingold algorithm.
#' @author P. Savary
#' @references \insertRef{hubert1985comparing}{graph4lg}
#' \insertRef{fruchterman1991graph}{graph4lg}
#' \insertRef{clauset2004finding}{graph4lg}
#' \insertRef{blondel2008fast}{graph4lg}
#' \insertRef{brandes2008modularity}{graph4lg}
#' \insertRef{pons2006computing}{graph4lg}
#' @examples
#' data(data_simul_genind)
#' data(pts_pop_simul)
#' mat_dist <- suppressWarnings(graph4lg::mat_geo_dist(data=pts_pop_simul,
#'       ID = "ID",
#'       x = "x",
#'       y = "y"))
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                       order(as.character(colnames(mat_dist)))]
#' graph <- gen_graph_thr(mat_w = mat_dist, mat_thr = mat_dist,
#'                             thr = 9500, mode = "larger")
#' plot_graph_modul(graph = graph, crds = pts_pop_simul)

plot_graph_modul <- function(graph,
                             algo = "fast_greedy",
                             nb_modul = NULL,
                             weight = "inv",
                             crds,
                             mode = "spatial",
                             weight_plot = TRUE,
                             width = "inv"){

  # Create a scaling function
  sc01 <- function(x){(x-min(x))/(max(x)-min(x))}

  # Custom palette with 26 colors
  palette <- c("#96A725","#1D72F5","#DF0101",
               "#FF9326","#A945FF","#0089B2",
               "#FDF060","#FFA6B2","#BFF217",
               "#60D5FD","#CC1577","#F2B950",
               "#7FB21D","#EC496F","#326397",
               "#B26314","#027368","#A4A4A4",
               "#610B5E", "red", "green", "darkblue",
               "orange", "blue", "grey", "black")

  # Check whether 'graph' is a graph
  if(!inherits(graph, "igraph")){
    stop("'graph' must be an object of class 'igraph'.")
  }

  # Check whether 'graph' is undirected.
  if(igraph::is.directed(graph)){
    stop("This function works for undirected graphs only")
  }

  # Store the number of nodes
  n_nodes <- length(igraph::V(graph))

  # If 'mode = 'spatial'', check whether the IDs from 'crds' corresponds to the
  # nodes names.
  if(mode == "spatial"){

    # Check whether the graph has node names
    if(is.null(igraph::V(graph)$name)){
      stop("Your graph must have node names.")
    }

    # If 'crds' is given, check whether it has as many elements
    # as there are nodes
    if(!exists("crds")){
      stop("You must provide the spatial coordinates of the graph nodes.")
    } else if( nrow(crds) != length(igraph::V(graph) ) ) {
      stop("'crds' must have the same number of rows as there
           are nodes in 'graph'")
      # Check whether 'crds' has valid column names
    } else if( !all( colnames(crds) == c("ID","x","y") ) ){
      stop("Column names of crds must be 'ID', 'x' and 'y'.")
      # Check whether IDs from 'crds' match with the graph nodes names
    } else if( !any( as.character(crds$ID) %in% igraph::V(graph)$name ) ){
      stop("The IDs of 'crds' elements are not the same as the names
           of the nodes in 'graph'")
      # Check whether spatial coordinates are numeric or integer
    } else if(!inherits(crds$x, c("integer", "numeric"))){
      stop("'x' must be of class 'numeric' or 'integer'")
    } else if(!inherits(crds$y, c("integer", "numeric"))){
      stop("'y' must be of class 'numeric' or 'integer'")
    } else {
      crds$ID <- as.character(crds$ID)
    }

    # If mode = aspatial, name the nodes with a character string
    # corresponding to their row number.
  } else if (mode == "aspatial") {
    igraph::V(graph)$name <- as.character(1:n_nodes)
  }  else {
    stop("'mode' must be 'spatial' or 'aspatial' (character strings).")
  }

  # Weights of the links during the calculation
  # If links are weighted
  if (igraph::is.weighted(graph)){
    # Inverse weight
    if(weight == "inv"){
      w <- 1/igraph::E(graph)$weight
      # Initial weight
    } else if (weight == "w"){
      w <- igraph::E(graph)$weight
      # Links are given 1 values as weights
    } else if (weight == "none"){
      w <- rep(1, length(igraph::E(graph)))
    } else {
      stop("'weight' must be either 'inv', 'w' or 'none'.")
    }
    # If links do not have weights, then links are given 1 values as weights
  } else {
    w <- rep(1, length(igraph::E(graph)))
    message("graph is not a weighted graph and its links were given 1 values
            as weights in the calculations.")
  }

  # Number of modules to create in both graphs
  if (inherits(nb_modul, "numeric")){
    n_m <- nb_modul
  } else if (is.null(nb_modul)){
    n_m <- NULL
    ##########################################################
    # Number of modules will be determined later
  } else {
    stop("'nb_modul' must be NULL, a numeric value or a numeric
         vector of length 2.")
  }

  # Creation of the modules
  if (algo == "fast_greedy"){
    # If the number of modules is not defined yet,
    # it is the number of modules
    # created by default by the algorithm.
    if(is.null(n_m)){
      n_m <- length(unique(igraph::cluster_fast_greedy(graph,
                                                       weights = w)$membership))
    }
    # We create the modules with the right algorithm, the right weighting
    # and we create the specified number of modules in the graph.
    x <- igraph::cut_at(igraph::cluster_fast_greedy(graph, weights = w),
                        no = n_m)

  } else if (algo == "louvain"){
    x <- igraph::cluster_louvain(graph, weights = w)$membership
    message("With this algorithm, 'nb_modul' parameter was not used
            and the number of modules is the default number
            as computed by the algorithm.")

  } else if (algo == "optimal"){
    x <- igraph::cluster_optimal(graph, weights = w)$membership
    message("With this algorithm, 'nb_modul' parameter was not used
            and the number of modules is the default number as computed
            by the algorithm.")

  } else if (algo == "walktrap"){
    if(is.null(n_m)){
      n_m <- length(unique(igraph::cluster_walktrap(graph,
                                                    weights = w)$membership))
    }
    x <- igraph::cut_at(igraph::cluster_walktrap(graph, weights = w), no = n_m)

  } else {
    stop("You must specify a correct 'algo' option.")
  }

  # Create a data.frame with the node partition into modules and their ID
  modul_df <- data.frame(ID = igraph::V(graph)$name,
                         modul = as.factor(as.vector(x)))
  modul_df$ID <- as.character(modul_df$ID)

  # If the graph nodes have spatial coordinates
  if (mode == "spatial") {

    # Create a df with the edge list
    graph_df <- data.frame(igraph::as_edgelist(graph))
    graph_df[, 1:2] <- lapply(graph_df[, 1:2], function(x){as.character(x)})

    # Attribute a weight to each link
    if (weight_plot == TRUE){
      graph_df$w <- igraph::E(graph)$weight
      names(graph_df) <- c("from", "to", "w")
    } else {
      names(graph_df) <- c("from", "to")
    }


    # Merge 'crds' and 'graph_df' to get the geographical coordinates
    graph_df <- merge(graph_df, crds, by.x = 'from', by.y = 'ID')
    graph_df <- merge(graph_df, crds, by.x = 'to', by.y = 'ID')

    # If 'mode = spatial', plot the nodes on a plane defined by latitude
    # and longitude axes
    # To that purpose, merge modul_df with crds.
    crds <- merge(crds, modul_df, by.x = 'ID', by.y = 'ID')

    # If links are weighted on the plot
    if (weight_plot == TRUE){
      # Give colnames to 'graph_df'
      names(graph_df) <- c("to", "from", "weight", "x", "y", "xend", "yend")

      # Transform 'w' according to 'width' option
      if (width == "w"){
        # Compute the scaled weight
        graph_df$w_sc <- sc01(graph_df$weight)
      } else if (width == "inv"){
        # Compute the inverse weight and scale it
        graph_df$w_sc <- 1/sc01(graph_df$weight)
        # Replace 'Inf' values by the largest values
        graph_df[which(graph_df$w_sc == "Inf" ),
                 'w_sc'] <- max(graph_df[which(graph_df$w_sc < Inf),
                                         "w_sc"]) + 2
        # Rescale again
        graph_df$w_sc <- sc01(graph_df$w_sc)
      } else {
        stop("You must specify a correct 'width' option.")
      }

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                                 xend = 'xend', yend = 'yend',
                                                 size = 'w_sc'),
                     color = "black")+
        scale_size_identity()+
        scale_color_manual(values = palette) +
        geom_point(data = crds, aes_string(x = 'x', y = 'y', color = 'modul'),
                   size = 6) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="Longitude", y="Latitude", color = "Module")

      # If the links are not weighted on the plot
    } else {
      # Give column names to 'graph_df'
      names(graph_df) <- c("to", "from", "x", "y", "xend", "yend")

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                                 xend = 'xend', yend = 'yend'),
                     color = "black")+
        scale_size_identity()+
        scale_fill_manual(values = palette) +
        geom_point(data = crds, aes_string(x = 'x', y = 'y', color = 'modul'),
                   size = 6) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="Longitude", y="Latitude", color = "Module")
    }


    # If the graph nodes do not have spatial coordinates
  } else if (mode == "aspatial"){

    # If 'mode = aspatial', then the nodes position is calculated with
    # Fruchterman and Reingold algorithm.
    # The links are weighted or not in the calculation,
    # according to 'weight_plot' option.

    # If links are weighted
    if (weight == TRUE){
      # Use the Fruchterman and Reingold algorithm to compute the node
      # coordinates giving a weight to the links equal to the inverse of
      # the distance they are weighted with in the graph
      crds <- igraph::layout_with_fr(graph,
                                     weights = 1/(igraph::E(graph)$weight))
    } else {
      # Use the Fruchterman and Reingold algorithm to compute
      # the node coordinates giving a weight of 1 to the links
      crds <- igraph::layout_with_fr(graph,
                                     weights = rep(1, length(igraph::E(graph))))
    }
    # Create the data.frame 'crds'
    crds <- data.frame(crds)
    names(crds) <- c("x", "y")

    # IDs of 'crds' elements are graph nodes names
    crds$ID <- igraph::V(graph)$name
    crds <- crds[, c('ID', 'x', 'y')]

    # Create a data.frame from the links of the graph
    graph_df <- data.frame(igraph::as_edgelist(graph))
    graph_df[, 1:2] <- lapply(graph_df[, 1:2], function(x){as.character(x)})

    # Add weights or not to the links
    if (weight_plot == TRUE){
      graph_df$w <- igraph::E(graph)$weight
      names(graph_df) <- c("from", "to", "w")
    } else {
      names(graph_df) <- c("from", "to")
    }

    # Merge 'graph_df' and 'crds' in order to get spatial coordinates
    # of the nodes
    graph_df <- merge(graph_df, crds, by.x = 'from', by.y = 'ID')
    graph_df <- merge(graph_df, crds, by.x = 'to', by.y = 'ID')

    # Add the modules in the crds data.frame
    crds <- merge(crds, modul_df, by.x = 'ID', by.y = 'ID')

    # If links are weighted
    if (weight_plot == TRUE){
      # Give colnames to 'graph_df'
      names(graph_df) <- c("to", "from", "weight", "x", "y", "xend", "yend")

      if (width == "w"){
        # Compute the scaled weight
        graph_df$w_sc <- sc01(graph_df$weight)
      } else if (width == "inv"){
        # Compute the inverse weight and scale it
        graph_df$w_sc <- 1/sc01(graph_df$weight)
        # Replace 'Inf' values by the largest values
        graph_df[which(graph_df$w_sc == "Inf" ),
                 'w_sc'] <- max(graph_df[which(graph_df$w_sc < Inf),
                                         "w_sc"]) + 2
        # Rescale again
        graph_df$w_sc <- sc01(graph_df$w_sc)
      } else {
        stop("You must specify a correct 'width' option.")
      }

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                                 xend = 'xend', yend = 'yend',
                                                 size = 'w_sc'),
                     color = "black")+
        scale_size_identity()+
        scale_fill_manual(values = palette) +
        geom_point(data = crds, aes_string(x = 'x', y = 'y', color = 'modul'),
                   size = 6) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="x", y="y", color = "Module")

      # If the links are not weighted
    } else {
      # Give column names to 'graph_df'
      names(graph_df) <- c("to", "from", "x", "y", "xend", "yend")

      # Create the plot
      g <- ggplot() +
        geom_segment(data = graph_df, aes_string(x = 'x', y = 'y',
                                                 xend = 'xend', yend = 'yend'),
                     color = "black")+
        scale_size_identity()+
        scale_fill_manual(values = palette) +
        geom_point(data = crds, aes_string(x = 'x', y = 'y', color = 'modul'),
                   size = 6) +
        geom_text(data = crds, aes_string(x = 'x', y = 'y', label = 'ID'),
                  size = 3, color = "black", fontface = "bold")+
        theme_bw()+
        labs(x="x", y="y", color = "Module")
    }

  } else {
    stop("You must specify a correct 'mode' option ('spatial' or 'aspatial')")
  }

  return(g)

}






