#' Plot graphs
#'
#' @description The function enables to plot graphs, whether spatial or not.
#'
#' @param graph A graph object of class \code{igraph}
#' @param crds (optional, default = NULL) If 'mode = 'spatial'', it is a
#' \code{data.frame} with the spatial coordinates of the graph nodes.
#' It must have three columns :
#' \itemize{
#' \item{ID: A character string indicating the name of the graph nodes.
#' The names must be the same as the node names of the graph of
#' class \code{igraph} (\code{igraph::V(graph)$name})}
#' \item{x: A numeric or integer indicating the longitude of the graph nodes.}
#' \item{y: A numeric or integer indicating the latitude of the graph nodes.}
#' }
#' This argument is not used when 'mode = 'aspatial'' and mandatory when 'mode =
#' 'spatial''.
#' @param mode A character string indicating whether the graph is
#' spatial ('mode = 'spatial'') or not ('mode = 'aspatial'' (default))
#' @param node_inter (optional, default = NULL) A character string indicating
#' whether the links of the graph are weighted by distances or by similarity
#' indices. It is only used when 'mode = 'aspatial'' to compute the node
#' positions with Fruchterman and Reingold algorithm. It can be equal to:
#' \itemize{
#' \item{'distance': Link weights correspond to distances. Nodes that are close
#' to each other will be close on the figure.}
#' \item{'similarity': Link weights correspond to similarity indices. Nodes that
#' are similar to each other will be close on the figure.}
#' }
#' @param link_width (optional, default = NULL) A character string indicating
#' how the width of the link is set on the figure. Their width can be:\itemize{
#' \item{inversely proportional to link weights ("inv_w", convenient with
#' distances, default)}
#' \item{proportional to link weights ("w")}
#' }
#' @param node_size (optional, default = NULL) A character string indicating
#' the graph node attribute used to set the node size on the figure. It must be
#' the name of a numeric or integer node attribute from the graph.
#' @param module (optional, default = NULL) A character string indicating
#' the graph node modules used to set the node color on the figure. It must be
#' the name of a node attribute from the graph with discrete values.
#' @param pts_col (optional, default = NULL) A character string indicating the
#' color used to plot the nodes (default: "#F2B950"). It must be a hexadecimal
#' color code or a color used by default in R. It cannot be used if 'module' is
#' specified.
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @details When the graph is not spatial ('mode = 'aspatial''),
#' the nodes coordinates are calculated with Fruchterman et Reingold algorithm.
#' The graph object \code{graph} of class \code{igraph} must have node names
#' (not necessarily in the same order as IDs in crds, given a merging is done).
#' @author P. Savary
#' @references \insertRef{fruchterman1991graph}{graph4lg}
#' @examples
#' data(pts_pop_ex)
#' data(data_ex_genind)
#' mat_w <- mat_gen_dist(data_ex_genind, dist = "DPS")
#' gp <- gen_graph_topo(mat_w = mat_w, topo = "mst")
#' g <- plot_graph_lg(graph = gp,
#'                              crds = pts_pop_ex,
#'                              mode = "spatial",
#'                              link_width = "inv_w")

plot_graph_lg <- function(graph,
                          crds = NULL,
                          mode = "aspatial",
                          node_inter = NULL,
                          link_width = NULL,
                          node_size = NULL,
                          module = NULL,
                          pts_col = NULL){


  # Check arguments
  ###########################################################################

  # Check whether graph is a graph
  if(!inherits(graph, "igraph")){
    stop("'graph' must be a graph object of class 'igraph'.")
  }

  # Check whether the graph has node names
  if(is.null(igraph::V(graph)$name)){
    stop("Your graph must have node names.")
  }

  # Check if mode is a valid argument
  if(!inherits(mode, "character")){
    stop("'mode' must be a character string.")
  } else if(!(mode %in% c("spatial", "aspatial"))){
    stop("'mode' must be either 'spatial' or 'aspatial'.")
  }


  # Check if link_width is a valid argument
  if(!is.null(link_width)){
    if(!inherits(link_width, "character")){
      stop("'link_width' must be a character string.")
    } else if(!(link_width %in% c("w", "inv_w"))){
      stop("'link_width' must be either 'w' or 'inv_w'.")
    } else {
      if(is.null(igraph::E(graph)$weight)){
        stop("Your graph must have link weights when 'link_width' is not NULL.")
      }
    }
  }


  # Check if node_size is a valid argument
  if(!is.null(node_size)){
    if(!inherits(node_size, "character")){
      stop("When specified, 'node_size' must be a character string.")
    } else if(is.null(igraph::get.vertex.attribute(graph, name = node_size))){
      stop("'node_size' doest not correspond to any graph node attribute.")
    }
  }

  # Check if module is a valid argument
  if(!is.null(module)){
    if(!inherits(module, "character")){
      stop("When specified, 'module' must be a character string.")
    } else if(is.null(igraph::get.vertex.attribute(graph, name = module))){
      stop("'module' doest not correspond to any graph node attribute.")
    }
  }

  # Check if pts_col is a valid argument
  if(!is.null(pts_col)){
    if(!inherits(pts_col, "character")){
      stop("When specified, 'pts_col' must be a character string.")
    } else if(!is.null(module)){
      stop("'pts_col' cannot be specified when 'module' is specified.")
    } else if(length(pts_col) != 1){
      stop("'pts_col' must be a character value, not a vector.")
    } else if(stringr::str_sub(pts_col, 1, 1) != "#"){
      stop("'pts_col' must be a color code.")
    }
  }

  # Check crds and compute when needed and create nodes data.frame
  ###########################################################################

  # When mode is spatial
  if (mode == "spatial"){

    # Warning if node_inter is specified
    if(!is.null(node_inter)){
      warning("When 'mode == 'spatial'', 'node_inter' argument is ignored")
    }

    # If 'crds' is given, check whether it has as many elements
    # as there are nodes
    if(is.null(crds)){
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


  } else if(mode == "aspatial"){

    # When mode is aspatial

    # Warning if crds is specified
    if(!is.null(crds)){
      warning("When 'mode == 'aspatial'', 'crds' argument is ignored")
    }

    # Check node_inter and compute crds with Fruchterman and Reingold
    if(inherits(node_inter, "character")){

      if(!(node_inter %in% c("distance", "similarity"))){
        stop("'node_inter' argument must be either 'distance', 'similarity'
           or NULL.")

      } else if (is.null(igraph::E(graph)$weight)){
        stop("Your graph must have link weights when 'node_inter' is not NULL.")
      } else if(node_inter == "distance"){

        # Vertices connected with a highly weighted link are placed closer
        # to each other. With 1/w, nodes connected with small weighted link
        # are closer.
        crds <- igraph::layout_with_fr(graph,
                                       weights = 1/(igraph::E(graph)$weight))
      } else if(node_inter == "similarity"){

        crds <- igraph::layout_with_fr(graph,
                                       weights = igraph::E(graph)$weight)

      }
    } else if(is.null(node_inter)){

      crds <- igraph::layout_with_fr(graph,
                                     weights = rep(1,
                                                   length(igraph::E(graph))))
    }

    # Create the data.frame 'crds'
    crds <- data.frame(crds)
    colnames(crds) <- c("x", "y")

    # IDs of 'crds' elements are graph nodes names
    crds$ID <- igraph::V(graph)$name
    crds <- crds[, c('ID', 'x', 'y')]

  }

  ## Node size
  if(is.null(node_size)){
    crds$n_size <- 6
  } else {
    n_size <- data.frame(ID = igraph::get.vertex.attribute(graph = graph,
                                                           name = "name"),
                         n_size = igraph::get.vertex.attribute(graph = graph,
                                                               name = node_size))
    crds <- merge(crds, n_size, by = "ID")

    if(length(unique(crds$n_size)) == 1){
      crds$n_size <- 6
    } else {
      crds$n_size <- sc01(crds$n_size)*4 + 4
    }
  }

  ## Modules
  if(is.null(module)){
    crds$module <- as.factor(1)
  } else {
    n_mod <- data.frame(ID = igraph::get.vertex.attribute(graph = graph,
                                                          name = "name"),
                        module = igraph::get.vertex.attribute(graph = graph,
                                                              name = module))
    crds <- merge(crds, n_mod, by = "ID")
    crds$module <- as.factor(crds$module)
  }

  # Create links data.frame
  ###########################################################################

  # Create a data.frame from the links of the graph
  graph_df <- data.frame(igraph::as_edgelist(graph))
  colnames(graph_df) <- c("from", "to")

  if(is.null(link_width)){
    graph_df$l_w <- 0.5
  } else if(link_width == "w"){
    if(length(unique(igraph::E(graph)$weight)) == 1){
      graph_df$l_w <- 0.5
    } else {
      graph_df$l_w <- sc01(igraph::E(graph)$weight)/2 + 0.5
    }
  } else if (link_width == "inv_w"){
    if(length(unique(igraph::E(graph)$weight)) == 1){
      graph_df$l_w <- 0.5
    } else {
      graph_df$l_w <- sc01(1/igraph::E(graph)$weight)/2 + 0.5
    }
  }

  # Merge 'graph_df' and 'crds' in order to get spatial
  # coordinates of the nodes
  graph_df <- merge(graph_df, crds[, c("ID", "x", "y")], by.x = "from", by.y = "ID")
  graph_df <- merge(graph_df, crds[, c("ID", "x", "y")], by.x = "to", by.y = "ID")


  # Give colnames to 'graph_df'
  colnames(graph_df) <- c("to", "from", "l_w", "x", "y", "xend", "yend")

  graph_df <- graph_df[, c("from", "to", "l_w", "x", "y", "xend", "yend")]



  # Create the plot
  ###########################################################################

  if(mode == "spatial"){
    xlab <- "Longitude"
    ylab <- "Latitude"
  } else {
    xlab <- "x"
    ylab <- "y"
  }

  if(is.null(pts_col)){
    pal <- mypalette
  } else {
    pal <- pts_col
  }

  g <- ggplot() +
    geom_segment(data = graph_df, aes(x = .data$x, y = .data$y,
                                      xend = .data$xend, yend = .data$yend,
                                      size = .data$l_w),
                 color = "black") +
    geom_point(data = crds, aes(x = .data$x, y = .data$y,
                                size = .data$n_size,
                                color = .data$module)) +
    geom_text(data = crds, aes(x = .data$x, y = .data$y,
                               label = .data$ID),
              size = 4, color = "black", fontface = "bold") +
    scale_size_identity() +
    scale_color_manual(values = pal) +
    theme_bw() +
    labs(x = xlab,
         y = ylab) +
    theme(legend.position = "none")


  return(g)


}


