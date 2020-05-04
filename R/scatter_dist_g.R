#' Plot scatterplots of distances to visualize the graph pruning intensity
#'
#' @description The function enables to plot scatterplots of the relationship
#' between two distances (often a genetic distance and a landscape distance
#' between populations or sample sites), while highlighting the population pairs
#' between which a link was conserved during the creation of a graph whose
#' nodes are populations (or sample sites). It thereby allows to visualize the
#' graph pruning intensity.
#'
#' @param mat_y A symmetric (complete) \code{matrix} with pairwise (genetic or
#' landscape) distances between populations or sample sites. These values
#' will be the point coordinates on the y axis. \code{mat_y} is the matrix
#' used to weight the links of the graph \code{x}, whose nodes correspond to
#' row and column names of \code{mat_y}.
#' @param mat_x A symmetric (complete) \code{matrix} with pairwise (genetic or
#' landscape) distances between populations or sample sites. These values
#' will be the point coordinates on the x axis.
#' \code{mat_x} and \code{mat_y} must have the same row and column names,
#' ordered in the same way.
#' @param graph A graph object of class \code{igraph}.
#' Its nodes must have the same names as the row and column of
#' \code{mat_y} and \code{mat_x} matrices. \code{x} must have weighted links.
#' Link weights have to be values from \code{mat_y} matrix. \code{graph} must
#' be an undirected graph.
#' @param thr_y (optional) A numeric value used to remove values from the
#' data before to plot. All values from \code{mat_y} above \code{thr_y}
#' are removed from the data.
#' @param thr_x (optional) A numeric value used to remove values from the
#' data before to plot. All values from \code{mat_x} above \code{thr_x}
#' are removed from the data.
#' @param pts_col_1 (optional) A character string indicating the color used to
#' plot the points associated to all populations or sample sites
#' pairs (default: "#999999"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @param pts_col_2 (optional) A character string indicating the color used to
#' plot the points associated to populations or sample sites pairs connected on
#' the graph (default: "black"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @details IDs in \code{mat_y} and \code{mat_x} must be the same and refer
#' to the same sampling sites or populations, and both matrices must be ordered
#' in the same way.
#' Matrices of genetic distance can be computed using
#' \code{\link{mat_gen_dist}}.
#' Matrices of landscape distance can be computed using
#' \code{\link{mat_geo_dist}} when the landscape distance needed is a
#' Euclidean geographical distance.
#' This function is based upon \code{\link{scatter_dist}} function.
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @author P. Savary
#' @examples
#' data(data_tuto)
#' mat_gen <- data_tuto[[1]]
#' mat_dist <- suppressWarnings(mat_geo_dist(data=pts_pop_simul,
#'       ID = "ID",
#'       x = "x",
#'       y = "y"))
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                      order(as.character(colnames(mat_dist)))]
#' x <- gen_graph_topo(mat_w = mat_gen, mat_topo = mat_dist, topo = "gabriel")
#' scat <- scatter_dist_g(mat_y = mat_gen, mat_x = mat_dist,
#'                        graph = x)

scatter_dist_g <- function(mat_y, mat_x,
                           graph,
                           thr_y = NULL, thr_x = NULL,
                           pts_col_1 = "#999999",
                           pts_col_2 = "black"){

  # Check whether 'mat_y' and 'mat_x' are symmetric matrices
  if(!inherits(mat_y, "matrix")){
    stop("'mat_y' must be an object of class 'matrix'.")
  } else if (!inherits(mat_x, "matrix")){
    stop("'mat_x' must be an object of class 'matrix'.")
  } else if (!Matrix::isSymmetric(mat_y)){
    stop("'mat_y' must be a symmetric pairwise matrix.")
  } else if (!Matrix::isSymmetric(mat_x)){
    stop("'mat_x' must be a symmetric pairwise matrix.")
  }

  # Check whether 'mat_y' and 'mat_x' have the same row and column names
  if(!all(row.names(mat_y) == row.names(mat_x))){
    stop("'mat_y' and 'mat_x' must have the same row names.")
  } else if(!all(colnames(mat_y) == colnames(mat_x))){
    stop("'mat_y' and 'mat_x' must have the same column names.")
  }

  # Check whether 'graph' is a graph object of class 'igraph', non-weighted,
  # non-directed, with node names and whose node names are the same as
  # row and column names of 'mat_x' and 'mat_y'
  if(!inherits(graph, "igraph")){
    stop("graph must be an object of class 'igraph'.")
  } else if (!(igraph::is.weighted(graph))){
    stop("graph must a weighted graph.")
  } else if (igraph::is.directed(graph)){
    stop("graph must a non-directed graph.")
  } else if (is.null(igraph::V(graph)$name)){
    stop("graph must have node names.")
  } else if (!all(igraph::V(graph)$name %in% row.names(mat_x))){
    stop("Node names from graph must be the same as the row and column names
        from 'mat_x'")
  }

  # Number of points
  nb_pts <- nrow(mat_y)
  # Number of nodes
  nb_nodes <- length(igraph::V(graph)$name)

  # Check whether the numbers of nodes and of points are equal
  if(nb_pts != nb_nodes){
    stop("The number of nodes in graph is not equal to the number of rows and
         columns in 'mat_y' and 'mat_x'.")
  }

  # Number of pairs (for a non-directed graph)
  nb_pairs <- nb_pts * (nb_pts - 1) / 2

  # Get the graph adjacency matrix
  graph_adj <- igraph::as_adjacency_matrix(graph, type = "both",
                                           attr = "weight",
                                           sparse = FALSE)

  # Create vectors with the link weights and the corresponding distances
  # in 'mat_x' and 'mat_y' (lower triangles of the matrices)
  if(all(row.names(graph_adj) == row.names(mat_y))){
    graph_val <- ecodist::lower(graph_adj)
  }
  y_val <- ecodist::lower(mat_y)
  x_val <- ecodist::lower(mat_x)


  # Remove values larger the thresholds if specified
  if (!is.null(thr_y)){
    if(is.numeric(thr_y)){
      if(thr_y < max(y_val)){
        y_val[which(y_val > thr_y)] <- NA
      } else {
        message("There was not any genetic distance larger than 'thr_y'.")
      }
    } else {
      stop("'thr_y' must be a numeric value.")
    }
  }
  # Remove values larger the thresholds if specified
  if (!is.null(thr_x)){
    if(is.numeric(thr_x)){
      if(thr_x < max(x_val)){
        x_val[which(x_val > thr_x)] <- NA
      } else {
        message("There was not any landscape distance larger than 'thr_x'.")
      }
    } else {
      stop("'thr_x' must be a numeric value.")
    }
  }

  # Create a data.frame to store the values
  dat <- data.frame(y_val = y_val, x_val = x_val, graph_val = graph_val)
  # Remove values when there is no link
  dat[which(dat$graph_val == 0), 'graph_val'] <- NA

  # Remove NA values from the data.frame
  if(any(is.na(dat$y_val))){
    dat <- dat[-which(is.na(dat$y_val)), ]
  }

  if(any(is.na(dat$x_val))){
    dat <- dat[-which(is.na(dat$x_val)), ]
  }

  # Create the plot
  scat <- ggplot() +
    geom_point(data = dat,
               aes_string(x = 'x_val', y = 'y_val'),
               color = pts_col_1, size = 1, shape = 16) +
    geom_smooth(data = dat,
                aes_string(x = 'x_val', y = 'y_val'),
                method = "loess", color = "black") +
    geom_point(data = dat[which(!is.na(dat$graph_val)),],
               aes_string(x = 'x_val', y = 'graph_val'),
               color = pts_col_2, size = 1.5, shape = 16) +
    labs(x = "x",
         y = "y") +
    theme_bw()

  return(scat)

}
