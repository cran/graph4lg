#' Add attributes to the nodes of a graph
#'
#' @description The function adds attributes to the nodes of a graph from
#' either an object of class \code{data.frame} or from a shapefile layer.
#' The node IDs in the input objects must be the same as in the graph object.
#'
#' @param graph A graph object of class \code{igraph}.
#' @param input A character string indicating the nature of the
#' input data from which come the attributes to add to the nodes.
#' \itemize{
#' \item{If 'input = "shp"', then attributes come from the attribute table of
#' a shapefile layer of type point.}
#' \item{If 'input = "df"', then attributes come from an object of class
#' \code{data.frame}}
#' }
#' In both cases, input attribute table or dataframe must have a column with
#' the exact same values as the node IDs.
#' @param data (only if 'input = "df"') The name of the object of
#' class \code{data.frame} with the attributes to add to the nodes.
#' @param dir_path (only if 'input = "shp"') The path (character string) to the
#' directory containing the shapefile layer of type point whose attribute
#' table contains the attributes to add to the nodes.
#' @param layer (only if 'input = "shp"') The name (character string) of the
#' shapefile layer of type point (without extension, ex.: "nodes" refers
#' to "nodes.shp" layer) whose attribute table contains the attributes
#' to add to the nodes.
#' @param index The name (character string) of the column with the nodes names
#' in the input data (column of the attribute table or of the dataframe).
#' @param include A character string (vector) indicating which columns of the
#' input data will be added as nodes' attributes.
#' By default, 'include = "all"', i.e. every column of the input data is added.
#' Alternatively, 'include' can be a vector with the names of the columns to add
#' (ex.: "c('x', 'y', 'pop_name')").
#' @details The graph can be created with the function
#' \code{\link{graphab_to_igraph}} from shapefile layers created with GRAPHAB.
#' Values of the metrics computed at the node level with GRAPHAB can then be
#' added to such a graph with this function.
#' @return A graph object of class \code{igraph}
#' @export
#' @author P. Savary
#' @examples
#' data("data_tuto")
#' graph <- data_tuto[[3]]
#' df_nodes <- data.frame(Id = igraph::V(graph)$name,
#'                        Area = runif(50, min = 10, max = 60))
#' graph <- add_nodes_attr(graph,
#'                         data = df_nodes,
#'                         input = "df",
#'                         index = "Id",
#'                         include = "Area")

add_nodes_attr <- function(graph,
                           input = "df",
                           data,
                           dir_path = NULL,
                           layer = NULL,
                           index = "Id",
                           include = "all"){

  # Check whether graph is a graph of class igraph
  if(!inherits(graph, "igraph")){
    stop("'graph' must be an object of class 'igraph'.")
  # and check if nodes have names
  } else if (is.null(igraph::V(graph)$name)){
    stop("'graph' must have nodes' names.")
  }

  # Create a vector with node names
  nds.names <- as.character(igraph::V(graph)$name)

  # Check whether input, data and index are compatible
  if(input == "df"){
    if(!inherits(data, "data.frame")){
      stop("'data' must be a data.frame when 'input = 'df''.")
    } else if (!(index %in% names(data))){
      stop("'index' must be the name of a column of 'data'.")
    }
  # Check whether input, dir_path and layer are compatible
  } else if (input == "shp"){

    if(any(c(is.null(dir_path),
             is.null(layer)))){
      stop("'dir_path' and 'layer' must be character strings when
           'input = 'shp''.")
    } else if(!all(c(inherits(dir_path, "character"),
                     inherits(layer, "character")))){
      stop("'dir_path' and 'layer' must be character strings when
           'input = 'shp''.")
    } else {
      # If 'dir_path' and 'layer' are well defined, open the GIS layer
      sink("aux")
      data <- rgdal::readOGR(dsn = dir_path, layer = layer)
      sink(NULL)
      # Get the attribute table of the layer as a data.frame
      data <- data.frame(data@data)
      # Check whether the attribute table contains a column named 'index'
      if (!(index %in% names(data))){
        stop("'index' must be the name of a column of the attribute table
             of 'layer'.")
      }
    }

  } else {
    stop("You must specify a correct 'input' option ('df' or 'shp').")
  }

  data.names <- as.character(data[, index])

  # Check whether node names are in data.names
  if(!all(nds.names %in% data.names)){
    stop("Column 'index' from input data must contain the nodes names
         of 'data'.")
  }

  # Get only the data relative to the nodes
  data <- data[which(data.names %in% nds.names), ]
  # Reorder data in the same order as the graph nodes
  data <- data[match(nds.names, data[, index]), ]
  # attrib are all the columns of data different from index
  attrib <- setdiff( names(data), index )

  if(inherits(include, "character")){
    # If include indicates more than one column
    # get the set of corresponding variables
    if(length(include) > 1){
      attrib <- attrib[which(attrib %in% include)]
    # If include = "all", attrib is not modified
    } else if (include == "all"){
      NULL
    # If include has one element, different from "all", attrib is reduced
    # to this element.
    } else {
      attrib <- attrib[which(attrib %in% include)]
    }

  } else {
    stop("'include' must be a character string or a vector
         of character strings.")
  }

  # If attrib does not contain anything, it stops.
  if(length(attrib) == 0){
    stop("Elements of 'include' must be attributes names from input data.")
  }

  # Add attributes to the graph's nodes from data and according to attrib
  for (i in 1:length(attrib)){
    graph <- igraph::set_vertex_attr(graph, attrib[i],
                                     value = data[, attrib[i]])
  }

  return(graph)

}
