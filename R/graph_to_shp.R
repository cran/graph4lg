#' Export a spatial graph to shapefile layers
#'
#' @description The function enables to export a spatial graph to
#' shapefile layers.
#'
#' @param graph A graph object of class \code{igraph}
#' @param crds (if 'mode = 'spatial'') A \code{data.frame} with the spatial
#' coordinates of the graph nodes. It must have three columns:
#' \itemize{
#' \item{ID: Name of the graph nodes (will be converted into character string).
#' The names must the same as the node names of the graph object of
#' class \code{igraph} (\code{igraph::V(graph)$name})}
#' \item{x: Longitude (numeric or integer) of the graph nodes in the coordinates
#' reference system indicated with the argument crds_crs.}
#' \item{y: Latitude (numeric or integer) of the graph nodes in the coordinates
#' reference system indicated with the argument crds_crs.}
#' }
#' @param mode Indicates which shapefile layers will be created
#' \itemize{
#' \item{If 'mode = 'both'' (default), then two shapefile layers are created,
#' one for the nodes and another for the links.}
#' \item{If 'mode = 'node'', a shapefile layer is created for the nodes only.}
#' \item{If 'mode = 'link'', a shapefile layer is created for the links only.}
#' }
#' @param metrics (not considered if 'mode = 'link'') Logical. Should graph
#' node attributes integrated in the attribute table of the node shapefile
#' layer? (default: FALSE)
#' @param crds_crs An integer indicating the EPSG code of the coordinates
#' reference system to use.
#' The projection and datum are given in the PROJ.4 format.
#' @param dir_path A character string corresponding to the path to the directory
#' in which the shapefile layers will be exported. If \code{dir_path = "wd"},
#' then the layers are created in the current working directory.
#' @param layer A character string indicating the suffix of the name of
#' the layers to be created.
#' @return Create shapefile layers in the directory specified with the parameter
#' 'dir_path'.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' data(data_tuto)
#' mat_w <- data_tuto[[1]]
#' gp <- gen_graph_topo(mat_w = mat_w, topo = "gabriel")
#' crds_crs <- 2154
#' crds <- pts_pop_simul
#' layer <- "graph_dps_gab"
#' graph_to_shp(graph = gp, crds = pts_pop_simul, mode = "both",
#'              crds_crs = crds_crs,
#'              layer = "test_fonct",
#'              dir_path = tempdir(),
#'              metrics = FALSE)
#'  }



graph_to_shp <- function(graph, crds, mode = "both", crds_crs,
                         layer, dir_path,
                         metrics = FALSE){

  # Check whether 'graph' is a graph object of class igraph
  if(!inherits(graph, "igraph")){
    stop("graph must be a graph objet of class igraph")
  }

  # Check whether 'crds' is an object of class data.frame
  if(!inherits(crds, "data.frame")){
    stop("crds must be an objet of class data.frame")
  }

  # Check whether 'crds_crs' is an integer
  if(!is.null(crds_crs)){
    if(!inherits(crds_crs, c("integer", "numeric"))){
      stop("crds_crs must be an integer or a numeric value")
    }
  }

  # Check whether 'layer' is a character string
  if(!inherits(layer, "character")){
    stop("layer must be a character string")
  }

  # Check whether 'dir_path' is a character string
  if(!inherits(dir_path, "character")){
    stop("dir_path must be a character string")
  }

  # If dir_path = "wd", then get the path to the current working directory
  if(dir_path == "wd"){
    dir_path <- getwd()
  }


  # Check whether 'dir_path' is the path to an existing directory
  if(!(file.exists(dir_path))){
    stop("dir_path must be the path to an existing directory")
  }

  # Check whether 'graph' has node names
  if(is.null(igraph::V(graph)$name)){
    stop("Your graph must have node names.")
  }

  # Check whether a correct 'mode' option was specified
  if(!any(c(mode, mode, mode) == c("link", "both", "node"))){
    stop("You must specify a correct 'mode' option")
  }

  # Check whether crds and the graph object are compatible
  if( nrow(crds) != length(igraph::V(graph) ) ) {
    stop("'crds' must have the same number of rows as there are
         nodes in 'graph'")
  } else if( !all( colnames(crds) == c("ID","x","y") ) ){
    stop("Column names of crds must be 'ID', 'x' and 'y'.")
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

  # Check whether 'metrics' is a logical
  if(!is.logical(metrics)){
    stop("'metrics' must be a logical")
  } else if(metrics){
    if(mode == "link"){
      message("'metrics' argument will not be used.")
    }
  }

  # If 'mode = 'link'' or 'mode = 'both'', then export the links
  if(any(c(mode, mode) == c("link", "both"))){

    # Create a data.frame from 'graph'
    graph_df <- data.frame(igraph::as_edgelist(graph))

    # Add the link weights to the data.frame
    # (weight or 1 if the graph is not weighted)
    if(!is.null(igraph::E(graph)$weight)){
      graph_df$w <- igraph::E(graph)$weight
    } else {
      graph_df$w <- 1
    }
    names(graph_df) <- c("from", "to", "w")

    # Merge 'graph_df' and 'crds' to get the spatial coordinates of the nodes
    graph_df <- merge(graph_df, crds, by.x = 1, by.y = 1)
    graph_df <- merge(graph_df, crds, by.x = 2, by.y = 1)

    names(graph_df) <- c("to", "from", "weight", "x", "y", "xend", "yend")

    # Create two data.frames with the coordinates of the start and end nodes
    # of the links
    begin.coord <- graph_df[, c("x", "y")]
    end.coord <- graph_df[, c("xend", "yend")]

    names(begin.coord) <- names(end.coord) <- c("lon", "lat")

    # Create list of simple feature geometries (linestrings)
    l_sf <- vector("list", nrow(begin.coord))
    for (i in seq_along(l_sf)){
      l_sf[[i]] <- sf::st_linestring(as.matrix(rbind(begin.coord[i, ],
                                                     end.coord[i,])))
    }
    # Create simple feature geometry list column
    l_sfc <- sf::st_sfc(l_sf, crs = crds_crs)

    # Convert to `sp` object
    lines_sp <- suppressWarnings(methods::as(l_sfc, "Spatial"))

    # Create a data.frame with edge attributes
    edge_att <- graph_df
    row.names(edge_att) <- paste( "ID",
                                  as.character(1:nrow(edge_att)), sep ="" )

    #lines_sp@lines
    # Create a spatial lines data.frame
    link_lay <- suppressWarnings(sp::SpatialLinesDataFrame(lines_sp, edge_att))


    if(paste0("link_", layer, ".shp") %in% list.files(dir_path)){

      message(paste0("A layer named ", paste0("link_", layer, ".shp"),
                     " already exists in the directory ", dir_path,
                     ". Please remove it or use another name"))

    } else {

      # Export the shapefile layer
      sf::st_write(obj = sf::st_as_sf(link_lay),
                   dsn = dir_path,
                   layer = paste0("link_", layer),
                   driver = "ESRI Shapefile",
                   delete_layer = TRUE)

      message(paste0("Layer link_", layer, ".shp was saved in the ",
                     "following directory: ", dir_path))
    }
  }

  # If 'mode = 'node'' or 'mode = 'both'', then export the nodes
  if(any(c(mode, mode) == c("node", "both"))){

    # If 'metrics = TRUE', then add the node attributes to the shp
    if(metrics){

      attr <- as.data.frame(igraph::get.vertex.attribute(graph))

      if(ncol(attr) != 1){
        crds2 <- merge(crds, attr, by.x = "ID", by.y = "name")
      } else {
        stop("The graph does not have any attribute to include
             apart from node names, already present in crds")
      }
      # Order 'crds2' as will be reordered crds to ensure concordance
      crds2 <- crds2[order(crds2$ID),]
    }

    # Order 'crds' according to the alphabetical order of the column 'ID'
    crds <- crds[order(crds$ID),]


    # Create a data.frame with the node coordinates
    xy <- crds[,c('x','y')]

    # Create a list with the spatial points coordinates
    mxy <- as.matrix(xy)

    # Create a list of point objects
    list_pts <- list()
    for(i in 1:nrow(xy)){
      list_pts[[i]] <- sf::st_point(mxy[i, ])
    }

    # Create the point layer
    node_lay <- sf::st_sfc(list_pts, crs = crds_crs)

    # If metrics, add crds2 with the attributes

    if(metrics){
      node_lay <- sf::st_sf(node_lay,
                            crds2)
    } else {
      node_lay <- sf::st_sf(node_lay,
                            crds)
    }


    if(paste0("node_", layer, ".shp") %in% list.files(dir_path)){

      message(paste0("A layer named ", paste0("node_", layer, ".shp"),
                     " already exists in the directory ", dir_path,
                     ". Please remove it or use another name"))

    } else {

      # Export the shapefile layer
      sf::st_write(obj = node_lay,
                   dsn = dir_path,
                   layer = paste0("node_", layer),
                   driver = "ESRI Shapefile",
                   delete_layer = TRUE)

      message(paste0("Layer node_", layer, ".shp was saved in the ",
                     "following directory: ", dir_path))
    }



  }


}
