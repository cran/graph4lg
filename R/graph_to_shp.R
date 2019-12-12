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
#' @param metrics (not possible if 'mode = 'link'') Logical. Should metrics
#' be calculated and integrated in the attribute table of the node shapefile
#' layer? (default: FALSE)
#' Metrics calculated are degrees, betweenness centrality and sum of
#' inverse weight (if links are weighted)
#' @param crds_crs A character string indicating the Coordinates
#' Reference System of the spatial coordinates of the nodes and of the
#' shapefile layers created.
#' The projection and datum are given in the PROJ.4 format.
#' @param dir_path A character string corresponding to the path to the directory
#' in which the shapefile layers will be exported. If \code{dir_path = "wd"},
#' then the layers are created in the current working directory.
#' @param layer_name A character string indicating the suffix of the name of
#' the layers to be created.
#' @return Create shapefile layers in the directory specified with the parameter
#' 'dir_path'.
#' @export
#' @author P. Savary
#' @examples
#' data(data_tuto)
#' mat_w <- data_tuto[[1]]
#' gp <- gen_graph_topo(mat_w = mat_w, topo = "gabriel")
#' crds_crs1 <- "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 "
#' crds_crs2 <- "+x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"
#' crds_crs <- paste(crds_crs1, crds_crs2, sep = "")
#' crds <- pts_pop_simul
#' layer_name <- "graph_dps_gab"
#' graph_to_shp(graph = gp, crds = pts_pop_simul, mode = "both",
#'              crds_crs = crds_crs,
#'              layer_name = "test_fonct",
#'              dir_path = tempdir(),
#'              metrics = TRUE)


graph_to_shp <- function(graph, crds, mode = "both", crds_crs,
                         layer_name, dir_path,
                         metrics = FALSE){

  # Check whether 'graph' is a graph object of class igraph
  if(!inherits(graph, "igraph")){
    stop("graph must be a graph objet of class igraph")
  }

  # Check whether 'crds' is an object of class data.frame
  if(!inherits(crds, "data.frame")){
    stop("crds must be an objet of class data.frame")
  }

  # Check whether 'crds_crs' is a character string
  if(!inherits(crds_crs, "character")){
    stop("crds_crs must be a character string")
  }

  # Check whether 'layer_name' is a character string
  if(!inherits(layer_name, "character")){
    stop("layer_name must be a character string")
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
    lines_sp <- methods::as(l_sfc, "Spatial")

    # Create a data.frame with edge attributes
    edge_att <- graph_df
    row.names(edge_att) <- paste( "ID",
                                  as.character(1:nrow(edge_att)), sep ="" )

    #lines_sp@lines
    # Create a spatial lines data.frame
    link_lay <- sp::SpatialLinesDataFrame(lines_sp, edge_att)

    # Add the CRS to 'link_lay'
    raster::crs(link_lay) <- crds_crs

    # Export the shapefile layer
    rgdal::writeOGR(link_lay,
                    dsn = dir_path,
                    layer = paste("link", layer_name, sep = "_"),
                    driver = "ESRI Shapefile")

  }

  # If 'mode = 'node'' or 'mode = 'both'', then export the nodes
  if(any(c(mode, mode) == c("node", "both"))){

    # Order 'crds' according to the alphabetical order of the column 'ID'
    crds <- crds[order(crds$ID),]

    # If 'metrics = TRUE', then compute three node-level metrics
    if (metrics){
      crds$btw <- igraph::betweenness(graph)
      crds$degree <- igraph::degree(graph)

      if(!is.null(igraph::E(graph)$weight)){
        crds$siw <- igraph::strength(graph,
                                     weights = 1/igraph::E(graph)$weight)
      }
    }
    # Create a data.frame with the node coordinates
    xy <- crds[,c('x','y')]
    # Create a spatial points data.frame
    node_lay <- sp::SpatialPointsDataFrame(coords = xy, data = crds,
                                         proj4string = sp::CRS(crds_crs))


    # Export the shapefile layer
    rgdal::writeOGR(node_lay,
                    dsn = dir_path,
                    layer = paste("node", layer_name, sep = "_"),
                    driver = "ESRI Shapefile")
  }

  message(paste("Created layer(s) were saved in the following directory: ",
              dir_path, sep = ""))
}
