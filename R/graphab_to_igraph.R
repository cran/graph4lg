#' Create landscape graphs from Graphab link set
#'
#' @description The function creates a landscape graph from a link set and a
#' set of habitat patches created with Graphab and converts it into a graph
#' object of class \code{igraph}. The graph has weighted links and
#' is undirected. Nodes attributes present in the habitat layer created with
#' Graphab project are included, including connectivity metrics when computed.
#'
#' @param proj_name A character string indicating the project name. It is also
#' the name of the directory in which proj_name.xml file is found. By default,
#' 'proj_name' is searched into the current working directory
#' @param linkset A character string indicating the name of the linkset used to
#' create the graph links. The linkset must have been created previously (see
#' the function \code{\link{graphab_link}}). It can be complete or planar. The
#' graph is given the topology of the selected link set.
#' @param habitat A character string indicating the name of a habitat type
#' created in the Graphab project.
#' @param nodes Deprecated parameter from graph4lg < 2.0.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory ('proj_name'). By default,
#' 'proj_name' is searched into the current working directory
#' @param weight A character string ("euclid" or "cost") indicating
#' whether to weight the links with Euclidean distance or
#' cost-distance (default) values.
#' @param fig Logical (default = FALSE) indicating whether to plot a figure of
#' the resulting spatial graph. The figure is plotted using function
#' \code{\link{plot_graph_lg}}. The plotting can be long if the graph has many
#' nodes and links.
#' @param crds Logical (default = FALSE) indicating whether to create an object
#' of class \code{data.frame} with the node centroid spatial coordinates. Such a
#' \code{data.frame} has 3 columns: 'ID', 'x', 'y'.
#' @return A graph object of class \code{igraph} (if crds = FALSE) or a
#' list of objects: a graph object of class \code{igraph} and a
#' \code{data.frame} with the nodes spatial coordinates (if crds = TRUE).
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' @examples
#' \dontrun{
#' proj_path <- system.file('extdata', package = 'graph4lg')
#' proj_name <- "graphab_example"
#' linkset <- "forest_link_planar"
#' habitat <- "forest"
#' graph <- graphab_to_igraph(proj_name = proj_name,
#'                            linkset = linkset,
#'                            habitat = habitat,
#'                            weights = "cost",
#'                            proj_path = proj_path,
#'                            crds = FALSE,
#'                            fig = FALSE)
#'                            }


graphab_to_igraph <- function(proj_name,
                              linkset,
                              habitat,
                              nodes = NULL,
                              weight = "cost",
                              proj_path = NULL,
                              fig = FALSE,
                              crds = FALSE){

  # Check whether the input fig and crds are logical
  if(!all(c(is.logical(fig),
            is.logical(crds)))){
    stop("Inputs 'fig' and 'crds' must be TRUE or FALSE")
  }


  #########################################
  # Check for project directory path
  if(!is.null(proj_path)){
    if(!dir.exists(proj_path)){
      stop(paste0(proj_path, " is not an existing directory or the path is ",
                  "incorrectly specified."))
    } else {
      proj_path <- normalizePath(proj_path)
    }
  } else {
    proj_path <- normalizePath(getwd())
  }

  #########################################
  # Check for proj_name class
  if(!inherits(proj_name, "character")){
    stop("'proj_name' must be a character string")
  } else if (!(paste0(proj_name, ".xml") %in%
               list.files(path = paste0(proj_path, "/", proj_name)))){
    stop("The project you refer to does not exist.
         Please use graphab_project() before.")
  }

  proj_end_path <- paste0(proj_path, "/", proj_name, "/", proj_name, ".xml")

  ## Check project version
  check_graphab_version(proj_path = proj_end_path)

  #########################################
  # Check for linkset class
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string specifying the name of the
         first link set involved in the comparison.")
  } else if (!(check_graphab_object(proj_path = proj_end_path,
                                    object_type = "linkset",
                                    name = linkset))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  #########################################
  # Check for habitat class and mode compatibility
  if(!inherits(habitat, "character")){
    stop("'habitat' must be a character string")
  } else if(!(check_graphab_object(proj_path = proj_end_path,
                                   object_type = "habitat",
                                   name = habitat))){
    stop("The habitat type you refer to does not exist.
         Please use graphab_habitat() before.")
  }

  #########################################
  # Check for nodes
  if(!is.null(nodes)){
    stop(paste0("Argument 'nodes' is deprecated and not used anymore ",
                "in graph4lg >= 2.0. You now need to provide a habitat type."))
  }

  #########################################
  # Load the habitat table with or without coordinates

  # Load the habitat patch layer
  patches <- suppressWarnings(
    sf::read_sf(paste0(proj_path, "/", proj_name, "/",
                       habitat, "/patches.gpkg"),
                as_tibble = FALSE))
  if(crds){
    # Get its centroid coordinates
    coords <- suppressWarnings(
      data.frame(
      sf::st_coordinates(
        sf::st_centroid(patches))))
    # Rename and re-arrange
    coords$ID <- patches$Id
    coords <- coords[, c("ID", "X", "Y")]
    colnames(coords) <- c("ID", "x", "y")
  }
  # Make patches aspatial from now on
  patches <- data.frame(
    sf::st_drop_geometry(patches[, -which(colnames(patches) == "the_geom")]))

  #########################################
  # Load the linkset
  links <- get_graphab_linkset(proj_name = proj_name,
                                  linkset = linkset,
                                  proj_path = proj_path)

  #########################################
  # Check whether the linkset includes only links of this habitat type
  if(!all(links$id1 %in% patches$Id)){
    stop(paste0("You must include a linkset defined only for this habitat. ",
                "You are providing a linkset connecting multiple habitats."))
  } else if(!all(links$id2 %in% patches$Id)){
    stop(paste0("You must include a linkset defined only for this habitat. ",
                "You are providing a linkset connecting multiple habitats."))
  }

  # If as many node ID in df_nodes and in df_links, then there is not any
  # isolated node and the graph can be created directly from the edge list
  # derived from df_links
  if(length(unique(patches$Id)) == length(unique(c(links$id1,
                                                    links$id2)))){
    edge_list <- as.matrix(links[, c('id1','id2')])
    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    # We add a weight to the links of the complete graph
    # which have many links with null weights
    if (weight == "cost"){
      igraph::E(graph)$weight <- links[, "dist"]
    } else if (weight == "euclid"){
      igraph::E(graph)$weight <- links[, "distm"]
    } else {
      stop("You must specify a correct 'weight' option ('cost' or 'euclid').")
    }

    # We name the nodes of the graph.
    # In V(graph), the nodes are in the increasing order of their number ID
    igraph::V(graph)$name <- as.character(1:length(igraph::V(graph)))

  } else {
    # There are isolated nodes to include in the graph even if they are not
    # in the edge list derived from df_links.
    # The function is then slower.

    # We create a vector with the number ID of the patches.
    veca <- as.character(1:nrow(patches))
    # We create a data.frame with all the unique possible combinations
    # of patches linked by a potential link
    df <- data.frame(expand.grid(veca, veca))
    df[, 1:2] <- lapply(df[, 1:2], function(x){as.numeric(as.character(x))})
    # We delete lines if ID1 <= ID2 to retain only unique combinations
    df <- df[- which(df$Var1 <= df$Var2),]
    df[, 1:2] <- lapply(df[, 1:2], function(x){as.character(x)})
    # The unique Id is given by "ID1-ID2" with ID1 < ID2 as in Graphab
    df$Id <- paste0(df$Var1, "-", df$Var2)
    df$distm <- df$dist <- rep(0, nrow(df))
    colnames(df)[1:2] <- c("id1", "id2")
    df <- df[, c("Id", "id1", "id2", "dist", "distm")]

    # df should have the same column names as df_links
    if(all(colnames(links) == colnames(df))){
      df <- df[-which(df$Id %in% links$Id), ]
      df <- rbind(df, links)
    } else {
      stop("Error probably due to unusual structure
           of the links spatial layer.")
    }

    # We extract the edgelist and create a complete unweighted graph
    edge_list <- as.matrix(df[, c('id1', 'id2')])
    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)

    # We add a weight to the links of the complete graph
    # which have many links with null weights
    if (weight == "cost"){
      igraph::E(graph)$weight <- df[, "dist"]
    } else if (weight == "euclid"){
      igraph::E(graph)$weight <- df[, "distm"]
    } else {
      stop("You must specify a correct 'weight' option ('cost' or 'euclid').")
    }

    # We extract the adjacency matrix of the first graph
    graph_mat <- igraph::as_adjacency_matrix(graph,
                                             type = "both",
                                             attr = "weight",
                                             sparse = FALSE)
    # We create another graph with as many nodes as the first one (number
    # of patches) but we only retained links with weights above 0.
    graph <- igraph::graph_from_adjacency_matrix(graph_mat,
                                                 weighted = TRUE,
                                                 mode = "undirected",
                                                 diag = FALSE)

    # In that case, the nodes of the graph are already named
  }

  igraph::V(graph)$name <- 1:length(igraph::V(graph))

  graph <- add_nodes_attr(graph = graph, input = "df",
                          data = patches, index = "Id")

  if(fig){
    if(crds){
      plot_spg <- plot_graph_lg(graph, mode = "spatial",
                                crds = coords,
                                node_size = "area",
                                link_width = "inv_w")
    } else {
      plot_spg <- plot_graph_lg(graph, mode = "aspatial",
                                link_width = "inv_w",
                                node_inter = "distance",
                                node_size = "area")
    }
    print(plot_spg)
  }

  if(crds){
    res <- list(graph, coords)
  } else {
    res <- graph
  }

  return(res)

}


