#' Create landscape graphs from Graphab link set
#'
#' @description The function creates a landscape graph from a link set created
#' with Graphab software or different functions of this package and converts
#' it into a graph object of class \code{igraph}.
#' The graph has weighted links and is undirected.
#' Nodes attributes present in the Graphab project are included, including
#' connectivity metrics when computed
#'
#' @param proj_name A character string indicating the project name. It is also
#' the name of the directory in which proj_name.xml file is found. By default,
#' 'proj_name' is searched into the current working directory
#' @param linkset A character string indicating the name of the linkset used to
#' create the graph links. The linkset must have been created previously (see
#' the function \code{\link{graphab_link}}). It can be complete or planar. The
#' graph is given the topology of the selected link set.
#' @param nodes A character string indicating whether the nodes of the created
#' graph are given all the attributes or metrics computed in Graphab or only
#' those specific to a given graph previously created with
#' \code{\link{graphab_graph}}
#' It can be:\itemize{
#' \item{\code{nodes = "patches"}(default): all the attributes and metrics of
#' the habitat patches are included as node attributes in \code{igraph} object.}
#' \item{\code{nodes = "graph_name"}(default): only the metrics of
#' the habitat patches computed from the graph 'graph_name' created with
#' \code{\link{graphab_graph}} are included as node attributes in
#' \code{igraph} object, along with some basic patch attributes.}
#' }
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
#' proj_path <- system.file('extdata',package='graph4lg')
#' proj_name <- "grphb_ex"
#' linkset <- "lkst1"
#' nodes <- "graph"
#' graph <- graphab_to_igraph(proj_name = proj_name,
#'                            linkset = "lkst1",
#'                            nodes = "graph",
#'                            links = links,
#'                            weights = "cost",
#'                            proj_path = proj_path,
#'                            crds = FALSE,
#'                            fig = FALSE)
#'                            }


graphab_to_igraph <- function(proj_name,
                              linkset,
                              nodes = "patches",
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
    chg <- 1
    wd1 <- getwd()
    setwd(dir = proj_path)
  } else {
    chg <- 0
    proj_path <- getwd()
  }

  #########################################
  # Check for proj_name class
  if(!inherits(proj_name, "character")){
    stop("'proj_name' must be a character string")
  } else if (!(paste0(proj_name, ".xml") %in% list.files(path = paste0("./", proj_name)))){
    stop("The project you refer to does not exist.
         Please use graphab_project() before.")
  }


  #########################################
  # Check for linkset class
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (length(list.files(path = paste0("./", proj_name),
                               pattern = "-links.csv")) == 0){
    stop("There is not any linkset in the project you refer to.
         Please use graphab_link() before.")
  } else if (!(paste0(linkset, "-links.csv") %in% list.files(path = paste0("./",
                                                                           proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  #########################################
  # Check for nodes
  if(!inherits(nodes, "character")){
    stop("'nodes' must be a character string.")
  } else if(!file.exists(paste0("./", proj_name, "/patches.csv"))){
    df_nodes <- foreign::read.dbf(file = paste0("./", proj_name, "/patches.dbf"))
  } else {
    df_nodes <- utils::read.csv(file = paste0("./", proj_name, "/patches.csv"))
  }

  #########################################
  # Select nodes columns if nodes == graph_name

  if(nodes != "patches"){

    if (length(list.files(path = paste0("./", proj_name), pattern = "-voronoi.shp")) == 0){
      stop("There is not any graph in the project you refer to.
         Please use graphab_graph() before.")
    } else if(!(paste0(nodes,
                       "-voronoi.shp") %in% list.files(path = paste0("./",
                                                                     proj_name)))){
      stop("The graph you refer to does not exist")

    } else {

      char_graph <- nchar(nodes)
      col_graph <- which(stringr::str_sub(colnames(df_nodes),
                                          -(char_graph + 1),
                                          -1) == paste0("_", nodes))
      df_nodes <- df_nodes[, c(1:4, col_graph)]

    }

  }
  #######################################
  # Get links

  df_links <- get_graphab_linkset(proj_name = proj_name,
                                  linkset = linkset)

  # If as many node ID in df_nodes and in df_links, then there is not any
  # isolated node and the graph can be created directly from the edge list
  # derived from df_links
  if(length(unique(df_nodes$Id)) == length(unique(c(df_links$ID1,
                                                    df_links$ID2)))){
    edge_list <- as.matrix(df_links[, c('ID1','ID2')])
    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    # We add a weight to the links of the complete graph
    # which have many links with null weights
    if (weight == "cost"){
      igraph::E(graph)$weight <- df_links[, "Dist"]
    } else if (weight == "euclid"){
      igraph::E(graph)$weight <- df_links[, "DistM"]
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
    veca <- as.character(1:nrow(df_nodes))
    # We create a data.frame with all the unique possible combinations
    # of patches linked by a potential link
    df <- data.frame(expand.grid(veca, veca))
    df[, 1:2] <- lapply(df[, 1:2], function(x){as.numeric(as.character(x))})
    # We delete lines if ID1 <= ID2 to retain only unique combinations
    df <- df[- which(df$Var1 <= df$Var2),]
    df[, 1:2] <- lapply(df[, 1:2], function(x){as.character(x)})
    # The unique Id is given by "ID1-ID2" with ID1 < ID2 as in Graphab
    df$Id <- paste0(df$Var1, "-", df$Var2)
    df$DistM <- df$Dist <- rep(0, nrow(df))
    colnames(df)[1:2] <- c("ID1", "ID2")
    df <- df[, c("Id", "ID1", "ID2", "Dist", "DistM")]

    # df should have the same column names as df_links
    if(all(colnames(df_links) == colnames(df))){
      df <- df[-which(df$Id %in% df_links$Id), ]
      df <- rbind(df, df_links)
    } else {
      stop("Error probably due to unusual structure
           of the links spatial layer.")
    }

    # We extract the edgelist and create a complete unweighted graph
    edge_list <- as.matrix(df[, c('ID1', 'ID2')])
    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)

    # We add a weight to the links of the complete graph
    # which have many links with null weights
    if (weight == "cost"){
      igraph::E(graph)$weight <- df[, "Dist"]
    } else if (weight == "euclid"){
      igraph::E(graph)$weight <- df[, "DistM"]
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
                          data = df_nodes, index = "Id")


  if(crds){

    #sp_patches <- suppressWarnings(rgdal::readOGR(dsn = paste0(getwd(), "/", proj_name),
    #                                             layer = "patches"))

    sp_patches <- suppressWarnings(sf::as_Spatial(sf::st_read(dsn = paste0(getwd(),
                                                                           "/", proj_name),
                                                              layer = "patches")))

    coords <- data.frame(cbind(sp_patches$Id,
                               sp::coordinates(sp_patches)))
    names(coords) <- c("ID", "x", "y")

  }

  if(fig){
    if(crds){
      plot_spg <- plot_graph_lg(graph, mode = "spatial",
                                crds = coords,
                                node_size = "Area",
                                link_width = "inv_w")
    } else {
      plot_spg <- plot_graph_lg(graph, mode = "aspatial",
                                link_width = "inv_w",
                                node_inter = "distance",
                                node_size = "Area")
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


