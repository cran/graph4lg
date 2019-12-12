#' Import landscape graphs from GRAPHAB software
#'
#' @description The function imports a landscape graph created with GRAPHAB
#' software and converts it into a graph object of class \code{igraph}.
#' The graph has weighted links and is undirected.
#' Nodes have spatial coordinates. Other nodes attributes can be included.
#' It takes shapefiles layers created with GRAPHAB as input.
#'
#' @param dir_path A character string indicating the path of the GRAPHAB project
#' directory. This directory normally contains several spatial layer files
#' in format .shp:
#' \itemize{
#' \item{the spatial layer of the habitat patches corresponding to the
#' nodes of the graph (usually named 'patches.shp').}
#' \item{(alternatively) an exported spatial layer of
#' the nodes (faster option).}
#' \item{the link spatial layer file used to import the graph.}
#' }
#' @param nodes A character string indicating the names of the node spatial
#' layer in format .shp (without extension, ex.: "nodes"
#' refers to "nodes.shp" layer).
#' This layer has been created with GRAPHAB and has therefore coordinates
#' in a projected coordinates reference system.
#' Default: nodes = "patches", referring to the spatial polygon layer of the
#' habitat patches.
#' @param links A character string indicating the name of the link spatial layer
#' in format .shp (without extension, ex.: "link_graph" refers to
#' "link_graph.shp" layer).
#' This layer has been created with GRAPHAB and has therefore coordinates
#' in a projected coordinates reference system. It includes in the attribute
#' tables between patches Euclidean as well as cost-distance. These
#' distances are used to weight the link.
#' @param weight A character string ("euc" or "cost") indicating
#' whether to weight
#' the links with Euclidean distance or cost-distance (default) values.
#' @param fig Logical (default = FALSE) indicating whether to plot a figure of
#' the resulting spatial graph. The figure is plotted using function
#' \code{\link{plot_graph_lg}}. The plotting can be long if the graph has many
#' nodes and links.
#' @param crds Logical (default = FALSE) indicating whether to create an object
#' of class \code{data.frame} with the nodes spatial coordinates. Such a
#' \code{data.frame} has 3 columns: 'ID', 'x', 'y'.
#' @details Nodes attributes can be added to the graph using the
#' function \code{add_nodes_attr}.
#' @return A graph object of class \code{igraph} (if crds = FALSE) or a
#' list of objects: a graph object of class \code{igraph} and a
#' \code{data.frame} with the nodes spatial coordinates (if crds = TRUE).
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' @examples
#' path <- system.file('extdata',package='graph4lg')
#' links <- "liens_rast2_1_11_01_19-links"
#' graph <- graphab_to_igraph(dir_path = path,
#'                            links = links,
#'                            fig = FALSE)

graphab_to_igraph <- function(dir_path,
                  nodes = "patches",
                  links,
                  weight = "cost",
                  fig = FALSE,
                  crds = FALSE){

  # Check whether the input are character strings
  if(!all(c(inherits(dir_path, "character"),
            inherits(nodes, "character"),
            inherits(links, "character"),
            inherits(weight, "character")))){
    stop("Inputs 'dir_path', 'nodes', 'links' and 'weight' must
         be character strings.")
  }

  #sink("aux")
  # Open the patches layer and get the attribute table as a data.frame
  patches <- rgdal::readOGR(dsn = dir_path, layer = nodes)
  nds_df <- data.frame(patches@data)

  # Open the links layer and get the attribute table as a data.frame
  links1 <- rgdal::readOGR(dsn = dir_path, layer = links)
  links1_df <- data.frame(links1@data)
  #sink(NULL)


  # If as many node ID in nds_df and in links1, then there is not any
  # isolated node and the graph can be created directly from the edge list
  # derived from links1
  if(length(unique(nds_df$Id)) == length(unique(c(links1_df$ID1,
                                                  links1_df$ID2)))){

    edge_list <- as.matrix(links1_df[, c('ID1','ID2')])
    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    # We add a weight to the links of the complete graph
    # which have many links with null weights
    if (weight == "cost"){
      igraph::E(graph)$weight <- links1_df[, "Dist"]
    } else if (weight == "euc"){
      igraph::E(graph)$weight <- links1_df[, "DistM"]
    } else {
      stop("You must specify a correct 'weight' option ('cost' or 'euc').")
    }

    # We name the nodes of the graph.
    # In V(graph), the nodes are in the increasing order of their number ID
    igraph::V(graph)$name <- as.character(1:length(igraph::V(graph)))

  } else {
    # There are isolated nodes to include in the graph even if they are not
    # in the edge list derived from links1.
    # The function is then slower.

    # We create a vector with the number ID of the patches.
    veca <- as.character(1:nrow(nds_df))
    # We create a data.frame with all the unique possible combinations
    # of patches linked by a potential link
    df <- data.frame(expand.grid(veca, veca))
    df[,1:2] <- lapply(df[, 1:2], function(x){as.numeric(as.character(x))})
    # We delete lines if ID1 <= ID2 to retain only unique combinations
    df <- df[- which(df$Var1 <= df$Var2),]
    df[,1:2] <- lapply(df[, 1:2], function(x){as.character(x)})
    # The unique Id is given by "ID1-ID2" with ID1 > ID2 as in GRAPHAB
    df$Id <- paste(df$Var1, "-", df$Var2, sep = "")
    df$DistM <- df$Dist <- rep(0, nrow(df))
    names(df)[1:2] <- c("ID1", "ID2")
    df <- df[, c(3,1:2,4:5)]

    # df should have the same column names as links1
    if(all(names(links1@data) == names(df))){
      df <- df[-which(df$Id %in% links1_df$Id), ]
      df <- rbind(df, links1_df)
    } else {
      stop("Error probably due to unusual structure
           of the links spatial layer.")
    }

    # We extract the edgelist and create a complete unweighted graph
    edge_list <- as.matrix(df[, c('ID1','ID2')])
    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    # We add a weight to the links of the complete graph
    # which have many links with null weights
    if (weight == "cost"){
      igraph::E(graph)$weight <- df[, "Dist"]
    } else if (weight == "euc"){
      igraph::E(graph)$weight <- df[, "DistM"]
    } else {
      stop("You must specify a correct 'weight' option ('cost' or 'euc').")
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


  coords <- data.frame(cbind(patches$Id, sp::coordinates(patches)))
  names(coords) <- c("ID", "x", "y")

  if(fig){
    plot_spg <- plot_graph_lg(graph, crds = coords)
    print(plot_spg)
  }


  if(crds){
    res <- list(graph, coords)
  } else {
    res <- graph
  }

  return(res)

}

