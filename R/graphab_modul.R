#' Create modules from a graph in the Graphab project
#'
#' @description The function creates modules from a graph by maximising
#' modularity
#'
#' @inheritParams graphab_metric
#' @param graph A character string indicating the name of the graph on which
#' the modularity index is computed. This graph has been created with Graphab
#' or using \code{\link{graphab_graph}} function and is associated
#' with a link set. Only the links present in the graph and their corresponding
#' weights will be used in the computation, together with patch areas.
#' @param nb (optional, default=NULL) An integer or numeric value indicating
#' the number of modules to be created. By default, it is the number that
#' maximises the modularity index.
#' @param return Logical (default=TRUE) indicating whether results are returned
#' to user.
#' @return If \code{return=TRUE}, the function returns a message indicating
#' whether the partition has been done. New options are being developed.
#' @details This function maximises a modularity index by searching for the
#' node partition involves a large number of links within modules and a small
#' number of inter-module links. Each link is given a weight in the computation,
#' such as the weight \eqn{w_{ij}} of the link between patches i and j is:
#' \deqn{w_{ij} = (a_{i} a_{j})^\beta e^{-\alpha d_{ij}}}.
#' This function does not allow users to convert automatically Euclidean
#' distances into cost-distances.
#' See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' graphab_modul(proj_name = "graphab_example",
#'                graph = "graph_forest",
#'                dist = 1000,
#'                prob = 0.05,
#'                beta = 1)
#' }


graphab_modul <- function(proj_name, # character
                           graph, # cost or euclid
                           dist, # dist threshold
                           prob = 0.05, # dispersal probability
                           beta = 1, # area weight
                           nb = NULL, # number of components
                           return = TRUE, #
                           proj_path = NULL, # if null getwd() otherwise a character path
                           parallel.java = NULL,
                           alloc_ram = NULL){

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

  ## Create proj_end_path proj_path/proj_name/proj_name.xml
  proj_end_path <- paste0(proj_path, "/", proj_name, "/", proj_name, ".xml")

  #######################################
  # Add '' to proj_path for cases with spaces in paths
  if(all(stringr::str_sub(proj_end_path, 1, 1) != "'",
         stringr::str_sub(proj_end_path, 1, 1) != "'",
         stringr::str_detect(string = proj_end_path,
                             pattern = " "))){
    proj_end_path_cmd <- paste0("'", proj_end_path, "'")
  } else {
    proj_end_path_cmd <- proj_end_path
  }

  ## Check project version
  check_graphab_version(proj_path = proj_end_path)

  #########################################
  # Check for graph class
  if(!inherits(graph, "character")){
    stop("'graph' must be a character string")
  } else if (!check_graphab_object(proj_path = proj_end_path,
                                   object_type = "graph",
                                   name = graph)){
    stop("The graph you refer to does not exist.
           Please use graphab_graph() to create it.")
  }

  #########################################
  # Check for dist prob beta
  if(!inherits(dist, c("numeric", "integer"))){
    stop("'dist' must be a numeric or integer value")
  } else if(!inherits(prob, c("numeric", "integer"))){
    stop("'prob' must be a numeric or integer value")
  } else if(!inherits(beta, c("numeric", "integer"))){
    stop("'beta' must be a numeric or integer value")
  } else if(beta < 0 || beta > 1){
    stop("'beta' must be between 0 and 1")
  } else if(prob < 0 || prob > 1){
    stop("'prob' must be between 0 and 1")
  }

  #########################################
  # Check for nb
  if(!is.null(nb)){
    if(!inherits(nb, c("numeric", "integer"))){
      stop("'nb' must be a numeric or integer value")
    }
  }


  #########################################
  # Check for return
  if(!is.logical(return)){
    stop("'return' must be a logical (TRUE or FALSE).")
  }

  #########################################
  # Check for parallel.java
  if(!is.null(parallel.java)){
    if(!inherits(parallel.java, c("numeric", "integer"))){
      stop("'parallel.java' must be a numeric or integer value.")
    }
  }

  #########################################
  # Check for Graphab
  gr <- get_graphab(res = FALSE, return = TRUE)

  if(gr == 1){
    message("Graphab has been downloaded")
  }

  #########################################
  # Get java path
  java.path <- Sys.which("java")

  #########################################
  # Get graphab path
  version <- "graphab-3.0.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg2_jar/", version)

  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

  if(!is.null(parallel.java)){
    cmd <- c(cmd, "-proc ", as.character(parallel.java))
  }

  cmd <- c(cmd,
           "--project", proj_end_path_cmd,
           "--usegraph", graph,
           "--cluster",
           paste0("d=", dist),
           paste0("p=", prob),
           paste0("beta=", beta))

  if(!is.null(nb)){
    cmd <- c(cmd, paste0("nb=", nb))
  }

  if(!is.null(alloc_ram)){
    if(inherits(alloc_ram, c("integer", "numeric"))){
      cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
    } else {
      stop("'alloc_ram' must be a numeric or an integer")
    }
  }

  #########################################
  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)


  if(length(rs) == 1){
    if(rs == 1){
      message("An error occurred")
    } else {
      message(paste0("Clusters have been computed in the project ",
                     proj_name))
    }
  } else {
    message(paste0("Clusters have been computed in the project ",
                   proj_name))
  }

  if(return){
    res <- "Voronoi polygons corresponding to module extents have been created
            in a shapefile layer in the project directory"
    return(res)
  }

}
