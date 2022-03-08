#' Computes corridors from least-cost paths already computed in
#' the Graphab project
#'
#' @description The function computes corridors around the least-cost paths
#' which have been computed in the Graphab project.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param graph A character string indicating the name of the graph with the
#' links from which the corridors are computed.
#' This graph has been created with Graphab or using \code{\link{graphab_graph}}
#' function and is associated with a link set.
#' Only the links present in the graph will be used in the computation.
#' @param maxcost An integer or numeric value indicating the maximum cost
#' distance from the least-cost paths considered for creating the corridors,
#' in cost distance units (except when \code{cost_conv = TRUE}).
#' @param format (optional, default = "raster") A character string indicating
#' whether the output is a raster file or a shapefile layer.
#' @param cost_conv FALSE (default) or TRUE. Logical indicating whether numeric
#' \code{thr} values are converted from cost-distance into Euclidean distance
#' using a log-log linear regression. See also \code{\link{convert_cd}}
#' function. Only used when \code{mode='neigh'}.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @details See more information in Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
#' Be careful, when capacity has been changed. The last changes are taken into
#' account for subsequent calculations in a project.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_corridor(proj_name = "grphb_ex",
#'                  graph = "graph",
#'                  maxcost = 1000,
#'                  format = "raster",
#'                  cost_conv = FALSE)
#' }

graphab_corridor <- function(proj_name,         # character
                             graph, # name of the graph
                             maxcost, # NULL or integer vector
                             format = "raster", # 'raster' (default) or 'vector'
                             cost_conv = FALSE, # FALSE (default) or TRUE
                             proj_path = NULL, # if NULL getwd() otherwise a character path
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

  proj_end_path <- paste0(proj_path, "/", proj_name, "/", proj_name, ".xml")


  #########################################
  # Check for graph class
  if(!inherits(graph, "character")){
    stop("'graph' must be a character string")
  } else if (!(paste0(graph, "-voronoi.shp") %in%
               list.files(path = paste0(proj_path, "/", proj_name)))){
    stop("The graph you refer to does not exist")
  } else if (length(list.files(path = paste0(proj_path,
                                             "/", proj_name),
                               pattern = "-voronoi.shp")) == 0){
    stop("There is not any graph in the project you refer to.
         Please use graphab_graph() before.")
  }

  #########################################
  # Check for maxcost
  if(!is.null(maxcost)){
    if(!inherits(maxcost, c("numeric", "integer"))){
      stop("'maxcost' must be a numeric or an integer threshold value.")
    }
  }

  #########################################
  # Check for cost_conv
  if(!is.logical(cost_conv)){
    stop("'cost_conv' must be a logical (TRUE or FALSE).")
  }

  #########################################
  # Check for format
  if(!inherits(format, "character")){
    stop("'format' must be a character string")
  } else if(!(format %in% c("raster", "vector"))){
    stop(paste0("'metric' must be equal to either 'raster' or 'vector'"))
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
  version <- "graphab-2.8.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)

  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
           "--project", proj_end_path)

  if(!is.null(graph)){
    cmd <- c(cmd, "--usegraph", graph)
  }

  cmd <- c(cmd, "--corridor")

  if(cost_conv){
    cmd <- c(cmd, paste0("maxcost={", maxcost, "}"))
  } else {
    cmd <- c(cmd, paste0("maxcost=", maxcost))
  }

  cmd <- c(cmd, paste0("format=", format))

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
      message(paste0("A ", format, " file with corridors has been created."))
    }
  } else {
    message(paste0("A ", format, " file with corridors has been created."))
  }

}
