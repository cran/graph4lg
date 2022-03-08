#' Create a graph in the Graphab project
#'
#' @description The function creates a graph from a link set in a Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param linkset (optional, default=NULL) A character string indicating the
#' name of the link set used to create the graph. If \code{linkset=NULL}, every
#' link set present in the project will be used to create a graph. Link sets
#' can be created with \code{\link{graphab_link}}.
#' @param name (optional, default=NULL) A character string indicating the
#' name of the graph created. If \code{name=NULL}, a name will be created. If
#' both \code{linkset=NULL} and \code{name=NULL}, then a graph will be created
#' for every link set present in the project and a name will be created every
#' time. In the latter case, a unique name cannot be specified. Link sets
#' can be created with \code{\link{graphab_link}}.
#' @param thr (optional, default=NULL) An integer or numeric value indicating
#' the maximum distance associated with the links of the created graph. It
#' allows users to create a pruned graph based on a distance threshold. Note that
#' when the link set used has a planar topology, the graph is necessarily a
#' pruned graph (not complete) and adding this threshold parameter can remove
#' other links. When the link set has been created with cost-distances, the
#' parameter is expressed in cost-distance units whereas when the link set is
#' based upon Euclidean distances, the parameter is expressed in meters.
#' @param cost_conv FALSE (default) or TRUE. Logical indicating whether numeric
#' \code{thr} values are converted from cost-distance into Euclidean distance
#' using a log-log linear regression. See also \code{\link{convert_cd}}
#' function.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @details By default, intra-patch distances are considered for metric
#' calculation. See more information in Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_graph(proj_name = "grphb_ex",
#'               linkset = "lcp",
#'               name = "graph")
#' }

graphab_graph <- function(proj_name,         # character
                          linkset = NULL, # cost or euclid
                          name = NULL, # character
                          thr = NULL, # threshold
                          cost_conv = FALSE, # FALSE (default) or true
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
  # Check for linkset class
  if(!is.null(linkset)){
    if(!inherits(linkset, "character")){
      stop("'linkset' must be a character string")
    } else if (!(paste0(linkset, "-links.csv") %in% list.files(path = paste0(proj_path, "/", proj_name)))){
      stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
    }
  } else if (length(list.files(path = paste0(proj_path,
                                             "/", proj_name),
                               pattern = "-links.csv")) == 0){

    stop("There is not any linkset in the project you refer to.
         Please use graphab_link() before.")

  } else {

    ngraph <- length(list.files(path = paste0(proj_path, "/", proj_name),
                                pattern = "-links.csv"))
    message(paste0(ngraph, " graph(s) will be created"))
  }

  #########################################
  # Check for name
  if(!is.null(name)){
    if(!inherits(name, "character")){
      stop("'name' must be a character string")
    }
  } else if (ngraph > 1){
    stop("You cannot specify a graph name when more than one graph is created")
  }

  #########################################
  # Check for thr
  if(!is.null(thr)){
    if(!inherits(thr, c("numeric", "integer"))){
      stop("'thr' must be a numeric or an integer threshold value.")
    }
  }

  #########################################
  # Check for cost_conv
  if(!is.logical(cost_conv)){
    stop("'cost_conv' must be a logical (TRUE or FALSE).")
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

  if(!is.null(linkset)){
    cmd <- c(cmd, "--uselinkset", linkset)
  }

  cmd <- c(cmd, "--graph")

  if(!is.null(name)){
    cmd <- c(cmd, paste0("name=", name))
  }

  if(!is.null(thr)){
    if(cost_conv){
      cmd <- c(cmd, paste0("threshold={", thr, "}"))
    } else {
      cmd <- c(cmd, paste0("threshold=", thr))
    }
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
      if(file.exists(paste0(proj_path, "/", proj_name, "/", name, "-voronoi.shp"))){
        message(paste0("Graph '", name, "' has been created in the project ",
                       proj_name))
      } else {
        message("The graph creation did not succeed.")
      }
    }
  } else {
    if(file.exists(paste0(proj_path, "/", proj_name, "/", name, "-voronoi.shp"))){
      message(paste0("Graph '", name, "' has been created in the project ",
                     proj_name))
    } else {
      message("The graph creation did not succeed.")
    }
  }

}





