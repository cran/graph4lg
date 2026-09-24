#' Merge graphs in a multiple-graph Graphab project
#'
#' @description The function creates one graph from several other graphs
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param graphs A character vector indicating the names of the graphs to be
#' merged. These graphs must already exist in the project. Graphs can be
#' created with \code{\link{graphab_graph}}.
#' @param name A character string indicating the name of the graph created
#' by merging the other graphs.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @param parallel.java An integer indicating how many computer cores are used
#' to run the .jar file. By default, \code{parallel.java = NULL}, and java sets
#' it according to local settings.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @details By default, intra-patch distances are considered for metric
#' calculation. See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_merge_graph(proj_name = "graphab_example",
#'               graphs = c("graph_forest_1", "graph_shrub_2"),
#'               name = "g_all")
#' }

graphab_merge_graph <- function(proj_name,         # character
                                graphs, # the graphs to merge
                                name = NULL, # character
                                proj_path = NULL, # if NULL getwd() otherwise a character path
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
  # Check for graphs class
  if(!inherits(graphs, "character")){
    stop("'graphs' must be a character vector.")
  } else if(!(length(graphs) >= 2)){
    stop("'graphs' must include the names of at least 2 graphs.")
  } else if(!all(unlist(lapply(graphs,
                               FUN = function(x){
                                 check_graphab_object(proj_path = proj_end_path,
                                                      object_type = "graph",
                                                      name = x)
                               })))){
    stop("'graphs' must only include graphs that exist in the project.")
  }

  #########################################
  # Check for name
  if(!inherits(name, "character")){
    stop("'name' must be a character string")
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
           "--mergegraph", paste0("name=", name),
           paste(graphs, collapse = ","))

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
    }
  }

  ## Check whether the new merged graph exists
  if(check_graphab_object(proj_path = proj_end_path,
                          object_type = "graph",
                          name = name)){
    message(paste0("Graph '", name, "' has been created in the project '",
                   proj_name, "'."))
  } else {
    message("An error occurred")
  }

}
