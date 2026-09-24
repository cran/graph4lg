#' Create a new habitat type made of meta-patches based on inter-patch
#' distances in an existing graph
#'
#' @description The function creates a new habitat type made of meta-patches
#' based on inter-patch distances in an existing pruned graph
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param graph A character string indicating the name of the graph from which
#' the meta-patches are created. This graph should have several components. Each
#' of them includes one or several habitat patches and they will each correspond
#' to a single habitat (meta-)patch in the new habitat type created by this
#' function. The graph should have been created using the
#' \code{\link{graphab_graph}} function and is associated with a link set. The
#' latter should not be complete.
#' @param mincapa An integer or numeric value indicating the minimum
#' capacity of the meta-patch to conserve in the new habitat type. This capacity
#' is computed as the sum of the capacities of all the patches included in a
#' meta-patch.
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
#' @details See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' Note that a name is automatically given to this new habitat type.
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' graphab_metapatch(proj_name = "graphab_example",
#'                  graph = "graph_forest_thr200")
#' }

graphab_metapatch <- function(proj_name,  # character
                              graph, # name of the graph
                              mincapa = NULL, # NULL or numeric value
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
  # Check for maxcost
  if(!is.null(mincapa)){
    if(!inherits(mincapa, c("numeric", "integer"))){
      stop("'mincapa' must be a numeric or an integer threshold value.")
    }
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
           "--project", proj_end_path_cmd)

  if(!is.null(graph)){
    cmd <- c(cmd, "--usegraph", graph)
  }

  cmd <- c(cmd, "--metapatch")

  if(!is.null(mincapa)){
    cmd <- c(cmd, paste0("mincapa=", mincapa))
  }

  if(!is.null(alloc_ram)){
    if(inherits(alloc_ram, c("integer", "numeric"))){
      cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
    } else {
      stop("'alloc_ram' must be a numeric or an integer")
    }
  }

  ## list habitats before metapatch
  exist_hab <- graphab_show(proj_path = proj_end_path)
  exist_hab <- exist_hab$Habitats

  #########################################
  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)

  if(length(rs) == 1){
    if(rs == 1){
      message("An error occurred")
    } else {
      new_hab <- graphab_show(proj_path = proj_end_path)
      new_hab <- new_hab$Habitats
      new_hab <- new_hab[-which(new_hab %in% exist_hab)]

      message(paste0("The new meta-patch habitat type called '", new_hab,
                     "' has been created in the project."))

    }
  } else {
    new_hab <- graphab_show(proj_path = proj_end_path)
    new_hab <- new_hab$Habitats

    new_hab <- new_hab[-which(new_hab %in% exist_hab)]

    message(paste0("The new meta-patch habitat type called '", new_hab,
                   "' has been created in the project."))

  }

}
