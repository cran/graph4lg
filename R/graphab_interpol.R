#' Creates a raster with interpolated connectivity metric values from metrics
#' already computed in the Graphab project
#'
#' @description The function creates a raster with interpolated connectivity
#' metric values from a metric already computed in the Graphab project.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param name A character string indicating the name of the raster to be
#' created after the interpolation.
#' @param reso An integer indicating the spatial resolution in meters of the
#' raster resulting from the metric interpolation.
#' @param linkset A character string indicating the name of the link set used
#' for the interpolation. It should be the one used to create the used graph
#' and the metric.
#' @param graph A character string indicating the name of the graph from which
#' the metric was computed and whose links are considered for a potential
#' multi-linkage with patches.
#' This graph has been created with Graphab or using \code{\link{graphab_graph}}
#' function and is associated with a link set.
#' @param var A character string indicating the name of the already computed
#' metric to be interpolated.
#' @param dist A numeric or integer value specifying the distance at which we
#' assume a probability equal to \code{prob} during the interpolation.
#' It is used to set \eqn{\alpha} for computing probabilities associated
#' with distances between each pixel and the neighboring patch(es) such that
#' probability between patch i and pixel j is \eqn{p_{ij}= e^{-\alpha d_{ij}}}.
#' @param prob A numeric or integer value specifying the probability
#' at distance \code{dist}. By default, \code{code=0.05}. It is used to set
#' \eqn{\alpha} (see param \code{dist} above).
#' @param thr (default NULL) If NULL, the value of each pixel is computed from
#' the value of the metric at the nearest habitat patch, weighted by a
#' probability depending on distance. If an integer, the value of each pixel
#' depends on the values of the metric taken at several of the nearest habitat
#' patches, up to a distance (cost or Euclidean distance, depending on the type
#' of linkset) equal to \code{thr}.
#' @param summed Logical (default = FALSE) only used if \code{thr} is not NULL,
#' and specifying whether multiple values are summed up (TRUE) or averaged
#' after being weighted by probabilities.
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
#' Be careful, when capacity has been changed. The last changes are taken into
#' account for subsequent calculations in a project.
#' This function does not work for multiple habitat graphs.
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' graphab_interpol(proj_name = "graphab_example",
#'                  name = "F_interp",
#'                  reso = 20,
#'                  linkset = "lcp",
#'                  graph = "graph_forest",
#'                  var = "F_d600_p0.5_beta1_graph",
#'                  dist = 600,
#'                  prob = 0.5)
#' }

graphab_interpol <- function(proj_name,   # character
                             name, # name of the created raster
                             reso, # integer
                             linkset, # name of the linkset
                             graph, # name of the graph
                             var, # name of the variable
                             dist, # dist threshold
                             prob = 0.05, # dispersal probability
                             thr = NULL, # NULL or integer depending on multi
                             summed = FALSE, # logical, only used if thr is integer
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
  # Check for name class
  if(!inherits(name, "character")){
    stop("'name' must be a character string")
  }

  #########################################
  # Check for linkset class
  if(!is.null(linkset)){
    if(!inherits(linkset, "character")){
      stop("'linkset' must be a character string")
    } else if (!check_graphab_object(proj_path = proj_end_path,
                                     object_type = "linkset",
                                     name = linkset)){
      stop("The linkset you refer to does not exist.
           Please use graphab_link() to create it.")
    }
  }

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
  # Check for var class
  if(!inherits(var, "character")){
    stop("'var' must be a character string")
  } else {

    ## Project content
    proj_content <- graphab_show(proj_path = proj_end_path)
    if(!(var %in% proj_content$Metrics)){
      stop("'var' must be the name of an already computed metric.")
    }
  }

  #########################################
  # Check for reso, dist, prob
  if(!inherits(reso, c("numeric", "integer"))){
    stop("'reso' must be a numeric or integer value")
  } else if(!inherits(dist, c("numeric", "integer"))){
    stop("'dist' must be a numeric or integer value")
  } else if(!inherits(prob, c("numeric", "integer"))){
    stop("'prob' must be a numeric or integer value")
  } else if(prob < 0 || prob > 1){
    stop("'prob' must be between 0 and 1")
  }

  #########################################
  # Check for thr
  if(!is.null(thr)){
    if(!inherits(thr, c("numeric", "integer"))){
      if(!is.null(thr)){
        stop("'thr' must be NULL, a numeric or an integer threshold value.")
      }
    }
  }

  #########################################
  # Check for summed
  if(!is.logical(summed)){
    stop("'summed' must be a logical (TRUE or FALSE).")
  }

  if(summed){
    if(is.null(thr)){
      message("The 'summed' argument is not used when thr = NULL.")
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

  cmd <- c(cmd, "--uselinkset", linkset, "--usegraph", graph, "--interp")
  cmd <- c(cmd,
           paste0("var=", var),
           paste0("d=", dist),
           paste0("p=", prob),
           paste0("name=", name),
           paste0("resol=", reso))

  if(!is.null(thr)){
    cmd <- c(cmd, paste0("multi=", thr))
    if(summed){
      cmd <- c(cmd, "ag=sum")
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
      message(paste0("A ", name, " interpolated raster file has been created."))
    }
  } else {
    message(paste0("A ", name, " interpolated raster file has been created."))
  }

}
