#' Add a point set to the Graphab project
#'
#' @description The function adds a spatial point set to the Graphab project.
#' These points will define a new habitat type with 0-capacity at first.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is.
#' @param name A character string indicating the name of the new habitat type
#' created by adding this point set to the project.
#' @param pointset Can be either:\itemize{
#' \item{A character string indicating the path (absolute or relative) to a
#' geopackage point layer with its extension (.gpkg).}
#' \item{A character string indicating the path to a .csv file with three
#' columns: ID, x and y, respectively indicating the point ID, longitude
#' and latitude.}
#' \item{A data.frame with three columns:
#' ID, x and y, respectively indicating the point ID, longitude and latitude.}
#' }
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @param parallel.java An integer indicating how many computer cores are used
#' to run the .jar file. By default, \code{parallel.java = NULL}, and java sets
#' it according to local settings.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @param linkset Deprecated parameter from graph4lg < 2.0
#' @param id Deprecated parameter from graph4lg < 2.0
#' @param return_val Deprecated parameter from graph4lg < 2.0
#' @details Point coordinates must be in the same coordinate reference system
#' as the habitat patches (and initial raster layer). See more information in
#' Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' graphab_pointset(proj_name = "graphab_example",
#'                  name = "pts_pop",
#'                  pointset = "pts_pop.gpkg")
#' }


graphab_pointset <- function(proj_name,
                             linkset = NULL, # deprecated
                             name,
                             pointset,
                             id = NULL, # deprecated
                             return_val = NULL, # deprecated
                             proj_path = NULL,
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

  # Check for deprecated parameters
  if(!is.null(linkset)){
    stop(paste0("Argument 'linkset' is deprecated and not used anymore ",
                "in graph4lg >= 2.0. Please use graph4lg <= 1.8 if you ",
                "need to use graphab_pointset() to get metric values ",
                "from the nearest patches or distances to patches. ",
                "Note that for this purpose, a new linkset can be created ",
                "between the point set and a given habitat type using ",
                "graphab_link()."))
  } else if(!is.null(id)){
    stop(paste0("Argument 'id' is deprecated and not used anymore ",
                "in graph4lg >= 2.0. Please use graph4lg <= 1.8 if you ",
                "need to use graphab_pointset() to get metric values ",
                "from the nearest patches or distances to patches. ",
                "Note that for this purpose, a new linkset can be created ",
                "between the point set and a given habitat type using ",
                "graphab_link()."))
  } else if(!is.null(return_val)){
    stop(paste0("Argument 'return_val' is deprecated and not used anymore ",
                "in graph4lg >= 2.0. Please use graph4lg <= 1.8 if you ",
                "need to use graphab_pointset() to get metric values ",
                "from the nearest patches or distances to patches. ",
                "Note that for this purpose, a new linkset can be created ",
                "between the point set and a given habitat type using ",
                "graphab_link()."))
  }

  ###############################
  # Check for pointset type and generate a .gpkg if needed

  # If character string : geopackage, csv
  if(inherits(pointset, "character")){

    # If geopackage in project root directory
    if(stringr::str_sub(pointset, -5, -1) == ".gpkg"){

      if(file.exists(pointset)){

        pts_gpkg <- normalizePath(pointset)
        p_type <- "gpkg"

      } else {
        stop(paste0("'pointset' geopackage layer '", pointset,
                    "' does not exist or is not found."))
      }

    } else if (stringr::str_sub(pointset, -4, -1) == ".csv"){

      # If csv
      if(file.exists(pointset)){

        pts <- utils::read.csv(file = pointset)
        if(!all(c("ID", 'x', 'y') %in% colnames(pts))){
          stop("The columns of pts must include 'ID', 'x' and 'y'.")
        } else {
          pts <- pts[, c('ID', 'x', 'y')]
        }
        p_type <- "csv"

      } else {
        stop(paste0("'pointset' .csv file '", pointset,
                    "' does not exist or is not found."))
      }
    }

  } else if(inherits(pointset, "SpatialPointsDataFrame")){
    # If SPDF
    stop(paste0("'pointset' cannot be a SpatialPointsDataFrame in ",
                "graph4lg >= 2.0."))

  } else if (inherits(pointset, "data.frame")){

    # If data.frame
    if(all(c("ID", 'x', 'y') %in% colnames(pointset))){
      p_type <- "df"
      pts <- pointset[, c('ID', 'x', 'y')]
    } else {
      stop("The columns of pts must include 'ID', 'x' and 'y'.")
    }

  } else {

    # else ERROR
    stop("'pointset' must be either a path to a geopackage layer, to a .csv
          file, or a data.frame.")
  }


  if(p_type %in% c("csv", "df")){

    # #### Get project CRS
    # # Copy the .xml file as a .txt file in temp files to open it
    # xml <- tempfile(pattern = ".txt")
    # file.copy(from = proj_end_path,
    #           to = xml)
    # file_data <- utils::read.table(xml)
    # line_crs <- file_data[which(stringr::str_detect(string = file_data[, 1],
    #                                                 pattern = "</wktCRS>")), 1]
    # project_crs <- as.numeric(stringr::str_sub(line_crs, 34, -18))

    # Create a spatial point layer and write to gpkg
    xy <- pts[, c('x', 'y')]
    # Create a list with the spatial points coordinates
    mxy <- as.matrix(xy)
    # Create a list of point objects
    list_pts <- list()
    for(i in 1:nrow(xy)){
      list_pts[[i]] <- sf::st_point(mxy[i, ])
    }
    # Create the point layer
    pts_geom <- sf::st_sfc(list_pts)
                           #crs = project_crs)
    pts_layer <- sf::st_sf(pts_geom,
                           pts)

    pts_gpkg <- tempfile(fileext = ".gpkg")

    # Export pts_layer as .gpkg layer
    sf::st_write(obj = pts_layer,
                 dsn = pts_gpkg,
                 layer = "pts_layer",
                 quiet = TRUE)

  }

  #########################################
  # Check for parallel.java
  if(!is.null(parallel.java)){
    if(!inherits(parallel.java, c("numeric", "integer"))){
      stop("'parallel.java' must be a numeric or integer value.")
    }
  }

  # Check for Graphab
  gr <- get_graphab(res = FALSE, return = TRUE)

  if(gr == 1){
    message("Graphab has been downloaded")
  }

  # Get java path
  java.path <- Sys.which("java")

  # Get graphab path
  version <- "graphab-3.0.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg2_jar/", version)

  # Command line
  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

  if(!is.null(parallel.java)){
    cmd <- c(cmd, "-proc ", as.character(parallel.java))
  }

  cmd <- c(cmd,
           "--project", proj_end_path_cmd,
           "--dataset", paste0("name=", name),
           paste0("file=", pts_gpkg))

  if(!is.null(alloc_ram)){
    if(inherits(alloc_ram, c("integer", "numeric"))){
      cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
    } else {
      stop("'alloc_ram' must be a numeric or an integer")
    }
  }

  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)

  ## Check whether an error occurred
  if(length(rs) == 1){
    if(rs == 1){
      message("An error occurred")
    }
  }

  ## Check whether the habitat exists
  if(check_graphab_object(proj_path = proj_end_path,
                          object_type = "habitat",
                          name = name)){
    message(paste0("Pointset '", name, "' has been created in the project '",
                   proj_name, "'."))
  } else {
    message("An error occurred")
  }

}
