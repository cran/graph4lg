#' Create a Graphab project
#'
#' @description The function creates a Graphab project from a raster file on
#' which habitat patches can be delimited.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml will be created.
#' @param raster A character string indicating the name of the .tif raster file
#' or of its path. If the path is not specified, the raster must be present in
#' the current working directory. Raster cell values must be in INT2S encoding.
#' @param nodata (optional, default=NULL) An integer or numeric value
#' specifying the code in the raster file associated with nodata value
#' (often corresponding to peripheric cells)
#' @param parallel.java An integer indicating how many computer cores are used
#' to run the .jar file. By default, \code{parallel.java = NULL}, and java sets
#' it according to local settings.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @details A habitat patch consists of the central pixel with its eight
#' neighbors if they are of the same value (8-connexity) and the path
#' geometry is not simplified. See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' @export
#' @author P. Savary, T. Rudolph
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' proj_name <- "graphab_example"
#' raster <- "rast_simul50.tif"
#' graphab_project(proj_name = proj_name,
#'                raster = raster)
#' }


graphab_project <- function(proj_name,
                            raster,
                            nodata = NULL,
                            parallel.java = NULL,
                            alloc_ram = NULL,
                            proj_path = NULL){


  ## Store the current directory path
  #wd1 <- getwd()

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

  #######################################
  # Add '' to proj_path for cases with spaces in paths
  if(all(stringr::str_sub(proj_path, 1, 1) != "'",
         stringr::str_sub(proj_path, 1, 1) != "'",
         stringr::str_detect(string = proj_path,
                             pattern = " "))){
    proj_path_cmd <- paste0("'", proj_path, "'")
  } else {
    proj_path_cmd <- proj_path
  }

  #########################################
  # Check for proj_name class
  if(!inherits(proj_name, "character")){
    stop("'proj_name' must be a character string")
  }

  #########################################
  # Check for raster class
  if(!inherits(raster, "character")){
    stop("'raster' must be a character string")
    # Is raster a tif file?
  } else if (stringr::str_sub(raster, -4, -1) != ".tif"){
    stop(paste0(raster, " must be a .tif raster file."))
    # Is raster an existing file
  } else if (!(file.exists(normalizePath(raster, mustWork = FALSE)))){
    stop(paste0(normalizePath(raster, mustWork = FALSE),
                " must be an existing .tif raster file."))
  }

  #######################################
  # Add '' to raster for cases with spaces in paths not transformed before
  if(all(stringr::str_sub(raster, 1, 1) != "'",
         stringr::str_sub(raster, 1, 1) != "'",
         stringr::str_detect(string = raster,
                             pattern = " "))){
    raster_cmd <- paste0("'", raster, "'")
  } else {
    raster_cmd <- raster
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

  cmd <- c(cmd, "--create", proj_name, raster_cmd)

  if(!is.null(nodata)){
    cmd <- c(cmd, paste0("nodata=", nodata))
  }

  cmd <- c(cmd, paste0("dir=", proj_path_cmd))

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
      if(file.exists(paste0(proj_path, "/", proj_name, "/", proj_name, ".xml"))){
        message(paste0("Graphab project ", proj_name," has been created in directory: ",
                       proj_path))
      } else {
        message("The project creation did not succeed.")
      }
    }
  } else {

    if(file.exists(paste0(proj_path, "/", proj_name, "/", proj_name, ".xml"))){
      message(paste0("Graphab project ", proj_name," has been created in directory: ",
                     proj_path))
    } else {
      message("The project creation did not succeed.")
    }
  }

}
