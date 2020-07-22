#' Create a Graphab project
#'
#' @description The function creates a Graphab project from a raster file on
#' which habitat patches can be delimited.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml will be created.
#' @param raster A character string indicating the name of the .tif raster file
#' or of its path. If the path is not specified, the raster must present in the
#' current working directory. Raster cell values must be in INT2S encoding.
#' @param habitat An integer or numeric value or vector indicating the
#' code.s (cell value.s) of the habitat cells in the raster file.
#' @param minarea (optional, default=0) An integer or numeric value specifiying
#' the minimum area in hectares for a habitat patch size to become a graph node.
#' @param nodata (optional, default=NULL) An integer or numeric value
#' specifying the code in the raster file associated with nodata value
#' (often corresponding to peripheric cells)
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
#' geometry is not simplified. See more information in Graphab 2.4 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.4-en.pdf}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' proj_name <- "grphb_ex"
#' raster <- "rast_ex.tif"
#' habitat <- 5
#' graphab_project(proj_name = proj_name,
#'                raster = raster,
#'                habitat = habitat)
#' }


graphab_project <- function(proj_name,
                            raster,
                            habitat,
                            minarea = 0,
                            nodata = NULL,
                            alloc_ram = NULL,
                            proj_path = NULL){


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
  }

  #########################################
  # Check for raster class
  if(!inherits(raster, "character")){
    stop("'raster' must be a character string")
    # Is raster present in the working directory?
  } else if (!(raster %in% list.files(pattern = ".tif"))){
    stop("'raster' must be a .tif raster file present in the working directory")
  }

  #########################################
  # Check for habitat class
  if(!inherits(habitat, c("numeric", "integer"))){
    stop("'habitat' must be an integer indicating the habitat code")
  }

  #########################################
  # Check for minarea class
  if(!inherits(minarea, c("numeric", "integer"))){
    stop("'minarea' must be an integer indicating minimum patch size")
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
  version <- "graphab-2.4.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)


  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab, "--create", proj_name, raster,
           paste0("habitat=", paste(habitat, collapse = ",")),
           paste0("minarea=", minarea), "con8")

  if(!is.null(nodata)){
    cmd <- c(cmd, paste0("nodata=", nodata))
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
    }
  } else if (!is.null(proj_path)){
    message(paste0("Graphab project ",proj_name," has been created in directory: ",
                   proj_path))
  } else {
    message(paste0("Graphab project ",proj_name," has been created in directory: ",
                   getwd()))
  }

  #########################################
  if(chg == 1){
    setwd(dir = wd1)
  }

}
