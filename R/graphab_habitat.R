#' Create a habitat type in a Graphab project
#'
#' @description The function creates a new habitat in the Graphab project in
#' two possible ways:\itemize{
#' \item{From the raster file defining the project, based on a code value.}
#' \item{From an external vector file (polygon or point type). NOT YET}
#' }
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param type A character string specifying whether the new habitat type
#' is defined from the raster layer defining the project (`type="raster"`,
#' default) or from an external vector layer (`type="vector"`) of type polygon
#' or point.
#' @param name A character string specifying the name (arbitrary) of the new
#' habitat type in the raster project.
#' @param rast_codes (default=NULL) If `type="raster"`, code value(s) on the
#' categorical raster layer defining the project corresponding to the type of
#' habitat. Must be an integer or a vector of integer values.
#' @param vec_layer A character string indicating the path (absolute or
#' relative) to a vector layer in shapefile (.shp) or geopackage (.gpkg) format.
#' @param vec_capa_field (default=NULL) If `type="vector"`, name of the
#' column of the vector layer attribute table setting the capacity of each
#' habitat patch.
#' @param minarea (optional, default=0) An integer or numeric value specifiying
#' the minimum area in hectares for a habitat patch size to become a graph node.
#' @param maxsize (optional, default=NULL, only if `type="raster"`). An integer
#' or numeric value specifying the maximum side length of the rectangular full
#' extent of each habitat patch in metric units. If this side length exceeds
#' \code{maxsize} m, then several patches are created.
#' @param con8 (optional, default=FALSE) A logical indicating whether a
#' neighborhood of 8 pixels (TRUE) is used for patch definition. By default,
#' \code{con8=4}, corresponding to 4 pixel neighborhood.
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
#' graphab_habitat(proj_name = proj_name,
#'                name = "forest",
#'                type = "raster",
#'                rast_codes = c(1,2))
#' }


graphab_habitat <- function(proj_name,
                            name,
                            type = "raster",
                            rast_codes = NULL,
                            vec_layer = NULL,
                            vec_capa_field = NULL,
                            minarea = 0,
                            maxsize = NULL,
                            con8 = FALSE,
                            parallel.java = NULL,
                            alloc_ram = NULL,
                            proj_path = NULL){


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
  # Check for name
  if(!inherits(name, "character")){
    stop("'name' must be a character string")
  }

  #########################################
  # Check for type
  if(!inherits(type, "character")){
    stop("'type' must be a character string")
  } else if(!any(type %in% c("raster", "vector"))){
    stop(paste0("'type' must be equal to either ",
                "'raster' or 'vector'."))
  }

  #########################################
  # Check for rast_codes
  if(type == "raster"){
    if(!inherits(rast_codes, c("numeric", "integer"))){
      stop(paste0("'rast_codes' must be an integer or a vector of integers",
                  " indicating the habitat code(s)."))
    }

    ## Check that the codes are in the raster
    raster_codes <- get_graphab_raster_codes(proj_name = proj_name,
                             mode = "all",
                             proj_path = proj_path)

    if(!all(rast_codes %in% raster_codes)){
      stop(paste0("All 'rast_codes' values must be integer values",
                  " corresponding to existing values in the project raster."))
    }
  }

  #########################################
  # Check for vec_layer
  if(!is.null(vec_layer)){
    # if not null, check type is not raster
    if(type == "raster"){
      message("'vec_layer' is specified but won't be used, as type = 'raster'.")
    } else {
      # if not null, check this is a character path
      if(!inherits(vec_layer, "character")){
        stop(paste0("'vec_layer' must be a character string providing the path",
                    " to a .shp or .gpkg vector layer."))
      }
    }
  } else if(type == "vector"){
    stop("'vec_layer' must be provided when type = 'vector'.")
  }

  #########################################
  # Check for vec_capa_field
  if(!is.null(vec_capa_field)){
    # if not null, check type is not raster
    if(type == "raster"){
      message("'vec_capa_field' is specified but won't be used, as type = 'raster'.")
    } else {
      # if not null, check this is a character path
      if(!inherits(vec_capa_field, "character")){
        stop(paste0("'vec_capa_field' must be a character string providing the ",
                    " an existing column name in the vector layer attr. table."))
      }
    }
  } else if(type == "vector"){
    stop("'vec_capa_field' must be provided when type = 'vector'.")
  }

  #########################################
  # Check that vector layer exists and vec_capa_field too
  if(type == "vector"){

    # If shp in project directory
    if(stringr::str_sub(vec_layer, -4, -1) == ".shp"){

      if(file.exists(vec_layer)){

        path_shape <- normalizePath(vec_layer)
        layer_shp <- stringr::str_sub(basename(path_shape), 1, -5)
        dir_shp <- dirname(path_shape)

        # Open the attribute table of the shapefile
        attr_table <- suppressWarnings(
            sf::read_sf(dsn = dir_shp,
                        layer = layer_shp,
                        query = paste0("SELECT * FROM ",
                                       layer_shp,
                                       " LIMIT 0"),
                        as_tibble = FALSE))

        # Check for a 'vec_capa_field' column in the layer attributes
        if(!(vec_capa_field %in% colnames(attr_table))){
          # Return an error if id is not a column of the attribute table
          stop(paste0(vec_layer, " shapefile layer must include an attribute",
                      " named ", vec_capa_field, "."))
        }

      } else {
        stop(paste0(vec_layer, " shapefile layer does not exist."))
      }

    } else if (stringr::str_sub(vec_layer, -5, -1) == ".gpkg"){
      # If gpkg
      if(file.exists(vec_layer)){
        attr_table <- suppressWarnings(
            sf::read_sf(vec_layer,
                        query = paste0("SELECT * FROM ",
                                       stringr::str_sub(basename(vec_layer),
                                                        1, -6),
                                       " LIMIT 0"),
                        as_tibble = FALSE))
        if(!(vec_capa_field %in% colnames(attr_table))){
          stop(paste0("The columns of ", vec_layer, " attribute table must",
                      " include ", vec_capa_field, "."))
        }
      } else {
        stop(paste0(vec_layer, " geopackage layer does not exist at the ",
                    "specified path."))
      }
    }

    # Add '' to vec_layer for cases with spaces in paths
    if(all(stringr::str_sub(vec_layer, 1, 1) != "'",
           stringr::str_sub(vec_layer, 1, 1) != "'",
           stringr::str_detect(string = vec_layer,
                               pattern = " "))){
      vec_layer_cmd <- paste0("'", vec_layer, "'")
    } else {
      vec_layer_cmd <- vec_layer
    }


  }

  #########################################
  # Check for minarea class
  if(!inherits(minarea, c("numeric", "integer"))){
    stop("'minarea' must be an integer indicating minimum patch size")
  }

  #########################################
  # Check for maxsize
  if(!is.null(maxsize)){
    if(type == "raster"){
      if(!inherits(maxsize, c("numeric", "integer"))){
        stop(paste0("'maxsize' must be an integer indicating maximum side length",
                    " of the rectangular extent of every habitat patch"))
      }
    } else {
      message("'maxsize' is specified but won't be used, as type = 'vector'.")
    }
  }

  #########################################
  # Check for con8 class
  if(!inherits(con8, c("logical"))){
    stop(paste0("'con8' must be a logical indicating whether a neighboorhood of ",
                " 8 pixels is used for patch definition (if TRUE). ",
                "Default=FALSE: 4 pixel neighboorhood."))
  }

  # Check for con8 class and type vector
  if(all(con8 == TRUE, type == "vector")){
    message(paste0("'con8' is set to TRUE but type = 'vector'",
                   " so 'con8' won't be used."))
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

  if(type == "raster"){

    cmd <- c(cmd, "--habitat",
             paste0("name=", name),
             paste0("codes=", paste(rast_codes, collapse = ",")))

    cmd <- c(cmd, paste0("minarea=", minarea))

    if(!is.null(maxsize)){
      cmd <- c(cmd, paste0("maxsize=", maxsize))
    }

    if(con8){
      cmd <- c(cmd, "con8")
    }

  } else {

    cmd <- c(cmd, "--vhabitat",
             paste0("name=", name),
             paste0("file=", vec_layer_cmd))
    cmd <- c(cmd, paste0("mincapa=", minarea))
    cmd <- c(cmd, paste0("capa=", vec_capa_field))

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
    message(paste0("Habitat type '", name, "' has been created in the project '",
                   proj_name, "'."))

    if(!file.exists(paste0(proj_path, "/",
                          proj_name, "/",
                          name, "/patches.gpkg"))){
      warning(paste0("The new habitat type does not include any patch. ",
                     "Please check the raster and the patch constraints."))
    }

  } else {
    message("An error occurred")
  }

}
