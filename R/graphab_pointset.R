#' Add a point set to the Graphab project
#'
#' @description The function adds a spatial point set to the Graphab project,
#' allowing users to identify closest habitat patch from each point and
#' get corresponding connectivity metrics.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is.
#' @param linkset A character string indicating the name of the link set used.
#' The link set is here used to get the defined cost values and compute the
#' distance from the point to the patches. Link sets can be created
#' with \code{\link{graphab_link}}.
#' @param pointset Can be either;\itemize{
#' \item{A character string indicating the path (absolute or relative) to a
#' shapefile point layer}
#' \item{A character string indicating the path to a .csv file with three
#' columns: ID, x and y, respectively indicating the point ID, longitude
#' and latitude.}
#' \item{A data.frame with three columns:
#' ID, x and y, respectively indicating the point ID, longitude and latitude.}
#' \item{A SpatialPointsDataFrame}
#' }
#' The point ID column must be 'ID' by default but can also be specified
#' by the \code{id} argument in all three cases.
#' @param id A character string indicating the name of the column in either
#' the .csv table, data.frame or attribute table, corresponding to the ID
#' of the points. By default, it should be 'ID'. This column is used for naming
#' the points when returning the output.
#' @param return_val Logical (default=TRUE) indicating whether the metrics
#' associated with closest habitat patches from the points are returned to
#' users.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @return If \code{return_val=TRUE}, the function returns a \code{data.frame}
#' with the properties of the nearest patch to every point in the point set,
#' as well as the distance from each point to the nearest patch.
#' @details Point coordinates must be in the same coordinate reference system
#' as the habitat patches (and initial raster layer). See more information in
#' Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_pointset(proj_name = "grphb_ex",
#'                graph = "graph",
#'                pointset = "pts.shp")
#' }


graphab_pointset <- function(proj_name,
                             linkset,
                             pointset,
                             id = "ID",
                             return_val = TRUE,
                             proj_path = NULL,
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

  # Check for linkset
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (length(list.files(path = paste0(proj_path, "/", proj_name),
                               pattern = "-links.csv")) == 0){
    stop("There is not any linkset in the project you refer to.
         Please use graphab_link() before.")
  } else if (!(paste0(linkset, "-links.csv") %in%
               list.files(path = paste0(proj_path, "/", proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  # Check for id
  if(!inherits(id, "character")){
    stop(paste0("'id' argument must be a character string specifying the",
                " name of the point ID column."))

  }

  # Check for pointset

  # If character string : shapefile, csv
  if(inherits(pointset, "character")){

    # If shp in project directory
    if(stringr::str_sub(pointset, -4, -1) == ".shp"){

      if(file.exists(pointset)){

        pts_shape <- normalizePath(pointset)
        p_type <- "shp"

        layer_shp <- stringr::str_sub(basename(pts_shape), 1, -5)
        dir_shp <- dirname(pts_shape)

        # Open the attribute table of the shapefile
        pts_table <- foreign::read.dbf(file = paste0(dir_shp, "/",
                                                     layer_shp,
                                                     ".dbf"))

        # Check for an id column in the layer attributes
        if(!(id %in% colnames(pts_table))){
          # Return an error if id is not a column of the attribute table
          stop(paste0("'pointset' shapefile layer must include an attribute",
                      " named ", id, "."))
        }


      } else {
        stop(paste0("'pointset' shapefile layer '", pointset,
                    "' does not exist."))
      }

    } else if (stringr::str_sub(pointset, -4, -1) == ".csv"){

      # If csv

      if(file.exists(pointset)){

        pts <- utils::read.csv(file = pointset)

        if(!all(c(id, 'x', 'y') %in% colnames(pts))){
          stop(paste0("The columns of pts must include '", id, "', 'x' and 'y'"))
        } else {
          pts <- pts[, c(id, 'x', 'y')]
        }

        p_type <- "csv"

        # Create a spatial points data.frame
        pts_layer <- sp::SpatialPointsDataFrame(coords = pts[ , c('x', 'y')],
                                                data = pts)

      }

    }

  } else if(inherits(pointset, "SpatialPointsDataFrame")){

    # If SPDF

    p_type <- "spdf"

    # Check for an id column in the spdf attributes
    if(!(id %in% colnames(pointset@data))){
      # Return an error if id is not a column of the attribute table
      stop(paste0("'pointset' SpatialPointsDataFrame must include an attribute",
                  " named ", id, "."))
    }

    pts_layer <- pointset

  } else if (inherits(pointset, "data.frame")){

    # If data.frame

    if(all(c(id, 'x', 'y') %in% colnames(pointset))){
      p_type <- "df"
      pts <- pointset[, c(id, 'x', 'y')]

      pts_layer <- sp::SpatialPointsDataFrame(coords = pts[ , c('x', 'y')],
                                              data = pts)

    } else {
      stop(paste0("The columns of pts must include '", id, "', 'x' and 'y'"))
    }

  } else {

    # else ERROR
    stop("'pointset' must be either a path to a shapefile layer, to a .csv file,
          a SpatialPointsDataFrame object or a data.frame.")
  }


  if(p_type %in% c("csv", "spdf", "df")){

    pts_shape <- tempfile(fileext = ".shp")

    layer_shp <- stringr::str_sub(basename(pts_shape), 1, -5)
    dir_shp <- dirname(pts_shape)

    # Export the shapefile layer
    sf::st_write(obj = sf::st_as_sf(pts_layer),
                 dsn = dir_shp,
                 layer = layer_shp,
                 driver = "ESRI Shapefile")

  }

  # Check for return_val
  if(!is.logical(return_val)){
    stop("'return_val' must be a logical (TRUE or FALSE).")
  }

  # Check for Graphab
  gr <- get_graphab(res = FALSE, return = TRUE)

  if(gr == 1){
    message("Graphab has been downloaded")
  }

  # Get java path
  java.path <- Sys.which("java")

  # Get graphab path
  version <- "graphab-2.8.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)

  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
           "--project", proj_end_path,
           "--uselinkset", linkset,
           "--pointset", pts_shape,
           paste0("id=", id))

  if(!is.null(alloc_ram)){
    if(inherits(alloc_ram, c("integer", "numeric"))){
      cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
    } else {
      stop("'alloc_ram' must be a numeric or an integer")
    }
  }

  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)

  name_pts <- paste0("Exo-",
                     layer_shp,
                     "_", linkset)

  if(return_val){

    df_pts <- foreign::read.dbf(file = paste0(proj_path, "/",
                                              proj_name, "/",
                                              name_pts, ".dbf"))
    df_pts <- df_pts[, which(colnames(df_pts) %in% c("Id", "IdPatch", "Cost"))]

    patches <- utils::read.csv(file = paste0(proj_path, "/",
                                             proj_name, "/patches.csv"))

    df_pts <- merge(df_pts, patches, by.x = "IdPatch", by.y = "Id")

    colnames(df_pts)[which(colnames(df_pts) == "IdPatch")] <- "Nearest_patch_ID"
    colnames(df_pts)[which(colnames(df_pts) == "Id")] <- "Point_ID"
    colnames(df_pts)[which(colnames(df_pts) == "Cost")] <- "Dist_to_patch"


    # if(p_type %in% c("csv", "spdf", "df")){
    #   if(stringr::str_detect(dir_shp, pattern = "Temp")){
    #     file.remove(list.files(path = dir_shp, pattern = layer_shp))
    #   }
    # }


    if(length(rs) == 1){
      if(rs == 1){
        message("An error occurred")
      } else {
        if(file.exists(paste0(proj_path, "/",
                              proj_name, "/",
                              name_pts, ".shp"))){
          message(paste0("Point set '", name_pts ,"' has been added to the project ",
                         proj_name))
        } else {
          message("The point set import did not succeed.")
        }
      }
    } else {
      if(file.exists(paste0(proj_path, "/",
                            proj_name, "/",
                            name_pts, ".shp"))){
        message(paste0("Point set '", name_pts ,"' has been added to the project ",
                       proj_name))
      } else {
        message("The point set import did not succeed.")
      }
    }

    return(df_pts)
  } else {

    # if(p_type %in% c("csv", "spdf", "df")){
    #   if(stringr::str_detect(dir_shp, pattern = "Temp")){
    #     file.remove(list.files(path = dir_shp, pattern = layer_shp))
    #   }
    # }

    if(length(rs) == 1){
      if(rs == 1){
        message("An error occurred")
      } else {
        if(file.exists(paste0(proj_path, "/",
                              proj_name, "/",
                              name_pts, ".shp"))){
          message(paste0("Point set '", name_pts ,"' has been added to the project ",
                         proj_name))
        } else {
          message("The point set import did not succeed.")
        }
      }
    } else {
      if(file.exists(paste0(proj_path, "/",
                            proj_name, "/",
                            name_pts, ".shp"))){
        message(paste0("Point set '", name_pts ,"' has been added to the project ",
                       proj_name))
      } else {
        message("The point set import did not succeed.")
      }
    }


  }



}




