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
#' and latitude}
#' \item{A data.frame with three columns:
#' ID, x and y, respectively indicating the point ID, longitude and latitude.}
#' \item{A SpatialPointsDataFrame}
#' }
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
#' Graphab 2.4 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.4-en.pdf}
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
                             return_val = TRUE,
                             proj_path = NULL,
                             alloc_ram = NULL){

  # Check for project directory path
  if(!is.null(proj_path)){
    chg <- 1
    wd1 <- getwd()
    setwd(dir = proj_path)
  } else {
    chg <- 0
    proj_path <- getwd()
  }

  # Check for proj_name class
  if(!inherits(proj_name, "character")){
    stop("'proj_name' must be a character string")
  } else if (!(paste0(proj_name, ".xml") %in% list.files(path = paste0("./", proj_name)))){
    stop("The project you refer to does not exist.
         Please use graphab_project() before.")
  }

  proj_end_path <- paste0(proj_name, "/", proj_name, ".xml")

  # Check for linkset
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (length(list.files(path = paste0("./", proj_name), pattern = "-links.csv")) == 0){
    stop("There is not any linkset in the project you refer to.
         Please use graphab_link() before.")
  } else if (!(paste0(linkset, "-links.csv") %in% list.files(path = paste0("./", proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  # Check for pointset

  # If character string : shapefile, csv
  if(inherits(pointset, "character")){

    # If shp in project directory
    if(stringr::str_sub(pointset, -4, -1) == ".shp"){

      if(file.exists(pointset)){

        pts_shape <- pointset
        p_type <- "shp"

        path_s <- strsplit(pts_shape, "/|\\\\")[[1]]

        n <- length(path_s)
        layer_shp <- stringr::str_sub(path_s[n], 1, -5)
        dir_shp <- stringr::str_sub(pts_shape, 1, -(nchar(layer_shp)+6))

      } else {

        stop(paste0("'pointset' shapefile layer '", pointset,
                    "' does not exist."))
      }

    } else if (stringr::str_sub(pointset, -4, -1) == ".csv"){

      # If csv

      if(file.exists(pointset)){

        pts <- utils::read.csv(file = pointset)

        if(!all(c('ID', 'x', 'y') %in% colnames(pts))){
          stop("The columns of pts must include 'ID', 'x' and 'y'")
        } else {
          pts <- pts[, c('ID', 'x', 'y')]
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


  } else if (inherits(pointset, "data.frame")){

    # If data.frame

    if(all(c('ID', 'x', 'y') %in% colnames(pointset))){
      p_type <- "df"
      pts <- pointset[, c('ID', 'x', 'y')]

      pts_layer <- sp::SpatialPointsDataFrame(coords = pts[ , c('x', 'y')],
                                              data = pts)

    } else {
      stop("The columns of pts must include 'ID', 'x' and 'y'")
    }

  } else {

    # else ERROR

    stop("'pointset' must be either a path to a shapefile layer, to a .csv file,
          a SpatialPointsDataFrame object or a data.frame.")
  }


  if(p_type %in% c("csv", "spdf", "df")){

    pts_shape <- tempfile(fileext = ".shp")

    path_s <- strsplit(pts_shape, "/|\\\\")[[1]]

    n <- length(path_s)
    layer_shp <- stringr::str_sub(path_s[n], 1, -5)
    dir_shp <- stringr::str_sub(pts_shape, 1, -(nchar(layer_shp)+6))

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
  version <- "graphab-2.4.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)

  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
           "--project", proj_end_path,
           "--uselinkset", linkset,
           "--pointset", pts_shape)

  if(!is.null(alloc_ram)){
    if(inherits(alloc_ram, c("integer", "numeric"))){
      cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
    } else {
      stop("'alloc_ram' must be a numeric or an integer")
    }
  }

  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)

  if(return_val){

    name_pts <- paste0("Exo-",
                       layer_shp,
                       "_", linkset)

    df_pts <- foreign::read.dbf(file = paste0(proj_name, "/", name_pts, ".dbf"))
    df_pts <- df_pts[, which(colnames(df_pts) %in% c("Id", "IdPatch", "Cost"))]

    patches <- utils::read.csv(file = paste0(proj_name, "/patches.csv"))

    df_pts <- merge(df_pts, patches, by.x = "IdPatch", by.y = "Id")

    colnames(df_pts)[which(colnames(df_pts) == "IdPatch")] <- "Nearest_patch_ID"
    colnames(df_pts)[which(colnames(df_pts) == "Id")] <- "Point_ID"
    colnames(df_pts)[which(colnames(df_pts) == "Cost")] <- "Dist_to_patch"

    return(df_pts)
  }

  if(p_type %in% c("csv", "spdf", "df")){
    if(stringr::str_detect(dir_shp, pattern = "Temp")){
      file.remove(list.files(path = dir_shp, pattern = layer_shp))
    }
  }

  if(chg == 1){
    setwd(dir = wd1)
  }

  if(length(rs) == 1){
    if(rs == 1){
      message("An error occurred")
    } else {
      message(paste0("Point set has been added to the project ",
                     proj_name))
    }
  } else {
    message(paste0("Point set has been added to the project ",
                   proj_name))
  }

}




