#' Compute Euclidean geographic distances between points
#'
#' @description The function computes Euclidean geographic distance between
#' points given their spatial coordinates either in a metric projected
#' Coordinate Reference System or in a polar coordinates system.
#'
#' @param data An object of class :\itemize{
#' \item \code{data.frame} with 3 columns: 2 columns with the point spatial
#' coordinates and another column with point IDs
#' \item \code{SpatialPointsDataFrame} }
#' @param ID (if \code{data} is of class \code{data.frame}) A character string
#' indicating the name of the column of \code{data} with the point IDs
#' @param x (if \code{data} is of class \code{data.frame}) A character string
#' indicating the name of the column of \code{data} with the point longitude
#' @param y (if \code{data} is of class \code{data.frame}) A character string
#' indicating the name of the column of \code{data} with the point latitude
#' @param crds_type A character string indicating the type of coordinate
#' reference system:\itemize{
#' \item{'proj' (default): a projected coordinate reference system}
#' \item{'polar': a polar coordinate reference system, such as WGS84}
#' }
#' @param gc_formula A character string indicating the formula used to compute
#' the Great Circle distance:\itemize{
#' \item{'vicenty'(default): Vincenty inverse formula for ellipsoids}
#' \item{'slc': Spherical Law of Cosines}
#' \item{'hvs': Harversine formula}
#' }
#' @return A pairwise matrix of geographic distances between points in meters
#' @export
#' @author P. Savary
#' @details When a projected coordinate reference system is used, it calculates
#' classical Euclidean geographic distance between two points using
#' Pythagora's theorem. When a polar coordinate reference system is used, it
#' calculates the Great circle distance between points using different methods.
#' Unless \code{method = "polar"}, when \code{data} is a \code{data.frame},
#' it assumes projected coordinates by default.
#' @examples
#' # Projected CRS
#' data(pts_pop_simul)
#' mat_dist <- mat_geo_dist(data=pts_pop_simul,
#'              ID = "ID",
#'              x = "x",
#'              y = "y")
#'
#' #Polar CRS
#' city_us <- data.frame(name = c("New York City", "Chicago",
#'                                "Los Angeles", "Atlanta"),
#'                       lat  = c(40.75170,  41.87440,
#'                                34.05420,  33.75280),
#'                       lon  = c(-73.99420, -87.63940,
#'                               -118.24100, -84.39360))
#' mat_geo_us <- mat_geo_dist(data = city_us,
#'                            ID = "name", x = "lon", y = "lat",
#'                            crds_type = "polar")



mat_geo_dist <- function(data,
                         ID = NULL,
                         x = NULL,
                         y = NULL,
                         crds_type = "proj",
                         gc_formula = "vicenty"){


  ### Check for crds_type
  if(!inherits(crds_type, "character")){
    stop("'crds_type' must be a character string")
  } else if(!(crds_type %in% c("proj", "polar"))){
    stop("'crds_type' must be 'proj' or 'polar'.")
  }

  ### Check for gc_formula
  if(!inherits(gc_formula, "character")){
    stop("'gc_formula' must be a character string")
  } else if(!(gc_formula %in% c("vicenty", "slc", "hvs"))){
    stop("'gc_formula' must be 'vicenty', 'slc' or 'hvs'.")
  }


  # If 'data' is a Spatial Points data.frame
  if(inherits(data, "SpatialPointsDataFrame")){

    # Check whether locations are not duplicated (and display a warning message
    # in such a case)
    if(any(duplicated(data@coords))){
      warning("At least 1 point location appears twice in
              the 'SpatialPointsDataFrame'.")
    }

    # Check whether the points have an ID
    if(is.null(ID) ) {
      stop("You have to specify the name of the ID column of the points in the
           attribute table, as an input 'id'")
    }

    # Check whether the argument 'x' is specified
    # because in this case it is useless
    if(!is.null(x) ) {
      warning("Unused argument 'x'")
    }

    # Check whether the argument 'y' is specified
    # because in this case it is useless
    if(!is.null(y) ) {
      warning("Unused argument 'y'")
    }


    #####################################################################
    #####################################################################

    # Check whether the data have projected coordinates
    if(stringr::str_sub(raster::crs(data), 7, 13) == "longlat"){

      if(crds_type == "proj"){
        stop("Your SpatialPointsDataFrame has polar coordinates
             you cannot use 'crds_type = 'proj'")
      } else {

        data <- cbind(data.frame(ID = data$ID), data.frame(data@coords))
        colnames(data) <- c("ID", "x", "y")

        df_lk <- expand.grid(data[, ID], data[, ID])
        df_lk <- merge(df_lk, data, by.x = "Var1", by.y = "ID")
        df_lk <- merge(df_lk, data, by.x = "Var2", by.y = "ID")
        colnames(df_lk) <- c("ID2", "ID1", "long1", "lat1", "long2", "lat2")
        df_lk <- df_lk[, c("ID1", "long1", "lat1", "ID2", "long2", "lat2")]

        df_lk$dist <- NA

        for(i in 1:nrow(df_lk)){

          df_lk[i, 'dist'] <- dist_great_circle(long1 = df_lk[i, "long1"],
                                                lat1 = df_lk[i, "lat1"],
                                                long2 = df_lk[i, "long2"],
                                                lat2 = df_lk[i, "lat2"],
                                                method = gc_formula)
        }

        mat <- graph4lg::df_to_pw_mat(data = df_lk, from = "ID1", to = "ID2", value = "dist")


      }

    } else {

      # Projected coordinates

      if(crds_type == "polar"){
        stop("Your SpatialPointsDataFrame has projected coordinates
             you cannot use 'crds_type = 'polar'")

      } else {

        # Check whether the data have projected coordinates in a common CRS
        if(!(stringr::str_sub(raster::crs(data), 7, 9) %in% c("lcc",
                                                              "utm", "mer"))){
          message("The CRS of your SpatialPointsDataFrame is not common,
              ensure it has projected coordinates.")
        }


        # Get the coordinates of the points
        coords <- data@coords

        # Compute the Euclidean distances between the points
        mat <- as.matrix(stats::dist(coords,
                                     method = "euclidean",
                                     diag = TRUE,
                                     upper = TRUE))
        # Give names to the rows and columns of the distance matrix
        name <- data@data[, ID]
        row.names(mat) <- colnames(mat) <- name


      }

    }

  } else if (inherits(data, "data.frame")){


    # Check whether 'x' is specified
    if(is.null(x) ) {
      stop("You must specify the name of the 'x' column,
             as an input 'x'")
    }

    # Check whether 'y' is specified
    if(is.null(y) ) {
      stop("You must specify the name of the 'y' column,
             as an input 'y'")
    }

    # Check whether the 'ID' column is specified
    if(is.null(ID) ) {
      stop("You have to specify the name of the ID column of the points in the
           data.frame, as an input 'ID'")
    }


    if(crds_type == "proj"){

      message("Coordinates were treated as projected coordinates. Check whether
              it is the case.")

      # Get the coordinates of the points
      coords <- data[, c(y, x)]

      # Compute the Euclidean distances between the points
      mat <- as.matrix(stats::dist(coords,
                                   method = "euclidean",
                                   diag = TRUE,
                                   upper = TRUE))
      # Give names to the rows and columns of the distance matrix
      name <- data[, ID]
      row.names(mat) <- colnames(mat) <- name


    } else {

      message("Coordinates were treated as polar coordinates. Check whether
              it is the case.")

      df_lk <- expand.grid(data[, ID], data[, ID])
      df_lk <- merge(df_lk, data, by.x = "Var1", by.y = ID)
      df_lk <- merge(df_lk, data, by.x = "Var2", by.y = ID)
      colnames(df_lk) <- c("ID2", "ID1", "long1", "lat1", "long2", "lat2")
      df_lk <- df_lk[, c("ID1", "long1", "lat1", "ID2", "long2", "lat2")]

      df_lk$dist <- NA

      for(i in 1:nrow(df_lk)){

        df_lk[i, 'dist'] <- dist_great_circle(long1 = df_lk[i, "long1"],
                                              lat1 = df_lk[i, "lat1"],
                                              long2 = df_lk[i, "long2"],
                                              lat2 = df_lk[i, "lat2"],
                                              method = gc_formula)
      }

      mat <- graph4lg::df_to_pw_mat(data = df_lk, from = "ID1", to = "ID2", value = "dist")

    }

  } else {
    stop("Input 'data' must be of class 'data.frame' or
         'SpatialPointsDataFrame'.")
  }

  return(mat)
}



