#' Compute Euclidean geographic distances between points
#'
#' @description The function computes Euclidean geographic distance between
#' points given their coordinates in a metric projected Coordinates
#' Reference System.
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
#' @return A pairwise matrix of geographic distances between points
#' @export
#' @author P. Savary
#' @details Calculate classical Euclidean geographic distance (in m)
#' between two points using Pythagora's theorem
#' @examples
#' data(pts_pop_simul)
#' mat_dist <- mat_geo_dist(data=pts_pop_simul,
#'              ID = "ID",
#'              x = "x",
#'              y = "y")

###############################################################################

###### mat_geo_dist function ###########################

mat_geo_dist <- function(data, ID = NULL, x = NULL, y = NULL){

  # If 'data' is a Spatial Points data.frame
  if(inherits(data, "SpatialPointsDataFrame")){

    # Check whether locations are not duplicated (and display a warning message
    # in such a case)
    if(any(duplicated(data@coords))){
      warning("At least 1 point location appears twice in
              the 'SpatialPointsDataFrame'.")
    }

    # Check whether the data have projected coordinates
    if(stringr::str_sub(raster::crs(data), 7, 13) == "longlat"){
      stop("Your SpatialPointsDataFrame must have
           projected (metric) coordinates")
    }

    # Check whether the data have projected coordinates in a common CRS
    if(!(stringr::str_sub(raster::crs(data), 7, 9) %in% c("lcc",
                                                          "utm", "mer"))){
      message("The CRS of your SpatialPointsDataFrame is not common,
              ensure it has projected coordinates.")
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
    stop("Input 'data' must be of class 'data.frame' or
         'SpatialPointsDataFrame'.")
  }

  return(mat)
}



