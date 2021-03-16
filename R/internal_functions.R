#' Compute Gini coefficient from a numeric vector
#'
#' @description The function computes Gini coefficient from a numeric vector
#'
#' @param x A numeric vector with positive values
#' @param unbiased A logical value indicating whether the computed coefficient
#' is biased or not. Unbiased value are equal to n/(n-1) times the biased ones.
#' @return A numeric value corresponding to the Gini coefficient of the numeric
#' vector
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' x <- c(10, 2, 5, 15)
#' gini <- gini_coeff(x)


gini_coeff <- function(x,
                       unbiased = TRUE){

  #####
  # Check x
  if(!inherits(x, c("integer", "numeric"))){
    stop("x must be a vector with numeric or integer values.")
  } else if(length(x) == 1){
    stop("x must include more than 1 value.")
  }

  #####
  # Check unbiased
  if(!is.logical(unbiased)){
    stop("unbiased must be either TRUE or FALSE.")
  }

  x <- x[!is.na(x)]
  x_ord <- x[order(x)]
  mean_x <- mean(x_ord)
  n <- length(x_ord)

  i_val <- 1:length(x_ord)

  terms <- c()
  for(i in 1:n){
    terms[i] <- i * (x_ord[i] - mean_x)
  }

  gini <- (2/(n^2 * mean_x)) * sum(terms)

  if(unbiased){
    gini <- (n/(n-1)) * gini
  }

  return(gini)
}

#' Compute the harmonic mean of a numeric vector
#'
#' @description The function computes the harmonic mean of a numeric vector
#'
#' @param x A numeric vector
#' @return A numeric value corresponding to the harmonic mean of the vector
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' x <- c(10, 2, 5, 15)
#' hm <- harm_mean(x)

harm_mean <- function(x){

  #####
  # Check x
  if(!inherits(x, c("integer", "numeric"))){
    stop("x must be a vector with numeric or integer values.")
  } else if(length(x) == 1){
    stop("x must include more than 1 value.")
  }

  hm <- 1/(mean(1/x))

  return(hm)

}




#' Extract patch areas from a categorical raster
#'
#' @description The function extracts patch areas from a categorical raster
#'
#' @param raster A RasterLayer object corresponding to a categorical raster layer
#' @param class An integer value or vector with the value(s) corresponding to
#' the code values of the raster layer within which points will be sampled.
#' are computed.
#' @param edge_size An integer value indicating the width of the edge
#' (in meters) of the raster layer which is ignored during the sampling
#' (default = 0). It prevents from sampling in the margins of the study area.
#' @param neighborhood An integer value indicating which cells are considered
#' adjacent when contiguous patches are delineated (it should be 8
#' (default, Queen's case) or 4 (Rook's case)). This parameter is ignored
#' when \code{by_patch = FALSE}.
#' @param surf_min An integer value indicating the minimum surface of a patch
#' considered for the sampling in number of raster cells. This parameter is used
#' whatever the \code{by_patch} argument is. Default is 0.
#' @return A data.frame with the areas of the patches
#' @export
#' @keywords internal
#' @author P. Savary

patch_areas <- function(raster,
                        class,
                        edge_size = 0,
                        neighborhood = 8,
                        surf_min = 0){

  # ____________________________
  # ____________________________
  # ------- Check the function arguments

  # raster
  if(!inherits(raster, c("RasterLayer"))){
    stop("'RasterLayer' must be either a 'RasterLayer' object.")
  }

  # edge_size
  if(!inherits(edge_size, c("integer", "numeric"))){
    stop("'edge_size' must be either a 'numeric' or 'integer' value.")
  }
  edge_size <- round(edge_size, digits = 0)

  # neighborhood
  if(!inherits(neighborhood, c("integer", "numeric"))){
    stop("'neighborhood' must be either a 'numeric' or 'integer' value.")
  } else if(!(neighborhood %in% c(4, 8))){
    stop("'neighborhood' must be equal to either 4 or 8.")
  }

  # surf_min
  if(!inherits(surf_min, c("integer", "numeric"))){
    stop("'surf_min' must be either a 'numeric' or 'integer' value.")
  }
  surf_min <- round(surf_min, digits = 0)


  # ____________________________
  # ____________________________
  # ------ Edge definition

  # Get the extent of the raster
  extent_r <- raster::extent(raster)

  if(edge_size == 0){

    # If edge_size is 0, then sample in the whole raster
    polyg_sample <- methods::as(extent_r, "SpatialPolygons")

  } else {

    # If edge_size != 0, then sample in the raster without its edges
    extent_r@xmin <- extent_r@xmin + edge_size
    extent_r@ymin <- extent_r@ymin + edge_size

    extent_r@xmax <- extent_r@xmax - edge_size
    extent_r@ymax <- extent_r@ymax - edge_size

    polyg_sample <- methods::as(extent_r, "SpatialPolygons")

  }

  # Define the CRS of the polygon in which we sample
  raster::crs(polyg_sample) <- raster::crs(raster)

  # Crop the raster with the polygon
  rast_without_edge <- raster::crop(raster, polyg_sample)

  # ____________________________
  # ____________________________
  # ------ Code values selection

  # Copy the raster and remove values other than the class code
  # (allows for multiple values in class)
  r_class <- rast_without_edge
  r_class[which(!(raster::values(r_class) %in% class))] <- NA

  # Check whether there remains some non-NA values
  # otherwise return an error.
  if(length(unique(raster::values(r_class))) == 1){
    if(is.na(unique(raster::values(r_class)))){
      stop("The 'class' value must be a class code value
           from 'raster'")
    }
  }

  # ____________________________
  # ____________________________
  # ------ Clump to define habitat patches

  # Clump
  r_clump <- raster::clump(r_class, directions = neighborhood)

  # Get clump values corresponding to the ID of each patch
  val_cl <- as.vector(r_clump)
  # Remove NA values outside patches
  val_cl <- val_cl[-which(is.na(val_cl))]
  # Summarise the data to get the number of pixel per patch
  val_tab <- data.frame(table(val_cl))
  val_tab <- val_tab[-which(val_tab$Freq < surf_min), ]

  colnames(val_tab) <- c("ID", "area")

  val_tab$area <- val_tab$area * raster::res(raster)[1] * raster::res(raster)[2]

  return(val_tab)

}




