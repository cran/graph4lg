#' Plot scatterplots of genetic distance vs landscape distance
#'
#' @description The function enables to plot scatterplots to visualize the
#' relationship between genetic distance (or differentiation) and landscape
#' distance (Euclidean distance, cost-distance, etc.)between populations or
#' sample sites.
#'
#' @param mat_gd A symmetric \code{matrix} with pairwise genetic distances
#' between populations or sample sites.
#' @param mat_ld A symmetric \code{matrix} with pairwise landscape distances
#' between populations or sample sites. These distances can be
#' Euclidean distances, cost-distances or resistance distances, among others.
#' @param method A character string indicating the smoothing method
#' used to fit a line on the scatterplot. Possible values are the same as
#' with function 'geom_smooth()' from \pkg{ggplot2} : 'lm', 'glm', 'gam',
#' 'loess' (default).
#' @param thr_gd (optional) A numeric value used to remove values from the
#' data before to plot. All genetic distances values above \code{thr_gd}
#' are removed from the data.
#' @param thr_ld (optional) A numeric value used to remove values from the
#' data before to plot. All landscape distances values above \code{thr_ld}
#' are removed from the data.
#' @param se Logical (optional, default = TRUE) indicating whether the
#' confidence interval around the smooth line is displayed.
#' @param smooth_col (optional) A character string indicating the color
#' used to plot the smoothing line (default: "blue"). It must be a hexadecimal
#' color code or a color used by default in R.
#' @param pts_col (optional) Character string indicating the color
#' used to plot the points (default: "#999999"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @details IDs in \code{mat_gd} and \code{mat_ld} must be the same and refer
#' to the same sampling sites or populations, and both matrices must be ordered
#' in the same way.
#' Matrix of genetic distance \code{mat_gd} can be computed using
#' \code{\link{mat_gen_dist}}.
#' Matrix of landscape distance \code{mat_ld} can be computed using
#' \code{\link{mat_geo_dist}} when the landscape distance needed is a
#' Euclidean geographical distance.
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @author P. Savary
#' @examples
#' data(data_tuto)
#' mat_dps <- data_tuto[[1]]
#' mat_dist <- suppressWarnings(mat_geo_dist(data = pts_pop_simul,
#'       ID = "ID",
#'       x = "x",
#'       y = "y"))
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                       order(as.character(colnames(mat_dist)))]
#' scatterplot_ex <- scatter_dist(mat_gd = mat_dps,
#'                               mat_ld = mat_dist)

scatter_dist <- function(mat_gd,
                         mat_ld,
                         method = "loess",
                         thr_gd = NULL, thr_ld = NULL,
                         se = TRUE,
                         smooth_col = "black",
                         pts_col = "#999999"){

  # Check whether 'mat_gd' and 'mat_ld' are symmetric matrices
  if(!inherits(mat_gd, "matrix")){
    stop("'mat_gd' must be an object of class 'matrix'.")
  } else if (!inherits(mat_ld, "matrix")){
    stop("'mat_ld' must be an object of class 'matrix'.")
  } else if (!Matrix::isSymmetric(mat_gd)){
    stop("'mat_gd' must be a symmetric pairwise matrix.")
  } else if (!Matrix::isSymmetric(mat_ld)){
    stop("'mat_ld' must be a symmetric pairwise matrix.")
  }

  # Check whether 'mat_gd' and 'mat_ld' have the same rows' and columns' names
  if(!all(row.names(mat_gd) == row.names(mat_ld))){
    stop("'mat_gd' and 'mat_ld' must have the same rows' names.")
  } else if(!all(colnames(mat_gd) == colnames(mat_ld))){
    stop("'mat_gd' and 'mat_ld' must have the same columns' names.")
  }

  # Check whether 'method' is a character string
  if (!inherits(method, "character")){
    stop("'method' must be a character string.")
  }

  # Number of points
  nb_pts <- nrow(mat_gd)
  # Number of pairs (for a non-directed graph)
  nb_pairs <- nb_pts * (nb_pts - 1) / 2
  # Create vectors with the genetic and landscape distances
  # from the lower triangles of the distance matrices
  gen <- ecodist::lower(mat_gd)
  land <- ecodist::lower(mat_ld)

  # Remove values larger the thresholds if specified
  if (!is.null(thr_gd)){
    if(is.numeric(thr_gd)){
      if(thr_gd < max(gen)){
        gen[which(gen > thr_gd)] <- NA
      } else {
        message("There was not any genetic distance larger than 'thr_gd'.")
      }
    } else {
      stop("'thr_gd' must be a numeric value.")
    }
  }

  # Remove values larger the thresholds if specified
  if (!is.null(thr_ld)){
    if(is.numeric(thr_ld)){
      if(thr_ld < max(land)){
        land[which(land > thr_ld)] <- NA
      } else {
        message("There was not any landscape distance larger than 'thr_ld'.")
      }
    } else {
      stop("'thr_ld' must be a numeric value.")
    }
  }

  # Create a data.frame to store the values
  dat <- data.frame(gen = gen, land = land)

  # Remove NA values from the data.frame
  if(any(is.na(dat$gen))){
    dat <- dat[-which(is.na(dat$gen)), ]
  }

  if(any(is.na(dat$land))){
    dat <- dat[-which(is.na(dat$land)), ]
  }

  # Number of observations conserved
  nb_conserv <- nrow(dat)

  # Display the number of observations conserved
  message(paste(nb_conserv, " out of ", nb_pairs, " values were used.",
                sep = ""))

  # Create the plot
  scat <- ggplot(data = dat, aes(x = land, y = gen)) +
    geom_point(color = pts_col, size = 1, shape = 16) +
    geom_smooth(method = method, se = se, color = smooth_col) +
    labs(x = "Landscape distance",
         y = "Genetic distance") +
    theme_bw()

  return(scat)

}
