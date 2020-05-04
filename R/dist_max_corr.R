#' Compute the distance at which the correlation between genetic distance
#' and landscape distance is maximal
#'
#' @description The function enables to compute the distance at which the
#' correlation between genetic distance and landscape distance is maximal,
#' using a method similar to that employed by van Strien et al. (2015).
#' Iteratively, distance threshold values are tested. For each value, all the
#' population pairs separated by a landscape distance larger than the threshold
#' are removed before the Mantel correlation coefficient between genetic
#' distance and landscape distance is computed.
#' The distance threshold at which the correlation is the strongest is then
#' identified. A figure showing the evolution of the correlation coefficients
#' when landscape distance threshold increases is plotted.
#'
#' @param mat_gd A symmetric \code{matrix} with pairwise genetic distances
#' between populations or sample sites.
#' @param mat_ld A symmetric \code{matrix} with pairwise landscape distances
#' between populations or sample sites. These distances can be
#' Euclidean distances, cost-distances or resistance distances, among others.
#' @param interv A numeric value indicating the interval between the different
#' distance thresholds for which the correlation coefficients are computed.
#' @param from (optional) The minimum distance threshold value at which the
#' correlation coefficient is computed.
#' @param to (optional) The maximum distance threshold value at which the
#' correlation coefficient is computed.
#' @param thr_gd (optional) A numeric value used to remove
#' genetic distance values from the data before the calculation.
#' All genetic distances values above 'thr_gd' are removed from the data.
#' This parameter can be used especially when there are outliers.
#' @param fig Logical (default = TRUE) indicating whether a figure is plotted.
#' @param pts_col (optional, if fig = TRUE) A character string indicating the
#' color used to plot the points (default: "#999999"). It must be a hexadecimal
#' color code or a color used by default in R.
#' @param line_col (optional, if fig = TRUE) A character string indicating the
#' color used to plot the line (default: "blue"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @details IDs in 'mat_gd' and 'mat_ld' must be the same and refer to the same
#' sampling sites or populations, and both matrices must be ordered
#' in the same way.
#' The correlation coefficient between genetic distance and landscape distance
#' computed is a Mantel correlation coefficient. If there are less than 50
#' pairwise values, the correlation is not computed, as in
#' van Strien et al. (2015). Such a method can be subject to criticism from
#' a strict statistical point of view given correlation coefficients computed
#' from samples of different size are compared.
#' The matrix of genetic distance 'mat_gd' can be computed using
#' \code{\link{mat_gen_dist}}.
#' The matrix of landscape distance 'mat_ld' can be computed using
#' \code{\link{mat_geo_dist}} when the landscape distance needed is a
#' Euclidean geographical distance.
#' Mantel correlation coefficients are computed using
#' the function \code{\link[vegan]{mantel}}.
#' @return A list of objects:
#' \itemize{
#' \item{The distance at which the correlation is the highest.}
#' \item{The vector of correlation coefficients at the different
#' distance thresholds}
#' \item{The vector of the different distance thresholds}
#' \item{A ggplot2 object to plot}
#' }
#' @import ggplot2
#' @export
#' @author P. Savary
#' @references \insertRef{van2015isolation}{graph4lg}
#' @examples
#' data("data_tuto")
#' mat_gen <- data_tuto[[1]]
#' mat_dist <- data_tuto[[2]]*1000
#' res_dmc <- dist_max_corr(mat_gd = mat_gen,
#'                          mat_ld = mat_dist,
#'                          from = 32000, to = 42000,
#'                          interv = 5000,
#'                          fig = FALSE)


dist_max_corr <- function(mat_gd, mat_ld,
                          interv,
                          from = NULL, to = NULL,
                          fig = TRUE,
                          thr_gd = NULL,
                          line_col = "black",
                          pts_col = "#999999"){

  # Check whether mat_gd and mat_ld are symmetric matrices
  if(!inherits(mat_gd, "matrix")){
    stop("'mat_gd' must be an object of class 'matrix'.")
  } else if (!inherits(mat_ld, "matrix")){
    stop("'mat_ld' must be an object of class 'matrix'.")
  } else if (!Matrix::isSymmetric(mat_gd)){
    stop("'mat_gd' must be a symmetric pairwise matrix.")
  } else if (!Matrix::isSymmetric(mat_ld)){
    stop("'mat_ld' must be a symmetric pairwise matrix.")
  }

  # Check whether mat_gd and mat_lad have same row names and column names
  # and are ordered in the same way
  if(!all(row.names(mat_gd) == row.names(mat_ld))){
    stop("'mat_gd' and 'mat_dist' must have the same row names.")
  } else if(!all(colnames(mat_gd) == colnames(mat_ld))){
    stop("'mat_gd' and 'mat_ld' must have the same column names.")
  }

  # Check whether interv is numeric
  if (!inherits(interv, "numeric")){
    stop("'interv' must be a numeric value.")
  }

  # Get maximum and minimum landscape distances values
  max_ld <- max(mat_ld, na.rm = TRUE)
  min_ld <- min(mat_ld, na.rm = TRUE)
  # Calculate the number of intervals
  n_interv_max <- (max_ld - min_ld) / interv

  # Check whether there are at least 5 distances thresholds
  if(n_interv_max < 5){
    warning("Less than 5 correlation coefficients can be computed given
            this 'interv' value.")
  }

  # If from and to are not specified, then they are
  # the minimum and maximum distance values in mat_ld, respectively
  if(all(c(is.null(from), is.null(to)))){
    from <- min_ld
    to <- max_ld
  } else if(!all(c(is.numeric(from), is.numeric(to)))){
    stop("'from' and 'to' must be numeric values.")
  }

  # Create a vector with the threshold values
  if(interv > from){
    # If interv > from, then the first threshold value is interv
    t1 <- interv
    # t2 is equal to the last threshold value
    t2 <- interv * floor(to / interv)
    vec_t <- seq(from = t1, to = t2,
                 by = interv)
  } else {
    # If interv < from, then the first threshold value is the first
    # multiple of interv larger than from
    t1 <- interv * ceiling(from / interv)
    t2 <- interv * floor(to / interv)
    vec_t <- seq(from = t1, to = t2,
                 by = interv)
  }

  # The last element of vec_t is the maximum landscape distance
  vec_t[length(vec_t) + 1] <- max_ld
  # There are nb_t different intervals
  nb_t <- length(vec_t)

  # Vector of correlation coefficients
  cc_val <- c()

  if (!is.null(thr_gd)){
    if (is.numeric(thr_gd)){
      if(thr_gd < max(mat_gd, na.rm = TRUE)){
        mat_gd2 <- mat_gd
        mat_gd2[mat_gd2 > thr_gd] <- NA
      } else {
        message("There was not any genetic distance larger than 'thr_gd'.")
      }
    } else {
      stop("'thr_gd' must be either NULL or a numeric value.")
    }
  } else {
    mat_gd2 <- mat_gd
  }


  for (i in (1:nb_t)){
    mat_l1 <- mat_ld
    mat_l1[mat_l1 > vec_t[i]] <- NA

    val <- ecodist::lower(mat_l1)

    if(length(val[!is.na(val)]) > 50){
      mant_res <- vegan::mantel(mat_l1, mat_gd2, na.rm = TRUE)
      cc_val[i] <- mant_res[[3]]
    } else {
      cc_val[i] <- NA
    }

  }

  id_max <- which(cc_val == max(cc_val, na.rm = TRUE))
  dist_max <- vec_t[id_max]
  dist_max <- dist_max[1]


  if(fig == TRUE){

    if (!inherits(pts_col, "character")){
      stop("'pts_col' must be a character string.")
    }


    if (!inherits(line_col, "character")){
      stop("'line_col' must be a character string.")
    }

    dat <- data.frame(vec_t = vec_t, cc_val = cc_val)
    if(any(is.na(dat$cc_val))){
      dat <- dat[-which(is.na(dat$cc_val)), ]
    }

    plot_dmc <- ggplot(data = dat, aes(x = vec_t, y = cc_val)) +
      geom_point(color = pts_col, size = 1, shape = 16) +
      geom_line(color = line_col) +
      labs(x = "Distance threshold",
           y = "Correlation coefficient") +
      theme_bw()
    print(plot_dmc)
  }

  #message(paste("Distance threshold at which correlation reaches
  # a maximum: ", dist_max, sep = ""))
  #message(paste("Maximum correlation coefficient between mat_gd and mat_ld:
  # ", cc_val[id_max], sep = ""))


  if (fig == TRUE){
    res_list <- list(dist_max, cc_val, vec_t, plot_dmc)
    names(res_list) <- c("distance at which correlation reaches a maximum",
                         "correlation coefficients",
                         "distance thresholds", "plot_dmc")
  } else {
    res_list <- list(dist_max, cc_val, vec_t)
    names(res_list) <- c("distance at which correlation reaches a maximum",
                         "correlation coefficients",
                         "distance thresholds")
  }

  return(res_list)

}
