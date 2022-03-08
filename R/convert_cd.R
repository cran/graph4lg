#' Fit a model to convert cost-distances into Euclidean distances
#'
#' @description The function fits a model to convert cost-distances into
#' Euclidean distances as implemented in Graphab software.
#'
#' @param mat_euc A symmetric \code{matrix} or \code{dist} object with
#' pairwise geographical Euclidean distances between populations or sample
#' sites. It will be the explanatory variable, and only values from the off
#' diagonal lower triangle will be used.
#' @param mat_ld A symmetric \code{matrix} or \code{dist} object with pairwise
#' landscape distances between populations or sample sites. These distances can
#' be cost-distances or resistance distances, among others. It will be the
#' explained variable, and only values from the off diagonal lower triangle
#' will be used.
#' @param method A character string indicating the method used to fit the model.
#' \itemize{
#' \item{If 'method = "log-log"' (default), then the model takes the
#' following form : log(ld) ~ A + B * log(euc)}
#' \item{If 'method = "lm"', then the model takes the following form :
#' ld ~ A + B * euc}
#' }
#' @param to_convert A numeric value or numeric vector with Euclidean distances
#' to convert into cost-distances.
#' @param fig Logical (default = TRUE) indicating whether a figure is plotted
#' @param line_col (if 'fig = TRUE') Character string indicating the color
#' used to plot the line (default: "blue"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @param pts_col (if 'fig = TRUE') Character string indicating the color
#' used to plot the points (default: "#999999"). It must be a hexadecimal color
#' code or a color used by default in R.
#' @details IDs in 'mat_euc' and 'mat_ld' must be the same and refer to the same
#' sampling site or populations, and both matrices must be ordered
#' in the same way.
#' Matrix of Euclidean distance 'mat_euc' can be computed using the function
#' \code{\link{mat_geo_dist}}.
#' Matrix of landscape distance 'mat_ld' can be computed using the function
#' \code{\link{mat_cost_dist}}.
#' Before the log calculation, 0 distance values are converted into 1,
#' so that they are 0 after this calculation.
#' @return A list of output (converted values, estimated parameters, R2)
#' and optionally a ggplot2 object to plot
#' @import ggplot2
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' @examples
#' data("data_tuto")
#' mat_ld <- data_tuto[[2]][1:10, 1:10] * 1000
#' mat_euc <- data_tuto[[1]][1:10, 1:10] * 50000
#' to_convert <- c(30000, 40000)
#' res <- convert_cd(mat_euc = mat_euc,
#'                   mat_ld = mat_ld,
#'                   to_convert = to_convert, fig = FALSE)


convert_cd <- function(mat_euc, mat_ld,
                       to_convert,
                       method = "log-log",
                       fig = TRUE,
                       line_col = "black",
                       pts_col = "#999999"){

  # Check whether mat_euc and mat_ld are symmetric matrices or dist objects
  if(!inherits(mat_euc, c("matrix", "dist"))){
    stop("'mat_euc' must be an object of class 'matrix' or 'dist'.")
  } else if (!inherits(mat_ld, c("matrix", "dist"))){
    stop("'mat_lc' must be an object of class 'matrix' or 'dist'.")
  } else if (inherits(mat_euc, "matrix")){
    if(!Matrix::isSymmetric(mat_euc)){
      stop("'mat_euc' must be a symmetric pairwise matrix.")
    }
  } else if (inherits(mat_ld, "matrix")){
    if (!Matrix::isSymmetric(mat_ld)){
      stop("'mat_ld' must be a symmetric pairwise matrix.")
    }
  } else if (inherits(mat_euc, "dist")){
    mat_euc <- as.matrix(mat_euc)
  } else if (inherits(mat_ld, "dist")){
    mat_ld <- as.matrix(mat_ld)
  }

  # Check whether mat_euc and mat_lad have same row names and column names
  # and are ordered in the same way
  if(!all(row.names(mat_euc) == row.names(mat_ld))){
    stop("'mat_euc' and 'mat_dist' must have the same row names.")
  } else if(!all(colnames(mat_euc) == colnames(mat_ld))){
    stop("'mat_euc' and 'mat_ld' must have the same column names.")
  }

  df <- data.frame(euc = ecodist::lower(mat_euc),
                   ld = ecodist::lower(mat_ld))

  if(method == "log-log"){

    if(nrow(df[which(df$euc == 0), ])){
      message("Euclidean distances equal to zero were replaced by 1 before
              computing the logarithm of the distances.")
      df[which(df$euc == 0), 'euc'] <- 1
    }

    if (nrow(df[which(df$ld == 0), ])){
      message("Landscape distances equal to zero were replaced by 1 before
              computing the logarithm of the distances.")

      df[which(df$ld == 0), 'ld'] <- 1
    }


    df$log_ld <- log(df$ld)
    df$log_euc <- log(df$euc)
    model <- stats::lm(data = df, formula = log_ld ~ log_euc)

    converted <- exp(stats::predict(model,
                                    newdata = data.frame(log_euc = log(to_convert))))


  } else if (method == "lm"){
    model <- stats::lm(data = df, formula = ld ~ euc)

    converted <- stats::predict(model,
                                newdata = data.frame(euc = to_convert))

  } else {
    stop("You must specify a correct 'method' option ('log-log' or 'lm').")
  }

  names(converted) <- to_convert

  param <- model$coefficients
  names(param) <- c("Intercept", "Slope")

  r_sq <- as.numeric(summary(model)[8])
  names(r_sq) <- "Multiple R-squared"

  if(fig == FALSE){
    list_res <- list(converted, param, r_sq)
    names(list_res) <- c("Converted values",
                         "Model parameters",
                         "Model multiple R-squared")
  } else if(fig == TRUE){

    df$predict <- stats::predict(model)

    # In the plot, the points are the initial values
    # The line is drawn from the values predicted by the model

    if(method == "log-log"){

      df2 <- data.frame(log_euc = log(to_convert),
                        log_conv = log(converted))

      plot <- ggplot() +
        geom_point(data = df,
                   aes(x = .data$log_euc, y = .data$log_ld), color = pts_col)+
        geom_line(data = df, aes(x = .data$log_euc, y = .data$predict),
                  color = line_col)+
        geom_point(data = df2,
                   aes(x = .data$log_euc, y = .data$log_conv),
                   color = "black", size = 3)+
        labs(x = "log ( Euclidean distance )",
             y = "log ( Landscape distance )") +
        theme_bw()

    } else {

      df2 <- data.frame(euc = to_convert,
                        conv = converted)

      plot <- ggplot() +
        geom_point(data = df,
                   aes(x = .data$euc, y = .data$ld),
                   color = pts_col)+
        geom_line(data = df,
                  aes(x = .data$euc, y = .data$predict),
                  color = line_col)+
        geom_point(data = df2,
                   aes(x = .data$euc, y = .data$conv),
                   color = "black", size = 3)+
        labs(x = "Euclidean distance",
             y = "Landscape distance") +
        theme_bw()

    }

    list_res <- list(converted, param, r_sq, plot)
    names(list_res) <- c("Converted values", "Model parameters",
                         "Model multiple R-squared", "Plot")

  } else {
    stop("'fig' must be logical (TRUE or FALSE).")
  }


  return(list_res)

}



