#' Convert degrees to radians
#'
#' @description The function converts degree to radians
#'
#' @param deg A coordinate in degrees
#' @return The coordinate in radians
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' deg2rad(40.75170)

deg2rad <- function(deg){
  rad <- deg*pi/180
  return(rad)
}


#' Calculate the Great-Circle distance between two points using the
#' Spherical Law of Cosines (slc)
#'
#' @description The function calculates the Great-Circle distance between two
#' points specified by radian latitude/longitude using the Spherical Law
#' of Cosines (slc)
#'
#' @param long1 Point 1 longitude in radians
#' @param lat1 Point 1 latitude in radians
#' @param long2 Point 2 longitude in radians
#' @param lat2 Point 2 latitude in radians
#' @return The distance between points 1 and 2 in meters
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' dist_gc_slc(long1 = -73.99420, lat1 = 40.75170,
#'             long2 = -87.63940, lat2 = 41.87440)

dist_gc_slc <- function(long1, lat1, long2, lat2) {

  R <- 6371
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R * 1000
  return(d)
}



#' Calculate the Great-Circle distance between two points using the
#' Harversine formula (hvs)
#'
#' @description The function calculates the Great-Circle distance between two
#' points specified by radian latitude/longitude using the
#' Harversine formula (hvs)
#'
#' @param long1 Point 1 longitude in radians
#' @param lat1 Point 1 latitude in radians
#' @param long2 Point 2 longitude in radians
#' @param lat2 Point 2 latitude in radians
#' @return The distance between points 1 and 2 in meters
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' dist_gc_hvs(long1 = -73.99420, lat1 = 40.75170,
#'             long2 = -87.63940, lat2 = 41.87440)


dist_gc_hvs <- function(long1, lat1, long2, lat2) {

  R <- 6371
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c * 1000
  return(d)
}



#' Calculate the Great-Circle distance between two points using the
#' Vincenty inverse formula for ellipsoids (vicenty)
#'
#' @description The function calculates the Great-Circle distance between two
#' points specified by radian latitude/longitude using the
#' Vincenty inverse formula for ellipsoids (vicenty)
#'
#' @param long1 Point 1 longitude in radians
#' @param lat1 Point 1 latitude in radians
#' @param long2 Point 2 longitude in radians
#' @param lat2 Point 2 latitude in radians
#' @return The distance between points 1 and 2 in meters
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' dist_gc_vicenty(long1 = -73.99420, lat1 = 40.75170,
#'             long2 = -87.63940, lat2 = 41.87440)


dist_gc_vicenty <- function(long1, lat1, long2, lat2) {

  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # ength of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid

  L <- long2-long1 # difference in longitude
  U1 <- atan((1-f) * tan(lat1)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)

  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL

  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
      (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                             B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  d <- b*A*(sigma-deltaSigma)

  return(d)
}



#' Compute the Great Circle distance between two points
#'
#' @description The function computes the Great Circle distance between two
#' two points defined by their longitudes and latitudes.
#'
#' @param long1 project name, project dir in which proj_name.xml is found
#' @param long2 raster.tif INT2S path or present in wd,
#' @param lat1 habitat code in the raster file
#' @param lat2 default 0, minimum habitat size in ha
#' @param method default NULL nodata code in the raster file
#' @keywords internal
#' @export
#' @author P. Savary
#' @examples
#' dist_great_circle(long1 = -73.99420,
#'                   lat1 = 40.75170,
#'                   long2 = -87.63940,
#'                   lat2 = 41.87440,
#'                   method = "vicenty")



dist_great_circle <- function(long1, long2, lat1, lat2,
                              method = "vicenty"){

  long1 <- deg2rad(long1)
  long2 <- deg2rad(long2)
  lat1 <- deg2rad(lat1)
  lat2 <- deg2rad(lat2)


  if(long1 == long2 & lat1 == lat2){
    d <- 0
  } else {

    if(method == "slc"){

      d <- dist_gc_slc(long1, lat1, long2, lat2)

    } else if(method == "hvs"){

      d <- dist_gc_hvs(long1, lat1, long2, lat2)

    } else if (method == "vicenty"){

      d <- dist_gc_vicenty(long1, lat1, long2, lat2)

    } else {
      stop("'method' must be 'slc', 'hvs' or 'vicenty'.")
    }
  }

  return(d)

}





