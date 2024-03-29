% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mat_cost_dist.R
\name{mat_cost_dist}
\alias{mat_cost_dist}
\title{Compute cost distances between points on a raster}
\usage{
mat_cost_dist(
  raster,
  pts,
  cost,
  method = "gdistance",
  return = "mat",
  direction = 8,
  parallel.java = 1,
  alloc_ram = NULL
)
}
\arguments{
\item{raster}{A parameter indicating the raster file on which cost distances
are computed. It can be:\itemize{
\item{A character string indicating the path to a raster file in format
.tif or .asc.}
\item{A \code{RasterLayer} object already loaded in R environment}
}
All the raster cell values must be present in the column 'code' from
\code{cost} argument.}

\item{pts}{A parameter indicating the points between which cost distances
are computed. It can be either: \itemize{
\item{A character string indicating the path to a .csv file. It must have
three columns:\itemize{
\item{ID: The ID of the points.}
\item{x: A numeric or integer indicating the longitude of the points.}
\item{y: A numeric or integer indicating the latitude of the points.}
}}
\item{A \code{data.frame} with the spatial coordinates of the points.
It must have three columns:\itemize{
\item{ID: The ID of the points.}
\item{x: A numeric or integer indicating the longitude of the points.}
\item{y: A numeric or integer indicating the latitude of the points.}
}}
\item{A \code{SpatialPointsDataFrame} with at least an attribute column
named "ID" with the point IDs.}
}
The point coordinates must be in the same spatial coordinate reference system
as the raster file.}

\item{cost}{A \code{data.frame} indicating the cost values associated to each
raster value. It must have two columns:\itemize{
\item{'code': raster cell values}
\item{'cost': corresponding cost values}
}}

\item{method}{A character string indicating the method used to compute the
cost distances. It must be:\itemize{
\item{'gdistance': uses the functions from the package \pkg{gdistance}
assuming that movement is possible in 8 directions from each cell, that
a geo-correction is applied to correct for diagonal movement lengths and that
raster cell values correspond to resistance (and not conductance).}
\item{'java': uses a .jar file which is downloaded on the user's machine if
necessary and if java is installed. This option substantially reduces
computation times and makes possible the parallelisation.}
}}

\item{return}{A character string indicating whether the returned object is a
\code{data.frame} (\code{return="df"}) or a pairwise
\code{matrix} (\code{return="mat"}).}

\item{direction}{An integer (4, 8, 16) indicating the directions in which
movement can take place from a cell. Only used when \code{method="gdistance"}.
By default, \code{direction=8}.}

\item{parallel.java}{An integer indicating how many computer cores are used
to run the .jar file. By default, \code{parallel.java=1}.}

\item{alloc_ram}{(optional, default = NULL) Integer or numeric value
indicating RAM gigabytes allocated to the java process when used. Increasing
this value can speed up the computations. Too large values may not be
compatible with your machine settings.}
}
\value{
The function returns:\itemize{
\item{If \code{return="mat"}, a pairwise \code{matrix} with cost-distance
values between points.}
\item{If \code{return="df"}, an object of type \code{data.frame} with three columns:
\itemize{
\item{from: A character string indicating the ID of the point of origin.}
\item{to: A character string indicating the ID of the point of destination.}
\item{cost_dist: A numeric indicating the accumulated cost-distance along
the least-cost path between point ID1 and point ID2.}
}
}
}
}
\description{
The function computes cost-distances associated to least cost
paths between point pairs on a raster with specified cost values.
}
\examples{
\dontrun{
x <- raster::raster(ncol=10, nrow=10, xmn=0, xmx=100, ymn=0, ymx=100)
raster::values(x) <- sample(c(1,2,3,4), size = 100, replace = TRUE)
pts <- data.frame(ID = 1:4,
                  x = c(10, 90, 10, 90),
                  y = c(90, 10, 10, 90))
cost <- data.frame(code = 1:4,
                   cost = c(1, 10, 100, 1000))
mat_cost_dist(raster = x,
              pts = pts, cost = cost,
              method = "gdistance")
}
}
\author{
P. Savary
}
