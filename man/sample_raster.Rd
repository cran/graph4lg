% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_raster.R
\name{sample_raster}
\alias{sample_raster}
\title{Sample points or patches on a categorical raster layer}
\usage{
sample_raster(
  raster,
  class,
  nb_pts,
  dist_min = 0,
  edge_size = 0,
  by_patch = TRUE,
  neighborhood = 8,
  surf_min = 0,
  prop_area = TRUE,
  step_max = 1000,
  output = "df",
  desc = TRUE
)
}
\arguments{
\item{raster}{A RasterLayer object corresponding to a categorical raster layer}

\item{class}{An integer value or vector with the value(s) corresponding to
the code values of the raster layer within which points will be sampled.}

\item{nb_pts}{An integer value indicating the number of points to be sampled}

\item{dist_min}{An integer value indicating the minimum distance separating
the sampled points (default = 0).}

\item{edge_size}{An integer value indicating the width of the edge of the
raster layer which is ignored during the sampling (default = 0). It prevents
from sampling in the margins of the study area.}

\item{by_patch}{A logical value indicating whether contiguous patches with
cells having the same code value are delineated prior to
sampling (default = TRUE). It prevents from sampling several points in the
same contiguous patch.}

\item{neighborhood}{An integer value indicating which cells are considered
adjacent when contiguous patches are delineated (it should be 8
(default, Queen's case) or 4 (Rook's case)). This parameter is ignored
when \code{by_patch = FALSE}.}

\item{surf_min}{An integer value indicating the minimum surface of a patch
considered for the sampling in number of raster cells. This parameter is used
whatever the \code{by_patch} argument is. Default is 0.}

\item{prop_area}{A logical value indicating whether sampling in large patches
is more likely (default = TRUE). If \code{by_patch = FALSE}, this parameter
is ignored. When \code{prop_area = TRUE}, the probability to sample a given
patch is proportional to its area.}

\item{step_max}{An integer value indicating how many sampling steps are
performed to identify a point set satisfying all the conditions before
returning an error.}

\item{output}{A character string indicating the type of returned output:
\itemize{
\item{'data.frame': A \code{data.frame} with three/four columns:
\itemize{
\item{ID: The point or patch centroid ID}
\item{x: The point or patch centroid longitude}
\item{y: The point or patch centroid latitude}
\item{area: The area of the sampled patch (only if \code{by_patch = TRUE})}
}}
\item{'pts_layer': A \code{SpatialPointsDataFrame} layer corresponding
to the sampled point (points or patch centroids)}
\item{'poly_layer': A \code{SpatialPolygonsDataFrame} layer corresponding
to the sampled patch polygons}
}}

\item{desc}{A logical value indicating whether the result should be
described or not (default = FALSE). If \code{desc = TRUE}, then the Gini
coefficient of the distances between points and of the patch areas (if
\code{by_patch = TRUE}) is computed with the \code{\link{gini_coeff}}.
An histogram of the link weights is also described.}
}
\value{
A list of object(s) with one or several elements according to the
\code{output} and \code{desc} arguments.
}
\description{
The function samples points or patches on a categorical raster
layer.
}
\author{
P. Savary
}
\keyword{internal}
