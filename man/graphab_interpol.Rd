% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphab_interpol.R
\name{graphab_interpol}
\alias{graphab_interpol}
\title{Creates a raster with interpolated connectivity metric values from metrics
already computed in the Graphab project}
\usage{
graphab_interpol(
  proj_name,
  name,
  reso,
  linkset,
  graph,
  var,
  dist,
  prob = 0.05,
  thr = NULL,
  summed = FALSE,
  proj_path = NULL,
  alloc_ram = NULL
)
}
\arguments{
\item{proj_name}{A character string indicating the Graphab project name.
The project name is also the name of the project directory in which the
file proj_name.xml is. It can be created with \code{\link{graphab_project}}}

\item{name}{A character string indicating the name of the raster to be
created after the interpolation.}

\item{reso}{An integer indicating the spatial resolution in meters of the
raster resulting from the metric interpolation.}

\item{linkset}{A character string indicating the name of the link set used
for the interpolation. It should be the one used to create the used graph
and the metric.}

\item{graph}{A character string indicating the name of the graph from which
the metric was computed and whose links are considered for a potential
multi-linkage with patches.
This graph has been created with Graphab or using \code{\link{graphab_graph}}
function and is associated with a link set.}

\item{var}{A character string indicating the name of the already computed
metric to be interpolated.}

\item{dist}{A numeric or integer value specifying the distance at which we
assume a probability equal to \code{prob} during the interpolation.
It is used to set \eqn{\alpha} for computing probabilities associated
with distances between each pixel and the neighboring patch(es) such that
probability between patch i and pixel j is \eqn{p_{ij}= e^{-\alpha d_{ij}}}.}

\item{prob}{A numeric or integer value specifying the probability
at distance \code{dist}. By default, \code{code=0.05}. It is used to set
\eqn{\alpha} (see param \code{dist} above).}

\item{thr}{(default NULL) If NULL, the value of each pixel is computed from
the value of the metric at the nearest habitat patch, weighted by a
probability depending on distance. If an integer, the value of each pixel
depends on the values of the metric taken at several of the nearest habitat
patches, up to a distance (cost or Euclidean distance, depending on the type
of linkset) equal to \code{thr}.}

\item{summed}{Logical (default = FALSE) only used if \code{thr} is not NULL,
and specifying whether multiple values are summed up (TRUE) or averaged
after being weighted by probabilities.}

\item{proj_path}{(optional) A character string indicating the path to the
directory that contains the project directory. It should be used when the
project directory is not in the current working directory. Default is NULL.
When 'proj_path = NULL', the project directory is equal to \code{getwd()}.}

\item{alloc_ram}{(optional, default = NULL) Integer or numeric value
indicating RAM gigabytes allocated to the java process. Increasing this
value can speed up the computations. Too large values may not be compatible
with your machine settings.}
}
\description{
The function creates a raster with interpolated connectivity
metric values from a metric already computed in the Graphab project.
}
\details{
See more information in Graphab 2.8 manual:
\url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
Be careful, when capacity has been changed. The last changes are taken into
account for subsequent calculations in a project.
}
\examples{
\dontrun{
graphab_interpol(proj_name = "grphb_ex",
                 name = "F_interp",
                 reso = 20,
                 linkset = "lcp",
                 graph = "graph",
                 var = "F_d600_p0.5_beta1_graph",
                 dist = 600,
                 prob = 0.5)
}
}
\author{
P. Savary
}
