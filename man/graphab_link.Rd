% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphab_link.R
\name{graphab_link}
\alias{graphab_link}
\title{Create a link set in the Graphab project}
\usage{
graphab_link(
  proj_name,
  distance = "cost",
  name,
  cost = NULL,
  topo = "planar",
  remcrosspath = FALSE,
  proj_path = NULL,
  alloc_ram = NULL
)
}
\arguments{
\item{proj_name}{A character string indicating the Graphab project name.
The project name is also the name of the project directory in which the
file proj_name.xml is. It can be created with \code{\link{graphab_project}}}

\item{distance}{A character string indicating whether links between patches
are computed based on:\itemize{
\item{Shortest cost distances: \code{distance='cost'} (default)}
\item{Straight Euclidean distances: \code{distance='euclid'}}
}
In the resulting link set, each link will be associated with its
corresponding cost-distance and the length of the least-cost path in meters
(if \code{distance='cost'}) or with its length in Euclidean distance
(if \code{distance='euclid'})}

\item{name}{A character string indicating the name of the created linkset.}

\item{cost}{This argument could be:\itemize{
\item{A \code{data.frame} indicating the cost values associated to each
raster cell value. These values refer to the raster used to create the
project with \code{graphab_project}. The data.frame must have two
columns:\itemize{
\item{'code': raster cell values}
\item{'cost': corresponding cost values}
}}
\item{The path to an external raster file in .tif format with cost values.}
}}

\item{topo}{A character string indicating the topology of the created
link set. It can be:\itemize{
\item{Planar (\code{topo='planar'} (default)): a planar set of links is
created. It speeds up the computation but will prevent from creating
complete graphs with \code{\link{graphab_graph}}.}
\item{Complete (\code{topo='complete'}): a complete set of links is created.
A link is computed between every pair of patches.}
}}

\item{remcrosspath}{(optional, default = FALSE) A logical indicating whether
links crossing patches are removed (TRUE).}

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
The function creates a link set between habitat patches in the
Graphab project.
}
\details{
By default, links crossing patches are not ignored nor broken into
two links. For example, a link from patches A to C crossing patch B
is created. It takes into account the distance inside patch B. It can be a
problem when computing BC index. See more information in Graphab 2.8 manual:
\url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
}
\examples{
\dontrun{
df_cost <- data.frame(code = 1:5,
                      cost = c(1, 10, 100, 1000, 1))
graphab_link(proj_name = "grphb_ex",
            distance = "cost",
            name = "lcp",
            cost = df_cost,
            topo = "complete")
}
}
\author{
P. Savary, T. Rudolph
}
