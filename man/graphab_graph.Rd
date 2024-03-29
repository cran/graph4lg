% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphab_graph.R
\name{graphab_graph}
\alias{graphab_graph}
\title{Create a graph in the Graphab project}
\usage{
graphab_graph(
  proj_name,
  linkset = NULL,
  name = NULL,
  thr = NULL,
  cost_conv = FALSE,
  proj_path = NULL,
  alloc_ram = NULL
)
}
\arguments{
\item{proj_name}{A character string indicating the Graphab project name.
The project name is also the name of the project directory in which the
file proj_name.xml is. It can be created with \code{\link{graphab_project}}}

\item{linkset}{(optional, default=NULL) A character string indicating the
name of the link set used to create the graph. If \code{linkset=NULL}, every
link set present in the project will be used to create a graph. Link sets
can be created with \code{\link{graphab_link}}.}

\item{name}{(optional, default=NULL) A character string indicating the
name of the graph created. If \code{name=NULL}, a name will be created. If
both \code{linkset=NULL} and \code{name=NULL}, then a graph will be created
for every link set present in the project and a name will be created every
time. In the latter case, a unique name cannot be specified. Link sets
can be created with \code{\link{graphab_link}}.}

\item{thr}{(optional, default=NULL) An integer or numeric value indicating
the maximum distance associated with the links of the created graph. It
allows users to create a pruned graph based on a distance threshold. Note that
when the link set used has a planar topology, the graph is necessarily a
pruned graph (not complete) and adding this threshold parameter can remove
other links. When the link set has been created with cost-distances, the
parameter is expressed in cost-distance units whereas when the link set is
based upon Euclidean distances, the parameter is expressed in meters.}

\item{cost_conv}{FALSE (default) or TRUE. Logical indicating whether numeric
\code{thr} values are converted from cost-distance into Euclidean distance
using a log-log linear regression. See also \code{\link{convert_cd}}
function.}

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
The function creates a graph from a link set in a Graphab project
}
\details{
By default, intra-patch distances are considered for metric
calculation. See more information in Graphab 2.8 manual:
\url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
}
\examples{
\dontrun{
graphab_graph(proj_name = "grphb_ex",
              linkset = "lcp",
              name = "graph")
}
}
\author{
P. Savary
}
