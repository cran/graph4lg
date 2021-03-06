% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphab_to_igraph.R
\name{graphab_to_igraph}
\alias{graphab_to_igraph}
\title{Create landscape graphs from Graphab link set}
\usage{
graphab_to_igraph(
  proj_name,
  linkset,
  nodes = "patches",
  weight = "cost",
  proj_path = NULL,
  fig = FALSE,
  crds = FALSE
)
}
\arguments{
\item{proj_name}{A character string indicating the project name. It is also
the name of the directory in which proj_name.xml file is found. By default,
'proj_name' is searched into the current working directory}

\item{linkset}{A character string indicating the name of the linkset used to
create the graph links. The linkset must have been created previously (see
the function \code{\link{graphab_link}}). It can be complete or planar. The
graph is given the topology of the selected link set.}

\item{nodes}{A character string indicating whether the nodes of the created
graph are given all the attributes or metrics computed in Graphab or only
those specific to a given graph previously created with
\code{\link{graphab_graph}}
It can be:\itemize{
\item{\code{nodes = "patches"}(default): all the attributes and metrics of
the habitat patches are included as node attributes in \code{igraph} object.}
\item{\code{nodes = "graph_name"}(default): only the metrics of
the habitat patches computed from the graph 'graph_name' created with
\code{\link{graphab_graph}} are included as node attributes in
\code{igraph} object, along with some basic patch attributes.}
}}

\item{weight}{A character string ("euclid" or "cost") indicating
whether to weight the links with Euclidean distance or
cost-distance (default) values.}

\item{proj_path}{(optional) A character string indicating the path to the
directory that contains the project directory ('proj_name'). By default,
'proj_name' is searched into the current working directory}

\item{fig}{Logical (default = FALSE) indicating whether to plot a figure of
the resulting spatial graph. The figure is plotted using function
\code{\link{plot_graph_lg}}. The plotting can be long if the graph has many
nodes and links.}

\item{crds}{Logical (default = FALSE) indicating whether to create an object
of class \code{data.frame} with the node centroid spatial coordinates. Such a
\code{data.frame} has 3 columns: 'ID', 'x', 'y'.}
}
\value{
A graph object of class \code{igraph} (if crds = FALSE) or a
list of objects: a graph object of class \code{igraph} and a
\code{data.frame} with the nodes spatial coordinates (if crds = TRUE).
}
\description{
The function creates a landscape graph from a link set created
with Graphab software or different functions of this package and converts
it into a graph object of class \code{igraph}.
The graph has weighted links and is undirected.
Nodes attributes present in the Graphab project are included, including
connectivity metrics when computed
}
\examples{
\dontrun{
proj_path <- system.file('extdata',package='graph4lg')
proj_name <- "grphb_ex"
linkset <- "lkst1"
nodes <- "graph"
graph <- graphab_to_igraph(proj_name = proj_name,
                           linkset = "lkst1",
                           nodes = "graph",
                           links = links,
                           weights = "cost",
                           proj_path = proj_path,
                           crds = FALSE,
                           fig = FALSE)
                           }
}
\references{
\insertRef{foltete2012software}{graph4lg}
}
\author{
P. Savary
}
