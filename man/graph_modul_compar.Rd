% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_modul_compar.R
\name{graph_modul_compar}
\alias{graph_modul_compar}
\title{Compare the partition into modules of two graphs}
\usage{
graph_modul_compar(
  x,
  y,
  mode = "graph",
  nb_modul = NULL,
  algo = "fast_greedy",
  node_inter = "distance",
  data = NULL
)
}
\arguments{
\item{x}{The first graph object
\itemize{
\item{If \code{mode = 'graph'} (default), \code{x} is a graph object of
class \code{igraph}.
Then, its nodes must have the same names as in graph \code{y}.}
\item{If \code{mode = 'data.frame'}, \code{x} refers to a column of
the \code{data.frame} 'data'.
Then \code{x} must be a character string indicating the name of the
column of 'data' with the modules' labels of the nodes in the first graph.
In that case, the column can be of class \code{numeric}, \code{character}
or \code{factor} but will be converted into a \code{numeric} vector
in any case.}
\item{If \code{mode = 'vector'}, \code{x} is a vector of
class \code{character}, \code{factor} or \code{numeric}.
In that case, it must have the same length as vector \code{y} and
will be converted into a \code{numeric} vector.}
}}

\item{y}{The second graph object
Same classes possible as for \code{x}. Must be of the same format as \code{x}}

\item{mode}{A character string indicating whether x and y are igraph objects,
vectors or columns from a data.frame. \code{mode} can be 'graph',
'data.frame' or 'vector'.}

\item{nb_modul}{(if x and y are igraph objects) A numeric or integer value
or a numeric vector with 2 elements indicating the number of modules to
create in both graphs.
\itemize{
\item{If \code{nb_modul} is a numeric value, then the same number of modules
are created in both graphs.}
\item{If \code{nb_modul} is a numeric vector of length 2, then the
numbers of modules created in graphs \code{x} and \code{y} are the
first and second elements of \code{nb_modul}, respectively.}
}}

\item{algo}{(if x and y are igraph objects) A character string indicating the
algorithm used to create the modules with \pkg{igraph}.
\itemize{
\item{If \code{algo = 'fast_greedy'} (default),
function \code{cluster_fast_greedy}
from \pkg{igraph} is used (Clauset et al., 2004).}
\item{If \code{algo = 'walktrap'} (default), function \code{cluster_walktrap}
from \pkg{igraph} is used (Pons et Latapy, 2006) with
4 steps (default options).}
\item{If \code{algo = 'louvain'}, function \code{cluster_louvain}
from \pkg{igraph} is used (Blondel et al., 2008).
In that case, the number of modules created in each graph is imposed.}
\item{If \code{algo = 'optimal'}, function \code{cluster_optimal}
from \pkg{igraph} is used (Brandes et al., 2008) (can be very long).
In that case, the number of modules created in each graph is imposed.}
}}

\item{node_inter}{(optional, if x and y are igraph objects,
default is 'none') A character string indicating whether the links of the
graph are weighted by distances or by similarity indices. It is only used
to compute the modularity index. It can be: \itemize{
\item{'distance': Link weights correspond to distances. Nodes that are close
to each other will more likely be in the same module.}
\item{'similarity': Link weights correspond to similarity indices. Nodes that
are similar to each other will more likely be in the same module. Inverse
link weights are then used to compute the modularity index.}
\item{'none': Links are not weighted for the computation, which is only
based on graph topology.}
}
Two different weightings can be used to create the modules of the two graphs.
\itemize{
\item{If \code{node_inter} is a character string, then the same link
weighting is used for both graphs.}
\item{If \code{node_inter} is a character vector of length 2, then
the link weighting used by the algorithm to create the modules of
graphs \code{x} and \code{y} is determined by the first and second elements
of \code{node_inter}, respectively.}
}}

\item{data}{(if x and y are columns from a data.frame) An object of class
data.frame with at least two columns and as many rows as there are nodes
in the graphs compared. The columns indicate the modules of each node in
2 different classifications.}
}
\value{
The value of the ARI
}
\description{
The function computes the Adjusted Rand Index (ARI) to
compare two graphs' partitions into modules or clusters more generally.
Both graphs must have the same number of nodes, but not necessarily the same
number of links. They must also have the same node names and in the
same order.
}
\details{
This index takes values between -1 and 1. It measures how often
pairs of nodes pertaining to the same module in one graph also pertain to
the same module in the other graph.
Therefore, large values indicate that both partitions are similar.
The Rand Index can be defined as the frequency of agreement between two
classifications into discrete classes. It is the number of times a pair of
elements are classified into the same class or in two different classes
in both compared classifications, divided by the total number of possible
pairs of elements. The Rand Index is between 0 and 1 but its maximum value
depends on the number of elements. Thus, another 'adjusted' index was
created, the Adjusted Rand Index. According to the Hubert et
Arabie's formula, the ARI is computed as follows:
\eqn{ARI=\frac{Index - Expected index}{Maximum index - Expected index}}
where the values of Index, Expected index and Maximum index are computed
from a contingency table.
This function uses \code{adjustedRandIndex} from package \pkg{mclust} which
applies the Hubert and Arabie's formula for the ARI.
This function works for undirected graphs only.
}
\examples{
data(data_ex_genind)
data(pts_pop_ex)
mat_dist <- suppressWarnings(graph4lg::mat_geo_dist(data=pts_pop_ex,
      ID = "ID",
      x = "x",
      y = "y"))
mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
                      order(as.character(colnames(mat_dist)))]
graph_obs <- gen_graph_thr(mat_w = mat_dist, mat_thr = mat_dist,
                            thr = 24000, mode = "larger")
mat_gen <- mat_gen_dist(x = data_ex_genind, dist = "DPS")
graph_pred <- gen_graph_topo(mat_w = mat_gen, mat_topo = mat_dist,
                            topo = "gabriel")
ARI <- graph_modul_compar(x = graph_obs, y = graph_pred)
}
\references{
\insertRef{dyer2004population}{graph4lg}
\insertRef{hubert1985comparing}{graph4lg}
\insertRef{clauset2004finding}{graph4lg}
\insertRef{blondel2008fast}{graph4lg}
\insertRef{brandes2008modularity}{graph4lg}
\insertRef{pons2006computing}{graph4lg}
}
\author{
P. Savary
}
