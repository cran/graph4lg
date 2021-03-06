% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_w_hist.R
\name{plot_w_hist}
\alias{plot_w_hist}
\title{Plot histograms of link weights}
\usage{
plot_w_hist(graph, fill = "#396D35", class_width = NULL)
}
\arguments{
\item{graph}{A graph object of class \code{igraph} whose links are weighted}

\item{fill}{A character string indicating the color used to fill
the bars (default: "#396D35"). It must be a hexadecimal color code or
a color used by default in R.}

\item{class_width}{(default values: NULL) A numeric or an integer specifying
the width of the classes displayed on the histogram. When it is not
specified, the width is equal to the difference between the minimum and
maximum values divided by 80.}
}
\value{
A ggplot2 object to plot
}
\description{
The function enables to plot histogram to visualize the
distribution of the link weights
}
\examples{
data(data_ex_genind)
mat_w <- mat_gen_dist(data_ex_genind, dist = "DPS")
gp <- gen_graph_topo(mat_w = mat_w, topo = "gabriel")
hist <- plot_w_hist(gp)
}
\author{
P. Savary
}
