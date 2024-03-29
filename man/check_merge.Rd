% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_merge.R
\name{check_merge}
\alias{check_merge}
\title{Check whether the option 'nomerge' was used when building the landscape
graph with Graphab}
\usage{
check_merge(proj_end_path)
}
\arguments{
\item{proj_end_path}{The path to the project .xml file.}
}
\value{
The function returns a logical indicating whether 'nomerge' was used.
If nomerge=TRUE, then it returns FALSE. If nomerge=FALSE, it returns TRUE.
}
\description{
The function checks whether the option 'nomerge' was used when
building the landscape graph with Graphab
}
\examples{
\dontrun{
proj_name <- "grphb_ex"
check_merge(proj_name = proj_name)
}
}
\author{
P. Savary
}
\keyword{internal}
