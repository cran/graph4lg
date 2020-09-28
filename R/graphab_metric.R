#' Compute connectivity metrics from a graph in the Graphab project
#'
#' @description The function computes connectivity metrics on a graph from a
#' link set in a Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is.
#' @param graph A character string indicating the name of the graph on which
#' the metric is computed. This graph has been created with Graphab
#' or using \code{\link{graphab_graph}} function and is associated
#' with a link set. Only the links present in the graph and their corresponding
#' weights will be used in the computation, together with patch areas.
#' @param metric A character string indicating the metric which will be computed
#' on the graph. This metric can be:\itemize{
#' \item{A global metric:\itemize{
#' \item{Probability of Connectivity (\code{metric = 'PC'}): Sum of products of
#' area of all pairs of patches weighted by their interaction probability,
#' divided by the square of the area of the study zone.
#' This ratio is the equivalent to the probability that two points randomly
#' placed in the study area are connected.}
#' \item{Integral Index of Connectivity (\code{metric = 'IIC'}): For the
#' entire graph: product of patch areas divided by the number of links
#' between them, the sum is divided by the square of the area of the study
#' zone. IIC is built like the PC index but using the inverse of a topological
#' distance rather than a negative exponential function of the distance
#' based on the link weight.}
#' }}
#' \item{A local metric:\itemize{
#' \item{Flux (\code{metric = 'F'}): For the focal patch i : sum of area
#' of patches other than i and weighted according to their minimum distance
#' to the focal patch through the graph. This sum is an indicator of the
#' potential dispersion from the patch i or, conversely to the patch i}
#' \item{Betweenness Centrality index (\code{metric = 'BC'}): Sum of the
#' shortest paths through the focal patch i, each path is weighted by the
#' product of the areas of the patches connected and of their interaction
#' probability. All possible paths between every pair of patches is
#' considered in this computation.}
#' \item{Interaction Flux (\code{metric = 'IF'}): Sum of products of the focal
#' patch area with all the other patches, weighted by their interaction
#' probability.}
#' \item{Degree (\code{metric = 'Dg'}): Number of edges connected to the
#' node i i.e. number of patches connected directly to the patch i.}
#' \item{Closeness Centrality index (\code{metric = 'CCe'}): Mean distance
#' from the patch i to all other patches of its component k.}
#' \item{Current Flux (\code{metric = 'CF'}): Sum of currents passing through
#' the patch i. \eqn{c_{i}^{j}} represents the current through the patch i when
#' currents are sent from all patches (except j) to the patch j.
#' The patch j is connected to the ground.}
#' }}
#' \item{A delta metric:\itemize{
#' \item{delta Probability of Connectivity (\code{metric = 'dPC'}): Rate of
#' variation between the value of PC index and the value of PC' corresponding
#' to the removal of the patch i. The value of \code{dPC} is decomposed
#' into three parts:\itemize{
#' \item{\eqn{dPC_{area}} is the variation induced by the area lost after removal;}
#' \item{\eqn{dPC_{flux}} is the variation induced by the loss of interaction
#' between the patch i and other patches;}
#' \item{\eqn{dPC_{connector}} is the variation induced by the modification of
#' paths connecting other patches and initially routed through i.}
#' }
#' }}
#' }}
#' For most metrics, the interaction probability is computed for each pair of
#' patches from the path that minimizes the distance d (or the cost) between
#' them. It then maximizes \eqn{{e}^{-\alpha d_{ij}}} for patches i and j.
#' To use patch capacity values different from the patch area, please use
#' directly Graphab software.
#' @param dist A numeric or integer value specifying the distance at which
#' dispersal probability is equal to \code{prob}. This argument is mandatory
#' for weighted metrics (PC, F, IF, BC, dPC, CCe, CF) but not used for others.
#' It is used to set \eqn{\alpha} for computing dispersal probabilities associated
#' with all inter-patch distances such that dispersal probability between
#' patches i and j is \eqn{p_{ij}= e^{-\alpha d_{ij}}}.
#' @param prob A numeric or integer value specifying the dispersal probability
#' at distance \code{dist}. By default, \code{code=0.05}. It is used to set
#' \eqn{\alpha} (see param \code{dist} above).
#' @param beta A numeric or integer value between 0 and 1 specifying the
#' exponent associated with patch areas in the computation of metrics
#' weighted by patch area. By default, \code{beta=1}. When \code{beta=0}, patch
#' areas do not have any influence in the computation.
#' @param cost_conv FALSE (default) or TRUE. Logical indicating whether numeric
#' \code{dist} values are converted from cost-distance into Euclidean distance
#' using a log-log linear regression. See also \code{\link{convert_cd}}
#' function.
#' @param return_val Logical (default = TRUE) indicating whether metric values
#' are returned in R (TRUE) or only stored in the patch attribute layer (FALSE)
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @return If \code{return_val=TRUE}, the function returns a \code{data.frame}
#' with the computed metric values and the corresponding patch ID when the
#' metric is local or delta metric, or the numeric value of the global metric.
#' @details The metrics are described in Graphab 2.4 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.4-en.pdf}
#' Graphab software makes possible the computation of other metrics.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_metric(proj_name = "grphb_ex",
#'                graph = "graph",
#'                metric = "PC",
#'                dist = 1000,
#'                prob = 0.05,
#'                beta = 1)
#' }



graphab_metric <- function(proj_name, # character
                           graph, # cost or euclid
                           metric, # character
                           dist = NULL, # dist threshold
                           prob = 0.05, # dispersal probability
                           beta = 1, # area weight
                           cost_conv = FALSE, # FALSE (default) or true
                           return_val = TRUE, #
                           proj_path = NULL, # if null getwd() otherwise a character path
                           alloc_ram = NULL){

  #########################################
  # Check for project directory path
  if(!is.null(proj_path)){
    chg <- 1
    wd1 <- getwd()
    setwd(dir = proj_path)
  } else {
    chg <- 0
    proj_path <- getwd()
  }

  #########################################
  # Check for proj_name class
  if(!inherits(proj_name, "character")){
    stop("'proj_name' must be a character string")
  } else if (!(paste0(proj_name, ".xml") %in% list.files(path = paste0("./", proj_name)))){
    stop("The project you refer to does not exist.
         Please use graphab_project() before.")
  }

  proj_end_path <- paste0(proj_name, "/", proj_name, ".xml")

  #########################################
  # Check for graph class
  if(!inherits(graph, "character")){
    stop("'graph' must be a character string")
  } else if (!(paste0(graph, "-voronoi.shp") %in% list.files(path = paste0("./", proj_name)))){
    stop("The graph you refer to does not exist")
  } else if (length(list.files(path = paste0("./", proj_name), pattern = "-voronoi.shp")) == 0){
    stop("There is not any graph in the project you refer to.
         Please use graphab_graph() before.")
  }

  #########################################
  # Check for metric and parameters

  list_all_metrics <- c("PC", "IIC", "dPC", "F", "BC", "IF", "Dg", "CCe", "CF")

  list_glob_metrics <- list_all_metrics[1:3]
  list_loc_metrics <- list_all_metrics[4:length(list_all_metrics)]

  list_dist_metrics <- c("PC", "F", "BC", "IF", "dPC")

  if(metric %in% list_dist_metrics){
    if(is.null(dist)){
      stop(paste0("To compute ", metric, ", specify a distance associated to
                  a dispersal probability (default=0.05)"))
    } else if(!inherits(dist, c("numeric", "integer"))){
      stop("'dist' must be a numeric or integer value")
    } else if(!inherits(prob, c("numeric", "integer"))){
      stop("'prob' must be a numeric or integer value")
    } else if(!inherits(beta, c("numeric", "integer"))){
      stop("'beta' must be a numeric or integer value")
    } else if(beta < 0 || beta > 1){
      stop("'beta' must be between 0 and 1")
    } else if(prob < 0 || prob > 1){
      stop("'prob' must be between 0 and 1")
    }
  }

  ### Special case with dPC
  if(metric == "dPC"){
    if(cost_conv){
      stop("Option 'cost_conv = TRUE' is not available with the metric dPC")
    }
  }

  # Special case of CF with beta
  if(metric == "CF"){
    if(beta < 0 || beta > 1){
      stop("'beta' must be between 0 and 1")
    }
  }

  if(!inherits(metric, "character")){
    stop("'metric' must be a character string")
  } else if(!(metric %in% list_all_metrics)){
    stop(paste0("'metric' must be ", paste(list_all_metrics, collapse = " or ")))
  } else if(metric %in% list_loc_metrics){
    level <- "patch"
  } else if(metric %in% list_glob_metrics){
    level <- "graph"
  }

  #########################################
  # Check for cost_conv
  if(!is.logical(cost_conv)){
    stop("'cost_conv' must be a logical (TRUE or FALSE).")
  }

  #########################################
  # Check for return_val
  if(!is.logical(return_val)){
    stop("'return_val' must be a logical (TRUE or FALSE).")
  }

  #########################################
  # Check for Graphab
  gr <- get_graphab(res = FALSE, return = TRUE)

  if(gr == 1){
    message("Graphab has been downloaded")
  }

  #########################################
  # Get java path
  java.path <- Sys.which("java")

  #########################################
  # Get graphab path
  version <- "graphab-2.4.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)

  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
           "--project", proj_end_path,
           "--usegraph", graph)

  if(level == "graph"){
    cmd <- c(cmd, "--gmetric", metric)
  } else if (level == "patch"){
    cmd <- c(cmd, "--lmetric", metric)
  }

  if(metric %in% list_dist_metrics){
    if(cost_conv){
      cmd <- c(cmd,
               paste0("d={", dist, "}"),
               paste0("p=", prob),
               paste0("beta=", beta))
    } else {
      cmd <- c(cmd,
               paste0("d=", dist),
               paste0("p=", prob),
               paste0("beta=", beta))
    }
  } else if (metric == "CF"){
    cmd <- c(cmd,
             paste0("beta=", beta))

  }

  if(metric == "dPC"){
    cmd[8] <- "--delta"
    cmd[13] <- "obj=patch"
  }

  if(!is.null(alloc_ram)){
    if(inherits(alloc_ram, c("integer", "numeric"))){
      cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
    } else {
      stop("'alloc_ram' must be a numeric or an integer")
    }
  }

  #########################################
  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)

  if(length(rs) == 1){
    if(rs == 1){
      message("An error occurred")
    } else {
      message(paste0("Metric '", metric, "' has been computed in the project ",
                     proj_name))
    }
  } else {
    message(paste0("Metric '", metric, "' has been computed in the project ",
                   proj_name))
  }

  if(return_val){
    if(level == "graph"){

      if(metric != "dPC"){

        val <- as.numeric(stringr::str_split(rs[2], pattern = ":")[[1]][2])
        vec_res <- c(paste0("Project : ", proj_name),
                     paste0("Graph : ", graph),
                     paste0("Metric : ", metric))

        if(metric %in% list_dist_metrics){
          vec_res <- c(vec_res,
                       paste0("Dist : ", dist),
                       paste0("Prob : ", prob),
                       paste0("Beta : ", beta))
        }

        vec_res <- c(vec_res,
                     paste0("Value : ", val))

        #print(paste(vec_res, collapse = ", "))
        res <- vec_res
      } else if (metric == "dPC") {

        #   fdpc <- list.files(path = proj_name, pattern = "delta-dPC")
        #   file.info(paste0("./", proj_name, "/", fdpc))[, 'mtime']
        #   file.info(paste0("./", proj_name, "/",
        #                    name_txt))

        name_txt <- paste0("delta-dPC_d", dist, "_p", prob, "_", graph, ".txt")
        res_dpc <- utils::read.table(file = paste0("./", proj_name, "/",
                                            name_txt),
                              header = TRUE)[-1, ]

        vec_res <- c(paste0("Project : ", proj_name),
                     paste0("Graph : ", graph),
                     paste0("Metric : ", metric),
                     paste0("Dist : ", dist),
                     paste0("Prob : ", prob),
                     paste0("Beta : ", beta))

        res <- list(vec_res, res_dpc)

      }

    } else if (level == "patch"){

      name_met <- gsub(stringr::str_split(rs[3], pattern = ":")[[1]][2],
                       pattern = " ", replacement = "")

      tab <- utils::read.csv(file = paste0(proj_name, "/patches.csv"))

      df_res <- tab[, c(1,
                        which(colnames(tab) == paste0(name_met, "_", graph)))]

      vec_res <- c(paste0("Project : ", proj_name),
                   paste0("Graph : ", graph),
                   paste0("Metric : ", metric))

      if(metric %in% list_dist_metrics){
        vec_res <- c(vec_res,
                     paste0("Dist : ", dist),
                     paste0("Prob : ", prob),
                     paste0("Beta : ", beta))
      } else if (metric == "CF"){
        vec_res <- c(vec_res,
                     paste0("Beta : ", beta))
      }

      res <- list(vec_res, df_res)
    }
    return(res)
  }

  #########################################
  if(chg == 1){
    setwd(dir = wd1)
  }

}
