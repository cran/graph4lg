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
#' \item{Equivalent Connectivity (\code{metric = 'EC'}): Square root of the
#' sum of products of capacity of all pairs of patches weighted by their
#' interaction probability. This is the size of a single patch (maximally
#' connected) that would provide the same probability of connectivity as the
#' actual habitat pattern in the landscape (Saura et al., 2011).}
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
#' @param multihab (optional, default = NULL) A character string indicating
#' whether the metric value should be decomposed across the different habitat
#' types and their pairwise combinations. If `multihab='all'`, all the
#' pairwise combinations are considered, including the within-habitat case.
#' If `multihab='inter'`, only the inter-habitat types combinations are
#' considered. Note that this argument is only required if the graph on which
#' you compute the metrics is based on several habitat types.
#' @param resfile (optional, default = NULL) A character string giving the
#' name of the text file storing the results of the computation. Must be
#' of the form 'file.txt'.
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
#' are returned in R (TRUE) or only stored in the patch attribute layers (FALSE)
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @param parallel.java An integer indicating how many computer cores are used
#' to run the .jar file. By default, \code{parallel.java = NULL}, and java sets
#' it according to local settings.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @return If \code{return_val=TRUE}, the function returns a \code{data.frame}
#' with the computed metric values and the corresponding patch ID when the
#' metric is local or delta metric, or the numeric value of the global metric.
#' If the metric is computed on a graph whose nodes belong to different habitat
#' categories, the returned table includes several habitat categories.
#' @details The metrics are described in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' Graphab software makes possible the computation of other metrics.
#' Be careful, when the same metric is computed several times, the option
#' \code{return=TRUE} is not returning the right columns. In these cases,
#' use \code{\link{get_graphab_metric}}.
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' graphab_metric(proj_name = "graphab_example",
#'                graph = "graph_forest",
#'                metric = "F",
#'                multihab = "all",
#'                dist = 1000,
#'                prob = 0.05,
#'                beta = 1,
#'                cost_conv = TRUE)
#' }

graphab_metric <- function(proj_name, # character
                           graph, # cost or euclid
                           metric, # character
                           multihab = NULL, # NULL, 'all' or 'inter'
                           resfile = NULL,
                           dist = NULL, # dist threshold
                           prob = 0.05, # dispersal probability
                           beta = 1, # area weight
                           cost_conv = FALSE, # FALSE (default) or true
                           return_val = TRUE, # return the metric values
                           proj_path = NULL, # if null getwd() otherwise a character path
                           parallel.java = NULL,
                           alloc_ram = NULL){

  #########################################
  # Check for project directory path
  if(!is.null(proj_path)){
    if(!dir.exists(proj_path)){
      stop(paste0(proj_path, " is not an existing directory or the path is ",
                  "incorrectly specified."))
    } else {
      proj_path <- normalizePath(proj_path)
    }
  } else {
    proj_path <- normalizePath(getwd())
  }

  #########################################
  # Check for proj_name class
  if(!inherits(proj_name, "character")){
    stop("'proj_name' must be a character string")
  } else if (!(paste0(proj_name, ".xml") %in%
               list.files(path = paste0(proj_path, "/", proj_name)))){
    stop("The project you refer to does not exist.
         Please use graphab_project() before.")
  }

  ## Create proj_end_path proj_path/proj_name/proj_name.xml
  proj_end_path <- paste0(proj_path, "/", proj_name, "/", proj_name, ".xml")

  #######################################
  # Add '' to proj_path for cases with spaces in paths
  if(all(stringr::str_sub(proj_end_path, 1, 1) != "'",
         stringr::str_sub(proj_end_path, 1, 1) != "'",
         stringr::str_detect(string = proj_end_path,
                             pattern = " "))){
    proj_end_path_cmd <- paste0("'", proj_end_path, "'")
  } else {
    proj_end_path_cmd <- proj_end_path
  }

  ## Check project version
  check_graphab_version(proj_path = proj_end_path)

  ## List existing metrics
  graphab_objects <- graphab_show(proj_path = proj_end_path)
  old_metrics <- graphab_objects[["Metrics"]]

  #########################################
  # Check for graph class
  if(!inherits(graph, "character")){
    stop("'graph' must be a character string")
  } else if (!check_graphab_object(proj_path = proj_end_path,
                                   object_type = "graph",
                                   name = graph)){
    stop("The graph you refer to does not exist.
           Please use graphab_graph() to create it.")
  }

  ### Check for multihab
  if(!is.null(multihab)){
    if(!inherits(multihab, "character")){
      stop("'multihab' must be a character string")
    }
  }

  ### Check for resfile
  if(!is.null(resfile)){

    if(!inherits(resfile, "character")){
      stop("'resfile' must be a character string")
    } else if(stringr::str_sub(resfile, -4, -1) != ".txt") {
      stop("'resfile' must be a text file name such as 'file.txt'.")
    }

    ## Check whether resfile already exists
    if(file.exists(paste0(proj_path, "/", proj_name,
                          "/", resfile))){
      warning(paste0("The file '", resfile, "' already exists and ",
                     "will be overwritten."))
    }

  }

  #########################################
  # Check for metric and parameters

  list_all_metrics <- c("PC", "EC", "IIC", "dPC",
                        "F", "BC", "IF", "Dg", "CCe", "CF")

  list_glob_metrics <- list_all_metrics[1:4]
  list_loc_metrics <- list_all_metrics[5:length(list_all_metrics)]

  list_dist_metrics <- c("PC", "EC",
                         "F", "BC",
                         "IF", "dPC")

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

  ### Metrics not available with multiple habitats
  if(metric %in% c("dPC", "IIC", "PC",
                   "CF", "Dg", "CCe")){
    if(!is.null(multihab)){
      stop("The metric is not available in the multiple habitat mode.")
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
    stop(paste0("'metric' must be ", paste(list_all_metrics,
                                           collapse = " or ")))
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
  # Check for parallel.java
  if(!is.null(parallel.java)){
    if(!inherits(parallel.java, c("numeric", "integer"))){
      stop("'parallel.java' must be a numeric or integer value.")
    }
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
  version <- "graphab-3.0.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(),
                            "/graph4lg2_jar/", version)

  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

  if(!is.null(parallel.java)){
    cmd <- c(cmd, "-proc ", as.character(parallel.java))
  }

  cmd <- c(cmd,
           "--project", proj_end_path_cmd,
           "--usegraph", graph)

  if(level == "graph"){
    cmd <- c(cmd, "--gmetric", metric)
    if(!is.null(resfile)){
      cmd <- c(cmd, paste0("resfile=", resfile))
    }
  } else if (level == "patch"){
    cmd <- c(cmd, "--lmetric", metric)
  }

  if(!is.null(multihab)){
    cmd <- c(cmd, paste0("mh=", multihab))
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
    cmd[which(cmd == "dPC")] <- "PC"
    cmd <- c(cmd, "--delta", "obj=patch")
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
    }
  }


  ## Check whether a new metric exists
  graphab_objects <- graphab_show(proj_path = proj_end_path)
  new_metrics <- graphab_objects[["Metrics"]]

  if(length(new_metrics) == length(old_metrics)){

    warning("An error occurred or the metric already existed.")

    if(return_val){
      if(any(level == "patch",
             metric == "dPC")){
       stop(paste0("You cannot return the value of an already existing metric ",
                   "with this function. Please use 'get_graphab_metric()."))
      } else if(level == "graph"){
        new_metric <- paste0(metric, "_", graph) # to generate an output
      }
    }

  } else {
    new_metric <- new_metrics[!(new_metrics %in% old_metrics)]
    name_new_metric <- stringr::str_sub(new_metric,
                                        1,
                                        nchar(new_metric) - nchar(graph) - 1)
    message(paste0("Metric '", new_metric,
                   "' has been computed in the project '",
                   proj_name, "'."))

  }

  ## Get the metric value
  if(return_val){

    if(level == "patch"){

      # Check the habitats
      existing_hab <- graphab_objects[["Habitats"]]
      split_hab <- stringr::str_split(existing_hab, " - ")

      hab_codes <- unlist(lapply(split_hab, "[", 1))
      hab_names <- unlist(lapply(split_hab, "[", 2))

      # Open all the 'patches' layer and search for the new metric values
      list_hab_df <- list()
      for(i in 1:length(hab_names)){

        # Get colnames of the relevant layer
        col_hab_i <- colnames(
            suppressWarnings(sf::read_sf(
              paste0(proj_path, "/", proj_name, "/",
                     hab_names[i], "/patches.gpkg"),
              query = "SELECT * FROM patches LIMIT 0",
              as_tibble = FALSE)))

        detect_test <- c(stringr::str_detect(string = col_hab_i,
                                             pattern = graph) &
                           stringr::str_detect(string = col_hab_i,
                                               pattern = name_new_metric))

        if(any(detect_test)){

          hab_i <- suppressWarnings(sf::st_drop_geometry(sf::read_sf(
              paste0(proj_path, "/", proj_name, "/",
                     hab_names[i], "/patches.gpkg"),
              as_tibble = FALSE)))

          hab_i <- hab_i[, which(colnames(hab_i) %in%
                                   c("idhab", "Id", "area",
                                     "perim", "capacity",
                                     col_hab_i[which(detect_test)]))]

          hab_i$habitat <- hab_names[i]

          hab_i <- hab_i[, c("habitat", "idhab", "Id",
                             "area", "perim", "capacity",
                             col_hab_i[which(detect_test)])]

          colnames(hab_i) <- c("habitat", "id_habitat", "id_patch",
                               "area", "perim", "capacity",
                               col_hab_i[which(detect_test)])

          list_hab_df[[i]] <- hab_i
        }
      }
      # Keep only the non-NULL list elements
      list_hab_df <- list_hab_df[!(unlist(lapply(list_hab_df, is.null)))]
      # Stack the table keeping their habitat name
      res <- do.call("rbind", list_hab_df)

    } else if(level == "graph"){

      if(metric == "dPC"){

        new_metric1 <- new_metric[1]
        name_new_metric1 <- name_new_metric[1]

        # Check the habitats
        existing_hab <- graphab_objects[["Habitats"]]
        split_hab <- stringr::str_split(existing_hab, " - ")

        hab_codes <- unlist(lapply(split_hab, "[", 1))
        hab_names <- unlist(lapply(split_hab, "[", 2))

        # Open all the 'patches' layer and search for the new metric values
        list_hab_df <- list()
        for(i in 1:length(hab_names)){

          # Get colnames of the relevant layer
          col_hab_i <- colnames(
              suppressWarnings(sf::read_sf(
                paste0(proj_path, "/", proj_name, "/",
                       hab_names[i], "/patches.gpkg"),
                query = "SELECT * FROM patches LIMIT 0",
                as_tibble = FALSE)))

          detect_test <- c(stringr::str_detect(string = col_hab_i,
                                               pattern = graph) &
                             stringr::str_detect(string = col_hab_i,
                                                 pattern = name_new_metric1))

          if(any(detect_test)){

            hab_i <- suppressWarnings(sf::st_drop_geometry(sf::read_sf(
                paste0(proj_path, "/", proj_name, "/",
                       hab_names[i], "/patches.gpkg"),
                as_tibble = FALSE)))

            hab_i <- hab_i[, which(colnames(hab_i) %in%
                                     c("idhab", "Id", "area",
                                       "perim", "capacity",
                                       col_hab_i[which(detect_test)]))]

            hab_i$habitat <- hab_names[i]

            hab_i <- hab_i[, c("habitat", "idhab", "Id",
                               "area", "perim", "capacity",
                               col_hab_i[which(detect_test)])]

            colnames(hab_i) <- c("habitat", "id_habitat", "id_patch",
                                 "area", "perim", "capacity",
                                 col_hab_i[which(detect_test)])

            list_hab_df[[i]] <- hab_i
          }
        }

        # Keep only the non-NULL list elements
        list_hab_df <- list_hab_df[!(unlist(lapply(list_hab_df, is.null)))]
        # Stack the table keeping their habitat name
        res <- do.call("rbind", list_hab_df)

      } else {

        if(!is.null(resfile)){
          # If resfile is specified, open the text file
          res_val <- utils::read.table(file = paste0(proj_path, "/",
                                                     proj_name, "/",
                                                     resfile),
                                       header = TRUE)
          # List of metric name and values
          res <- list(new_metric,
                      res_val)
          names(res) <- c("Metric name", "Metric value table")
        } else {
          # If resfile is NULL, the text file is named as the metric
          res_val <- utils::read.table(file = paste0(proj_path, "/",
                                                     proj_name, "/",
                                                     metric, ".txt"),
                                       header = TRUE)
          # List of metric name and values
          res <- list(new_metric,
                      res_val)
          names(res) <- c("Metric name", "Metric value table")
        }
      }
    }
    return(res)
  }
}
