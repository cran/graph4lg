#' Get metrics computed at the patch/node level in the Graphab project
#'
#' @description The function gets the values of metrics computed at the
#' patch/node level in the Graphab project.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is.
#' @param habitat (default=NULL) A character string indicating the
#' name of the habitat type for which the metric values should be returned.
#' If NULL, there is no selection on the habitat type.
#' @param graph (default=NULL) A character string indicating the
#' name of the graph for which the metric values should be returned.
#' If NULL, there is no selection on the graph name.
#' @param metric (default=NULL) A character string indicating the
#' name of the metric for which the values should be returned.
#' If NULL, there is no selection on the metric name.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When \code{proj_path = NULL}, the project directory is equal
#' to \code{getwd()}.
#' @return A data.frame or list of data.frames with the values of metrics
#' computed at the patch/node level for a given habitat type or a given graph.
#' The data.frame indicate the type of habitat to which each node/patch
#' corresponds, and the name of the metric includes that of the graph on which
#' it was computed.
#' @details Only local metrics computed for patches/nodes, for a given habitat
#' type and from one or several graphs can be returned. Global metrics are not
#' allowed. When the name of a metric is given, only the values of this metric
#' are returned. At least one of the parameters \code{habitat}, \code{habitat},
#' or \code{metric} must not be NULL.
#' See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' get_graphab_metric(proj_name = "graphab_example",
#'                    habitat = "forest")
#' }

get_graphab_metric <- function(proj_name, # character
                               habitat = NULL,
                               graph = NULL,
                               metric = NULL,
                               proj_path = NULL){ # if null getwd() otherwise a character path

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

  proj_end_path <- paste0(proj_path, "/", proj_name, "/", proj_name, ".xml")

  ## Check project version
  check_graphab_version(proj_path = proj_end_path)

  #########################################
  # Get the project information
  project_info <- graphab_project_desc(proj_name = proj_name,
                                       proj_path = proj_path)
  project_habitat <- project_info[["Habitats"]][[1]]
  project_graph <- project_info[["Graphs"]]
  project_metric <- project_info[["Metrics"]]

  ## Stop the function if no metric exists.
  if(is.null(project_metric)){
    stop("No metric has been computed in this project yet.")
  }

  #########################################
  ## Check that at least one argument is not NULL and that they are correct
  if(all(c(is.null(habitat),
           is.null(graph),
           is.null(metric)))){
    stop(paste0("At least one of the 'habitat', 'graph', and 'metric' ",
                "parameters must be different from NULL."))
  }

  if(!is.null(habitat)){ ## check habitat
    if(!inherits(habitat, "character")){
      stop("'habitat' must be either NULL or a character string.")
    } else if (length(habitat) != 1){
      stop("You must provide a single habitat name.")
    } else if (!check_graphab_object(proj_path = proj_end_path,
                                     object_type = "habitat",
                                     name = habitat)){
      stop("The habitat you refer to does not exist.
           Please use graphab_habitat() before.")
    } else if(!file.exists(paste0(proj_path, "/",
                                  proj_name, "/",
                                  habitat, "/patches.gpkg"))){
      stop(paste0("The habitat type '", habitat, "' does not include any patch. ",
                  "Metrics cannot be found."))
    }
  }

  if(!is.null(graph)){ ## check graph
    if(!inherits(graph, "character")){
      stop("'graph' must be either NULL or a character string.")
    } else if (length(graph) != 1){
      stop("You must provide a single graph name.")
    } else if (!check_graphab_object(proj_path = proj_end_path,
                                     object_type = "graph",
                                     name = graph)){
      stop("The graph you refer to does not exist.
           Please use graphab_graph() before.")
    }
  }

  if(!is.null(metric)){ ## check metric
    if(!inherits(metric, "character")){
      stop("'metric' must be either NULL or a character string.")
    } else if (length(metric) != 1){
      stop("You must provide a single metric name.")
    } else if (!(metric %in% project_metric$metric)){
      stop("The metric you refer to does not exist.
           Please use graphab_metric() before.")
    }
  }

  #########################################
  ## Cover the different cases, starting with metric, then graph and habitat

  # If metric is provided
  ##########################
  if(!is.null(metric)){

    # Warning that other arguments are ignored.
    if(any(c(!is.null(graph), !is.null(habitat)))){
      warning("When 'metric' is provided, 'habitat' and 'graph' are ignored.")
    }

    # Get the metric short name (without graph) and the graph name
    metric_graph <- project_metric[which(project_metric$metric == metric),
                                   "graph"]
    metric_short <- stringr::str_sub(metric, 1, -(nchar(metric_graph)+2))

    # Open all the 'patches' layer and search for the metric values
    list_hab_df <- list()
    for(i in 1:length(project_habitat)){

      # Get colnames of the relevant layer
      col_hab_i <- colnames(
        suppressWarnings(sf::read_sf(
          paste0(proj_path, "/", proj_name, "/",
                 project_habitat[i], "/patches.gpkg"),
          query = "SELECT * FROM patches LIMIT 0",
          as_tibble = FALSE)))

      detect_test <- c(stringr::str_detect(string = col_hab_i,
                                           pattern = metric_graph) &
                         stringr::str_detect(string = col_hab_i,
                                             pattern = metric_short))

      if(any(detect_test)){

        hab_i <- suppressWarnings(
          sf::st_drop_geometry(
            sf::read_sf(
              paste0(proj_path, "/", proj_name, "/",
                     project_habitat[i], "/patches.gpkg"),
              as_tibble = FALSE)))

        hab_i <- hab_i[, which(colnames(hab_i) %in%
                                 c("idhab", "Id", "area",
                                   "perim", "capacity",
                                   col_hab_i[which(detect_test)]))]

        hab_i$habitat <- project_habitat[i]

        hab_i <- hab_i[, c("habitat", "idhab", "Id",
                           "area", "perim", "capacity",
                           col_hab_i[which(detect_test)])]

        colnames(hab_i) <- c("habitat", "id_habitat", "id_patch",
                             "area", "perim", "capacity",
                             col_hab_i[which(detect_test)])

        list_hab_df[[i]] <- hab_i
      }
    }

    if(length(list_hab_df) == 0){
      stop(paste0("You provided a global metric but this function only ",
                  "returns local metric values."))
    } else {
      # Keep only the non-NULL list elements
      list_hab_df <- list_hab_df[!(unlist(lapply(list_hab_df, is.null)))]
      # Stack the table keeping their habitat name
      res <- do.call("rbind", list_hab_df)
    }

    ##########################
  } else if(!is.null(graph)){

    # Warning that habitat argument is ignored.
    if(!is.null(habitat)){
      warning("When 'graph' is provided, 'habitat' is ignored.")
    }

    # Check that metrics have been comptued on the given graph
    if(!any(project_metric$graph == graph)){
      stop(paste0("No metric has been computed on the graph '", graph, "'."))
    }

    # Open all the 'patches' layer and search for the graph name in columns
    list_hab_df <- list()
    for(i in 1:length(project_habitat)){

      # Get colnames of the relevant layer
      col_hab_i <- colnames(
        suppressWarnings(sf::read_sf(
          paste0(proj_path, "/", proj_name, "/",
                 project_habitat[i], "/patches.gpkg"),
          query = "SELECT * FROM patches LIMIT 0",
          as_tibble = FALSE)))

      # Search for the graph name in the metric names
      to_search_in <- stringr::str_sub(col_hab_i, -nchar(graph), -1)
      detect_test <- !(is.na(stringr::str_match(string = to_search_in,
                                                pattern = graph)))

      if(any(detect_test)){

        hab_i <- suppressWarnings(
          sf::st_drop_geometry(
            sf::read_sf(
              paste0(proj_path, "/", proj_name, "/",
                     project_habitat[i], "/patches.gpkg"),
              as_tibble = FALSE)))

        hab_i <- hab_i[, which(colnames(hab_i) %in%
                                 c("idhab", "Id", "area",
                                   "perim", "capacity",
                                   col_hab_i[which(detect_test)]))]

        hab_i$habitat <- project_habitat[i]

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

    ##########################
  } else { ### For habitat types

    # Open the provided habitat table
    hab_df <- suppressWarnings(
      sf::st_drop_geometry(
        sf::read_sf(
          paste0(proj_path, "/", proj_name, "/",
                 habitat, "/patches.gpkg"),
          as_tibble = FALSE)))

    hab_df$habitat <- habitat

    hab_df <- hab_df[, c("habitat", "idhab", "Id",
                         "area", "perim", "capacity",
                         colnames(hab_df)[6:(ncol(hab_df)-1)])]
    colnames(hab_df)[1:6] <- c("habitat", "id_habitat", "id_patch",
                               "area", "perim", "capacity")

    res <- hab_df
  }
  #########################################
  return(res)
}
