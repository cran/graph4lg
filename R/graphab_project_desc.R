#' Describe the objects of a Graphab project
#'
#' @description The function describes the objects of a Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml has been created.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @return It returns a large list describing the objects in several categories,
#' including:\itemize{
#' \item{Source raster: the values of the source raster file.}
#' \item{Habitats: the list of created habitats, their code, their name,
#' their type (raster, vector), and, if raster, their code in the raster.}
#' \item{Linksets: the list of created linksets and their type. If computed
#' in cost distance units, a table of costs or the name of the external cost
#' raster file is returned.}
#' \item{Graphs: the list of graphs that have been created.}
#' \item{Metrics: the list of computed metrics and their associated graphs.}
#' }
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_project_desc(proj_name = "graphab_example")
#' }


graphab_project_desc <- function(proj_name,
                                 proj_path = NULL){


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

  ### Source raster codes
  raster_codes <- get_graphab_raster_codes(proj_name = proj_name,
                                           mode = 'all',
                                           proj_path = proj_path)

  ### List objects with graphab_show
  graphab_objects <- graphab_show(proj_path = proj_end_path)

  # Get the habitats
  existing_hab <- graphab_objects[["Habitats"]]
  split_hab <- stringr::str_split(existing_hab, " - ")
  hab_codes <- unlist(lapply(split_hab, "[", 1))
  hab_names <- unlist(lapply(split_hab, "[", 2))

  # Get the linksets
  existing_link <- graphab_objects[["Linksets"]]

  # Get the graphs
  existing_graph <- graphab_objects[["Graphs"]]

  # Get the metrics
  existing_metric <- graphab_objects[["Metrics"]]

  #########################################################
  # Copy the .xml file as a .txt file in temp files to open it
  xml <- tempfile(pattern = ".txt")
  file.copy(from = proj_end_path,
            to = xml)
  file_data <- utils::read.table(xml)

  #### Find information about each object type and summarize it

  # For habitats
  if(!is.null(existing_hab)){

    # where it starts and ends in the .xml file
    type_begin <- "<habitats>"
    type_end <- "</habitats>"

    # Get object type lines
    first_line <- min(which(file_data[, 1] == type_begin))
    last_line <- min(which(file_data[, 1] == type_end))

    hab_df <- data.frame(id_hab = NA, name = NA,
                         type = NA, rast_codes = NA)
    hab_df <- hab_df[-1, ]

    for(i in 1:length(existing_hab)){

      target <- paste0("<string>", hab_names[i], "</string>")
      # Search for the target name in these lines
      line_target <- which(file_data[first_line:last_line,
                                     1] == target)

      # line_target + first_line - 1 (start) + 1: line describing the source
      if(file_data[first_line + line_target, 1] == "<Habitat>"){

        # Check whether several patch codes are given and return them
        # get the limits of this habitat specification
        first_hab_line <- first_line + line_target
        last_hab_line <- min(which(file_data[first_hab_line:nrow(file_data),
                                             1] == "</Habitat>")) + first_hab_line - 1

        # find the patch codes block

        begin_patch_codes <- first_hab_line +
          which(file_data[(first_hab_line+1):(last_hab_line-1),
                          1] == "<patchCodes>")
        end_patch_codes <- first_hab_line +
          which(file_data[(first_hab_line+1):(last_hab_line-1),
                          1] == "</patchCodes>")

        lines_patch_codes <- seq(begin_patch_codes + 1,
                                 end_patch_codes - 1,
                                 1)

        rast_codes <- unlist(lapply(file_data[lines_patch_codes, 1],
                                    FUN = function(x){
                                      x <- stringr::str_sub(x, 6, -7)
                                    }))

        rast_codes <- paste(rast_codes, collapse = ", ")

        hab_df <- rbind(hab_df,
                        data.frame(id_hab = hab_codes[i],
                                   name = hab_names[i],
                                   type = "Raster",
                                   rast_codes = rast_codes))

      } else if(stringr::str_detect(file_data[first_line + line_target, 1],
                                    pattern = "VectorHabitat>")){
        hab_df <- rbind(hab_df,
                        data.frame(id_hab = hab_codes[i],
                                   name = hab_names[i],
                                   type = "Vector",
                                   rast_codes = NA))
      }
    }

    # Reorder by habitat ID
    hab_df <- hab_df[order(hab_df$id_hab), ]

    # Make a result list
    hab_list <- list(hab_names, hab_df)
    names(hab_list) <- c("Habitat names", "Habitat table")

  } else {
    hab_list <- NULL
  }

  # For linksets ##########################################
  if(!is.null(existing_link)){

    linkset_df <- data.frame(name = NA, type = NA,
                             topology = NA, dist_max = NA,
                             inter = NA, ext_cost = NA)
    linkset_df <- linkset_df[-1, ]

    linkset_cost <- list()

    # where it starts and ends in the .xml file
    type_begin <- "<linksets>"
    type_end <- "</linksets>"

    # Get object type lines
    first_line <- min(which(file_data[, 1] == type_begin))
    last_line <- min(which(file_data[, 1] == type_end))

    for(i in 1:length(existing_link)){

      # Find the line at which this linkset starts
      target <- paste0("<string>", existing_link[i], "</string>")
      line_target <- which(file_data[first_line:last_line,
                                     1] == target)

      # Define the limits of this linkset block
      first_link_line <- first_line + line_target - 1
      last_link_line <- min(which(file_data[first_link_line:nrow(file_data),
                                            1] == "</entry>")) +
        first_link_line - 1

      # Check whether it is a Euclidean or cost-distance linkset
      if(any(file_data[first_link_line:last_link_line,
                       1] == "<typeDist>EUCLID</typeDist>")){
        type_dist <- "Euclidean"
      } else {
        type_dist <- "Cost"
      }

      # Check the topology
      if(any(file_data[first_link_line:last_link_line,
                       1] == "<topology>COMPLETE</topology>")){
        topology <- "Complete"
      } else {
        topology <- "Planar"
      }

      # Check for dist max
      line_distmax <- which(stringr::str_sub(file_data[first_link_line:last_link_line,
                                                       1], 1, 9) == "<distMax>")
      # get the dist max
      dist_max <- as.numeric(stringr::str_sub(file_data[first_link_line +
                                                          line_distmax -1,
                                                        1],
                                              10, -11))

      # Check for interhabitat link
      if(any(file_data[first_link_line:last_link_line,
                       1] == "<inter>false</inter>")){
        inter <- FALSE
      } else {
        inter <- TRUE
      }


      ## If COST, check for cost values or external raster file

      if(type_dist == "Cost"){

        # Check whether it is an external cost file or a cost table
        if(any(stringr::str_sub(file_data[first_link_line:last_link_line,
                                          1], 1, 13) == "<extCostFile>")){

          # Get the line with the path to the file
          line_extcost <- which(stringr::str_sub(file_data[first_link_line:last_link_line,
                                                           1],
                                                 1, 13) == "<extCostFile>") +
            first_link_line - 1

          # Get the external cost file
          ext_cost <- basename(stringr::str_sub(file_data[line_extcost, 1],
                                                14, -15))
          df_cost <- ext_cost

        } else {

          # Cost from raster codes
          ext_cost <- NA

          ### Get the cost values
          # Find the lines where they are defined (without the <costs> lines)
          first_cost_line <- which(file_data[first_link_line:last_link_line,
                                             1] == "<costs>") +
            first_link_line # first cost value for 0

          last_cost_line <- which(file_data[first_link_line:last_link_line,
                                            1] == "</costs>") +
            first_link_line - 2

          cost_lines <- file_data[first_cost_line:last_cost_line, 1]

          cost_values <- unlist(
            lapply(
              cost_lines,
              FUN = function(x){
                as.numeric(stringr::str_sub(x, 9, -10))
              }))

          # Create a data.frame to export
          df_cost_raw <- data.frame(code = 0:max(raster_codes),
                                    cost = cost_values)
          df_cost <- df_cost_raw[which(df_cost_raw$code %in% raster_codes), ]

        }

      } else {

        # if Euclidean, ext_cost is NA
        ext_cost <- NA
        df_cost <- "Euclidean distance - no cost parameter"
      }

      # Fill linkset_df
      linkset_df <- rbind(linkset_df,
                          data.frame(name = existing_link[i],
                                     type = type_dist,
                                     topology = topology,
                                     dist_max = dist_max,
                                     inter = inter,
                                     ext_cost = ext_cost))

      # Store the cost parameter
      linkset_cost[[i]] <- df_cost

    }

    # Name the items of the result list
    names(linkset_cost) <- existing_link
    linkset_list <- list(existing_link, linkset_df, linkset_cost)
    names(linkset_list) <- c("Linkset names", "Linkset description",
                             "Linkset cost parameters")

  } else {
    linkset_list <- NULL
  }

  # For graphs ##########################################
  if(!is.null(existing_graph)){
    graph_list <- existing_graph
  } else {
    graph_list <- NULL
  }

  # For metrics ##########################################
  if(!is.null(existing_metric)){

    df_graph_metric <- data.frame(graph = NA,
                                  metric = NA)
    df_graph_metric <- df_graph_metric[-1, ]

    for(i in 1:length(existing_graph)){

      graph_i <- existing_graph[i]

      id_metric <- which(stringr::str_sub(existing_metric,
                                          -nchar(graph_i),
                                          -1) == graph_i)

      if(length(id_metric) > 0){
        df_graph_metric <- rbind(df_graph_metric,
                                 data.frame(graph = rep(graph_i, length(id_metric)),
                                            metric = existing_metric[id_metric]))
      }
    }
    metric_list <- df_graph_metric
  } else {
    metric_list <- NULL
  }


  # Gather all results
  res <- list(paste0("Unique codes of the source raster: ", # raster codes
                     paste(raster_codes, collapse = ", ")),
              hab_list,
              linkset_list,
              graph_list,
              metric_list)

  names(res) <- c("Source raster", "Habitats",
                  "Linksets", "Graphs", "Metrics")

  return(res)

}



