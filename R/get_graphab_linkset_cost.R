#' Get cost values associated with a linkset in a Graphab project
#'
#' @description The function extracts the cost values associated with a
#' linkset in a Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml will be created.
#' @param linkset (optional, default=NULL) A character string indicating the
#' name of the link set used to create the graph. Link sets can be created
#' with \code{\link{graphab_link}}.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @return The function returns a data.frame with the cost values corresponding
#' to every raster code value.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' proj_name <- "grphb_ex"
#' get_graphab_linkset_cost(proj_name = proj_name,
#'                linkset = "lkst1")
#' }


get_graphab_linkset_cost <- function(proj_name,
                                     linkset,
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


  #########################################
  # Check for linkset class
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (!(paste0(linkset, "-links.csv") %in% list.files(path = paste0("./", proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  #########################################################
  # Copy the .xml file as a .txt file in temp files to open it

  xml <- tempfile(pattern = ".txt")
  file.copy(from = proj_end_path,
            to = xml)
  file_data <- utils::read.table(xml)

  #########################################################
  # Get the linkset parameters

  # Line numbers with linkset names
  lines_linkset_names <- which(file_data[, 1] == "<Linkset>") + 1
  # Linkset names
  names_linkset <- stringr::str_sub(file_data[lines_linkset_names, 1], 7, -8)
  # Line number of the linkset with same name as linkset argument
  line_linkset <- lines_linkset_names[which(names_linkset == linkset)]

  # Type of dist of the linkset
  type_dist <- stringr::str_sub(file_data[line_linkset + 2, 1], 13, -14)

  if(type_dist == "1"){

    message(paste0("Linkset ", linkset, " is a Euclidean linkset without ",
                   "associated cost values"))

    df_cost <- data.frame(code = NA,
                          cost = NA)


  } else if(type_dist == "2"){

    # Get the code values
    # codes <- get_graphab_raster_codes(proj_name = proj_name,
    #                                   mode = "all")
    codes <- graph4lg::get_graphab_raster_codes(proj_name = proj_name,
                                                mode = "all",
                                                proj_path = proj_path)

    # Get the cost values

    # Line numbers with linkset names
    lines_costs <- which(file_data[, 1] == "<costs>")
    lines_end_costs <- which(file_data[, 1] == "</costs>")

    # First line number with <costs> + 2 is the first cost value of the linkset
    # The closest after, + 1 to get cost values, but + 2 because first cost is 0
    first_cost_line <- min(lines_costs[lines_costs > line_linkset]) + 2

    # Last line number with cost values
    # first line with </costs> after first cost line and -1 to get the last cost
    last_cost_line <- min(lines_end_costs[lines_end_costs > first_cost_line]) - 1

    # Extract the costs
    cost_values <- file_data[first_cost_line:last_cost_line, 1]

    # Get the numeric values
    cost_values <- unlist(lapply(cost_values,
                                 FUN = function(x){stringr::str_sub(x, 9, -10)}))
    cost_values <- as.numeric(cost_values)

    if(length(codes) != length(cost_values)){
      # if(any(cost_values == 0)){
      #   cost_values <- cost_values[-which(cost_values == 0)]
      # }
      message("The number of cost values does not strictly correspond ",
              "to the number of code values. Cost values were probably ",
              "given for absent code values.")
      df_cost <- data.frame(code = NA,
                            cost = NA)

    } else {
      # Create a data.frame to export
      df_cost <- data.frame(code = codes,
                            cost = cost_values)
    }

  }

  return(df_cost)


}
