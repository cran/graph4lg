#' Get linkset computed in the Graphab project
#'
#' @description The function gets a linkset computed in the Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is.
#' @param linkset A character string indicating the name of the link set
#' whose properties are imported. The link set has been created with Graphab
#' or using \code{\link{graphab_link}} function.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @return A data.frame with the link properties (from, to, cost-distance,
#' Euclidean distance), but not their spatial information (no spatial lines).
#' @details See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}.
#' This function works if \code{link{get_graphab}} function works correctly.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' get_graphab_linkset(proj_name = "graphab_example",
#'                linkset = "forest_link_planar")
#' }


get_graphab_linkset <- function(proj_name,
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

  ## Check project version
  check_graphab_version(proj_path = proj_end_path)

  #########################################
  # Check for linkset class and the existence of the linkset
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (!check_graphab_object(proj_path = proj_end_path,
                                   object_type = "linkset",
                                   name = linkset)){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  ## Link all the project files
  proj_files <- list.files(path = paste0(proj_path, "/", proj_name),
                           recursive = TRUE, full.names = TRUE)

  ## Find the file which matches the linkset .pgkg file name
  linkset_file <- proj_files[which(!is.na(
    stringr::str_match(string = proj_files,
                       pattern = paste0("/", linkset, "-links.gpkg"))))]

  ## Load the data from the linkset file
  df <- suppressWarnings(
    sf::st_drop_geometry(
      sf::read_sf(linkset_file,
                  as_tibble = FALSE)
    )
  )

  return(df)
}
