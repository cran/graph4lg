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
#' Euclidean distance)
#' @details See more information in Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}.
#' This function works if \code{link{get_graphab}} function works correctly.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' get_graphab_linkset(proj_name = "grphb_ex",
#'                linkset = "lkst1")
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


  #########################################
  # Check for linkset class
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (length(list.files(path = paste0(proj_path, "/", proj_name),
                               pattern = "-links.csv")) == 0){
    stop("There is not any linkset in the project you refer to.
         Please use graphab_link() before.")
  } else if (!(paste0(linkset, "-links.csv") %in% list.files(path = paste0(proj_path,
                                                                           "/", proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  df <- foreign::read.dbf(file = paste0(proj_path, "/",
                                        proj_name, "/",
                                        linkset, "-links.dbf"))


  return(df)
}
