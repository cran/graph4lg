#' Get metrics computed at the node in the Graphab project
#'
#' @description The function gets the metrics computed at the node-level in
#' the Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @return A data.frame with metrics computed at the patch level.
#' @details The imported metrics describe the patches and have been computed
#' from the different graphs created in the Graphab project.
#' See more information in Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' get_graphab_metric(proj_name = "grphb_ex")
#' }

get_graphab_metric <- function(proj_name, # character
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

  #########################################
  # Check for patches.csv
  if(!file.exists(paste0(proj_path, "/", proj_name, "/patches.csv"))){
    df <- foreign::read.dbf(file = paste0(proj_path, "/",
                                          proj_name, "/patches.dbf"))
  } else {
    df <- utils::read.csv(file = paste0(proj_path, "/",
                                        proj_name, "/patches.csv"))
  }

  return(df)
}
