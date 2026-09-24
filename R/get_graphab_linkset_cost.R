#' Get cost values associated with a linkset in a Graphab project
#'
#' @description The function extracts the cost parameters associated with a
#' linkset in a Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml has been created.
#' @param linkset A character string indicating the name of the link set used
#' to create the graph. Link sets can be created with \code{\link{graphab_link}}.
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
#' proj_name <- "graphab_example"
#' get_graphab_linkset_cost(proj_name = proj_name,
#'                linkset = "forest_link_planar")
#' }


get_graphab_linkset_cost <- function(proj_name,
                                     linkset,
                                     proj_path = NULL){

  #########################################
  # Get the project information
  project_info <- graphab_project_desc(proj_name = proj_name,
                                       proj_path = proj_path)
  project_linkset <- project_info[["Linksets"]]

  #########################################
  # Check for linkset class
  if(!inherits(linkset, "character")){
    stop("'linkset' must be a character string")
  } else if (!(linkset %in% project_linkset[["Linkset names"]])){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  ###### Get the parameters
  linkset_param <- project_linkset[["Linkset cost parameters"]][[linkset]]

  return(linkset_param)

}
