#' Create a link set in the Graphab project
#'
#' @description The function creates a link set between habitat patches in the
#' Graphab project.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param distance A character string indicating whether links between patches
#' are computed based on:\itemize{
#' \item{Shortest cost distances: \code{distance='cost'} (default)}
#' \item{Straight Euclidean distances: \code{distance='euclid'}}
#' }
#' In the resulting link set, each link will be associated with its
#' corresponding cost-distance and the length of the least-cost path in meters
#' (if \code{distance='cost'}) or with its length in Euclidean distance
#' (if \code{distance='euclid'})
#' @param name A character string indicating the name of the created linkset.
#' @param cost This argument could be:\itemize{
#' \item{A \code{data.frame} indicating the cost values associated to each
#' raster cell value. These values refer to the raster used to create the
#' project with \code{graphab_project}. The data.frame must have two
#' columns:\itemize{
#' \item{'code': raster cell values}
#' \item{'cost': corresponding cost values}
#' }}
#' \item{The path to an external raster file in .tif format with cost values.}
#' }
#' @param topo A character string indicating the topology of the created
#' link set. It can be:\itemize{
#' \item{Planar (\code{topo='planar'} (default)): a planar set of links is
#' created. It speeds up the computation but will prevent from creating
#' complete graphs with \code{\link{graphab_graph}}.}
#' \item{Complete (\code{topo='complete'}): a complete set of links is created.
#' A link is computed between every pair of patches.}
#' }
#' @param remcrosspath (optional, default = FALSE) A logical indicating whether
#' links crossing patches are removed (TRUE).
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @details By default, links crossing patches are not ignored nor broken into
#' two links. For example, a link from patches A to C crossing patch B
#' is created. It takes into account the distance inside patch B. It can be a
#' problem when computing BC index. See more information in Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
#' @export
#' @author P. Savary, T. Rudolph
#' @examples
#' \dontrun{
#' df_cost <- data.frame(code = 1:5,
#'                       cost = c(1, 10, 100, 1000, 1))
#' graphab_link(proj_name = "grphb_ex",
#'             distance = "cost",
#'             name = "lcp",
#'             cost = df_cost,
#'             topo = "complete")
#' }


graphab_link <- function(proj_name,         # character
                         distance = "cost", # cost or euclid
                         name, # character
                         cost = NULL, # NULL, data.frame code cost or ext file
                         topo = "planar", # planar or complete
                         remcrosspath = FALSE,
                         proj_path = NULL, # if null getwd() otherwise a character path
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

  proj_end_path <- paste0(proj_path, "/", proj_name, "/", proj_name, ".xml")

  #########################################
  # Check for distance
  if(!inherits(distance, "character")){
    stop("'distance' must be a character string")
  } else if (!(distance %in% c("cost", "euclid"))){
    stop("'distance' must be equal to 'cost' or 'euclid'")
  }


  #########################################
  # Check for remcrosspatch
  if(!inherits(remcrosspath, "logical")){
    stop("'remcrosspath' must be a logical.")
  }

  ###########################################################################################
  # Check cost argument

  if(distance == "cost"){

    # Scenario 1: 'cost' argument is a lookup table (data.frame)
    if(inherits(cost, "data.frame")){

      if(!all(c("code", "cost") %in% colnames(cost))){
        stop("The columns of cost must include 'code' and 'cost'")
      } else if (any(is.na(as.numeric(cost$code)))){
        stop("'code' column must include numeric values")

      } else if (any(is.na(as.numeric(cost$cost)))){
        stop("'cost' column must include numeric values")

      }

      if(inherits(cost$code, c("factor", "character"))){
        cost$code <- as.numeric(as.character(cost$code))
      }

      if(inherits(cost$cost, c("factor", "character"))){
        cost$cost <- as.numeric(as.character(cost$cost))
      }

      rast_codes <- graph4lg::get_graphab_raster_codes(proj_name = proj_name,
                                                       mode = 'all',
                                                       proj_path = proj_path)

      if(!all(rast_codes %in% cost$code)){
        stop("'code' column must include all the raster code values.")
      }

      # Cost values argument
      ncode <- nrow(cost)

      vec_cost <- c()
      for(i in 1:ncode){
        vec_cost <- c(vec_cost, paste0(cost[i, "code"], "=", cost[i, 'cost']))
      }


      # Scenario 2: 'cost' argument is a filename referencing a
      # cost surface raster (GeoTiff)
    } else {

      if(inherits(cost, "character")){

        if(stringr::str_sub(cost, start=-4L) == '.tif'){
          extcost <- cost
          if(!(file.exists(normalizePath(extcost, mustWork = FALSE)))){
            stop(paste0(extcost,
                        " must be an existing cost surface raster file ('.tif')"))
          }
        }

      } else {
        stop("'cost' must be a data.frame or a cost surface raster file ('.tif')")
      }

    }

    # If cost is not NULL with euclid option
  } else if (!is.null(cost)) {

    message("'cost' argument is ignored with 'distance = euclid'")

  }

  #########################################
  # Check for name
  if(!inherits(name, "character")){
    stop("'name' must be a character string")
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
  version <- "graphab-2.8.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)

  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
           "--project", proj_end_path,
           "--linkset",
           paste0("distance=", distance),
           paste0("name=", name))

  if(topo == "complete"){
    cmd <- c(cmd, "complete")
  }

  if(remcrosspath){
    cmd <- c(cmd, "remcrosspath")
  }


  if (distance == "cost") {
    if(inherits(cost, "data.frame")) {
      cmd <- c(cmd, vec_cost)
    } else {
      cmd <- c(cmd, paste0('extcost=', extcost))
    }
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
      if(file.exists(paste0(proj_path, "/", proj_name, "/", name, "-links.shp"))){
        message(paste0("Link set '", name, "' has been created in the project ",
                       proj_name))
      } else {
        message("The link set creation did not succeed.")
      }
    }
  } else {

    if(file.exists(paste0(proj_path, "/", proj_name, "/", name, "-links.shp"))){
      message(paste0("Link set '", name, "' has been created in the project ",
                     proj_name))
    } else {
      message("The link set creation did not succeed.")
    }
  }

}
