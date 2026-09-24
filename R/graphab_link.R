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
#' @param habitat (optional, default=NULL) A character string indicating the
#' name of the habitat types whose patches will be connected by the least-cost
#' paths. If \code{name='all'}, all existing habitats are considered.
#' If a single habitat is indicated (e.g., `habitat = 'forest'`), then
#' only this one is used. If several habitats are indicated, they are merged
#' before the computation (e.g., `habitat = c('forest', 'grass')`). If NULL,
#' each is separately considered.
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
#' @param inter (optional, default = FALSE) A logical indicating whether only
#' links among different habitat types are computed (if TRUE, no "intra-type"
#' links).
#' @param maxcost (optional, default = NULL) A value specifying the maximum
#' accumulated cost of a least-cost path. The computation stops beyond.
#' @param remcrosspath (optional, default = FALSE) A logical indicating whether
#' links crossing patches are removed (TRUE).
#' @param parallel.java An integer indicating how many computer cores are used
#' to run the .jar file. By default, \code{parallel.java = NULL}, and java sets
#' it according to local settings.
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
#' problem when computing BC index. See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' @export
#' @author P. Savary, T. Rudolph
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' df_cost <- data.frame(code = 1:5,
#'                       cost = c(1, 10, 100, 1000, 1))
#' graphab_link(proj_name = "graphab_example",
#'             habitat = "forest",
#'             distance = "cost",
#'             name = "lcp",
#'             cost = df_cost,
#'             topo = "complete")
#' }


graphab_link <- function(proj_name,         # character
                         distance = "cost", # cost or euclid
                         name, # character
                         habitat = NULL, #
                         cost = NULL, # NULL, data.frame code cost or ext file
                         topo = "planar", # planar or complete
                         inter = FALSE,
                         maxcost = NULL,
                         remcrosspath = FALSE,
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

  #########################################
  # Check for distance
  if(!inherits(distance, "character")){
    stop("'distance' must be a character string")
  } else if (!(distance %in% c("cost", "euclid"))){
    stop("'distance' must be equal to 'cost' or 'euclid'")
  }

  #########################################
  # Check for topo
  if(!inherits(topo, "character")){
    stop("'topo' must be a character string")
  } else if (!(topo %in% c("complete", "planar", "planarcost"))){
    stop("'topo' must be equal to 'complete', 'planar' or 'planarcost'.")
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
          } else {
            # Add '' to extcost for cases with spaces in paths
            if(all(stringr::str_sub(extcost, 1, 1) != "'",
                   stringr::str_sub(extcost, 1, 1) != "'",
                   stringr::str_detect(string = extcost,
                                       pattern = " "))){
              extcost <- paste0("'", extcost, "'")
            }
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
  # Check for inter
  if(!inherits(inter, "logical")){
    stop("'inter' must be a logical.")
  }

  #########################################
  # Check for maxcost
  if(!is.null(maxcost)){
    if(!inherits(maxcost, c("integer", "numeric"))){
      stop("'maxcost' must be a numeric value.")
    }
  }

  #########################################
  # Check for name
  if(!inherits(name, "character")){
    stop("'name' must be a character string")
  }

  #########################################
  # Check for habitat
  if(!is.null(habitat)){
    if(!inherits(habitat, "character")){
      stop("'habitat' must be a character string")
    } else if(any(habitat != "all")){
      if(!all(unlist(lapply(habitat,
                            FUN = function(x){
                              check_graphab_object(proj_path = proj_end_path,
                                                   object_type = "habitat",
                                                   name = x)
                            })))){
        stop("'habitat' must only include habitats that exist in the project.")
      }

      ## Check that habitats are not empty
      for(h in 1:length(habitat)){
        if(!file.exists(paste0(proj_path, "/",
                               proj_name, "/",
                               habitat[h], "/patches.gpkg"))){
          stop(paste0("The habitat type '", habitat[h], "' does not include any patch. ",
                         "Please check the raster and the patch constraints."))
        }

      }

    }
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
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg2_jar/", version)

  #########################################
  # Command line

  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

  if(!is.null(parallel.java)){
    cmd <- c(cmd, "-proc ", as.character(parallel.java))
  }

  cmd <- c(cmd,
           "--project", proj_end_path_cmd)

  # Add habitat specification
  if(!is.null(habitat)){
    if(length(habitat) > 1) {
      cmd <- c(cmd, "--mergehabitat", paste(habitat, collapse = ","))
    } else if(habitat == "all"){
      cmd <- c(cmd, "--mergehabitat", "ALL")
    } else {
      cmd <- c(cmd, "--usehabitat", habitat)
    }
  }

  cmd <- c(cmd, "--linkset",
           paste0("distance=", distance),
           paste0("name=", name),
           paste0("topo=", topo))

  if(inter){
    cmd <- c(cmd, "inter")
  }

  if(!is.null(maxcost)){
    cmd <- c(cmd, paste0("maxcost=", maxcost))
  }

  if(remcrosspath){
    cmd <- c(cmd, "remcrosspath")
  }

  if (distance == "cost") {
    if(inherits(cost, "data.frame")) {
      cmd <- c(cmd, vec_cost)
    } else {
      cmd <- c(cmd, paste0("extcost=", extcost))
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

  ## Check whether an error occurred
  if(length(rs) == 1){
    if(rs == 1){
      message("An error occurred")
    }
  }

  ## Check whether the linkset exists
  if(check_graphab_object(proj_path = proj_end_path,
                          object_type = "linkset",
                          name = name)){
    message(paste0("Link set '", name, "' has been created in the project '",
                   proj_name, "'."))
  } else {
    message("An error occurred")
  }

}
