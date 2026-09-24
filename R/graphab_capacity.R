#' Computes custom capacities of patches in the Graphab project
#'
#' @description The function computes custom capacities of patches
#' in the Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param habitat A character string indicating the name of the habitat type
#' whose patches will get a new capacity value. When the habitat type is defined
#' from the source raster codes, all \code{mode} options are available. When
#' it is of type 'vector', \code{mode} can only be `ext_file` or `neigh`.
#' Only one habitat type can be treated at a time.
#' @param mode A character string indicating the way capacities are
#' computed. It must be either:\itemize{
#' \item{\code{mode='area'}(default): The capacity of the patches is computed
#' as the area of each habitat patch. The argument \code{exp} makes it
#' possible to raise area to a power given by an exponent.}
#' \item{\code{mode='ext_file'}: The capacity of the patches is given by an
#' external .csv file. See argument \code{ext_file} below.}
#' \item{\code{mode='neigh'}: The capacity is computed as a function of the
#' neighboring raster cells of each habitat patch. The number of cells
#' with a value given by the \code{codes} argument is summed up to the
#' distance \code{thr}. This number can be weighted according to the
#' \code{weight} argument.}
#' }
#' @param patch_codes (default=NULL) Deprecated parameter from verson 1.8.
#' @param exp An integer value specifying the power to which patch area are
#' raised when \code{mode='area'}. When not specified, \code{exp=1} by default.
#' @param ext_file A character string specifying the name of the .csv file in
#' which patch capacities are stored. It must be located either in the working
#' directory or in the directory defined by \code{proj_path}. It must have
#' as many rows as there are patches of the considered habitat type. Its column
#' names must include 'Id' and 'Capacity'. The 'Id' column must correspond to
#' the patches' ID in the 'patches' layer (see \code{\link{get_graphab_metric}}).
#' The 'Capacity' column must contain the corresponding patch capacities to
#' assign each patch.
#' @param thr (optional, default=NULL) An integer or numeric value indicating
#' the maximum distance in cost distance units (except when
#' \code{cost_conv = TRUE}) at which cells are considered for computing the
#' capacity when \code{mode='neigh'}.
#' @param linkset (optional, default=NULL) A character string indicating the
#' name of the link set used to take distance into account when computing
#' the capacity. Only used when \code{mode='neigh'}. Link sets can be
#' created with \code{\link{graphab_link}}.
#' @param codes An integer value or a vector of integer values specifying the
#' codes of the raster cells taken into account when computing the capacity in
#' the neighbourhood of the patches, when \code{mode='neigh'}.
#' @param cost_conv FALSE (default) or TRUE. Logical indicating whether numeric
#' \code{thr} values are converted from cost-distance into Euclidean distance
#' using a log-log linear regression. See also \code{\link{convert_cd}}
#' function. Only used when \code{mode='neigh'}.
#' @param weight A logical indicating whether the cells are weighted by a
#' weight decreasing with the distance from the patches (TRUE) or not (FALSE).
#' The weights follow a negative exponential decline such that
#' wi = exp(-alpha*di), where wi is the weight of cell i, di its distance from
#' the patch and alpha a parameter determined such that wi = 0.05 when di = thr.
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
#' @details See more information in Graphab 3.0 manual:
#' \url{https://thema.umlp.fr/productions/software/graphab/download/manual-3.0-en.pdf}
#' Be careful, when capacity has been changed. The last changes are taken into
#' account for subsequent calculations in a project.
#' @export
#' @author P. Savary
#' @references \insertRef{foltete2012software}{graph4lg}
#' \insertRef{foltete2021graphab}{graph4lg}
#' \insertRef{savary2024multiple}{graph4lg}
#' @examples
#' \dontrun{
#' graphab_capacity(proj_name = "graphab_example",
#'                  mode = "area")
#' }

graphab_capacity <- function(proj_name,         # character
                             habitat, # character
                             mode = "area", # character
                             patch_codes = NULL, # NULL or integer vector
                             exp = NULL, # integer
                             ext_file = NULL, # character
                             thr = NULL, # threshold NULL or numerical vector
                             linkset = NULL, # cost or euclid
                             codes = NULL, # NULL or integer vector
                             cost_conv = FALSE, # FALSE (default) or TRUE
                             weight = FALSE, # default FALSE, but TRUE for link weighting
                             proj_path = NULL, # if NULL getwd() otherwise a character path
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

  ## Get the project information
  project_info <- graphab_project_desc(proj_name = proj_name,
                                       proj_path = proj_path)
  project_habitat <- project_info[["Habitats"]]
  project_linkset <- project_info[["Linksets"]]
  project_raster <- project_info[["Source raster"]]

  #########################################
  # Check for parallel.java
  if(!is.null(parallel.java)){
    if(!inherits(parallel.java, c("numeric", "integer"))){
      stop("'parallel.java' must be a numeric or integer value.")
    }
  }

  # Get graphab path
  version <- "graphab-3.0.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg2_jar/", version)

  #########################################
  # Check for habitat class and mode compatibility
  if(!inherits(habitat, "character")){
    stop("'habitat' must be a character string")
  } else if(!(habitat %in% project_habitat[["Habitat names"]])){
    stop("The habitat type you refer to does not exist.
         Please use graphab_habitat() before.")
  } else {

    # Get habitat info and type
    habitat_df <- project_habitat[["Habitat table"]]
    habitat_df <- habitat_df[which(habitat_df$name == habitat), ]

    # If vector habitat, mode area is not possible
    if(all(c(habitat_df$type == "Vector", mode == "area"))){
      stop("mode='area' cannot be used with vector habitat types.")
    }

  }

  #########################################
  # Check for patch_codes class
  if(!is.null(patch_codes)){
    stop(paste0("Argument 'patch_codes' is deprecated and not used anymore ",
                "in graph4lg >= 2.0. Please use graph4lg <= 1.8 if you ",
                "really need to use an older version of graphab_capacity(). ",
                "Please note that most functionalities have been conserved, ",
                "yet with a new syntax."))
  }

  #########################################
  # Check for weight and cost_conv
  if(!inherits(weight, "logical")){
    stop("'weight' must be a logical.")
  } else if(!inherits(cost_conv, "logical")){
    stop("'cost_conv' must be a logical.")
  }


  ######## Commands distinguished by modes
  ########################################################
  if(mode == "area"){

    # Check for not null parameters and return a message if not used
    if(!is.null(linkset)){
      message("Argument 'linkset' is not used when 'mode='area''.")
    } else if(!is.null(codes)){
      message("Argument 'codes' is not used when 'mode='area''.")
    } else if(!is.null(thr)){
      message("Argument 'thr' is not used when 'mode='area''.")
    } else if(!is.null(ext_file)){
      message("Argument 'ext_file' is not used when 'mode='area''.")
    } else if(weight){
      message("Argument 'weight' is not used when 'mode='area''.")
    } else if(cost_conv){
      message("Argument 'cost_conv' is not used when 'mode='area''.")
    }

    #### Command line

    cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

    if(!is.null(parallel.java)){
      cmd <- c(cmd, "-proc ", as.character(parallel.java))
    }

    cmd <- c(cmd,
             "--project", proj_end_path_cmd,
             "--usehabitat", habitat,
             "--capa",  "area")

    # Add the exponent associated with area when needed
    if(!is.null(exp)){
      if(!inherits(exp, c("numeric", "integer"))){
        stop("'exp' argument must be a numeric or integer value")
      } else if (length(exp) > 1){
        stop("'exp' argument must be a numeric or integer value")
      } else {
        cmd <- c(cmd, paste0("exp=", exp))
      }
    }

    # Add the habitat raster codes in the command with a weight of 1
    # (space in the df value for multiple codes, need to split and paste)
    cmd <- c(cmd,
             paste0(paste(unlist(stringr::str_split(habitat_df$rast_codes,
                                                    ", ")),
                          collapse = ","),
                    "=1"))

    ########################################################
  } else if (mode == "ext_file"){

    # Check for not null parameters and return a message if not used
    if(!is.null(linkset)){
      message("Argument 'linkset' is not used when 'mode='ext_file''.")
    } else if(!is.null(codes)){
      message("Argument 'codes' is not used when 'mode='ext_file''.")
    } else if(!is.null(thr)){
      message("Argument 'thr' is not used when 'mode='ext_file''.")
    } else if(!is.null(exp)){
      message("Argument 'exp' is not used when 'mode='ext_file''.")
    } else if(weight){
      message("Argument 'weight' is not used when 'mode='ext_file''.")
    } else if(cost_conv){
      message("Argument 'cost_conv' is not used when 'mode='ext_file''.")
    }

    ############
    # Check ext_file

    if(is.null(ext_file)){
      stop("'ext_file' argument must be specified when 'mode='ext_file''.")

    } else if(!inherits(ext_file, "character")){
      stop("'ext_file' argument must be a character string specifying the
           path to an existing .csv file.")
    } else if(!file.exists(normalizePath(ext_file, mustWork = FALSE))){
      stop(paste0(normalizePath(ext_file, mustWork = FALSE),
                  " is not an existing .csv file."))
    } else {

      # Open ext_file to check for column names
      capa_file <- utils::read.csv(file = ext_file)

      patches <- suppressWarnings(
        sf::st_drop_geometry(
          sf::read_sf(
            paste0(proj_path, "/",
                   proj_name, "/",
                   habitat, "/patches.gpkg"),
            query = "SELECT * FROM patches",
            as_tibble = FALSE)))
      nb_patches <- nrow(patches)

      # Check for column names
      if(!all(c("Id", "Capacity") %in% colnames(capa_file))){
        stop(paste0("Column names of ", ext_file, " must include ",
                    "'Id' and 'Capacity'."))
      } else if(nrow(capa_file) != nb_patches){
        stop(paste0(ext_file, " must include as many rows as there ",
                    "are patches in the project."))
      } else if(!all(patches$Id %in% capa_file$Id)){
        stop(paste0(ext_file, " must include all the Id of ",
                    "the patches in the project."))
      }

      #######################################
      # Add '' to ext_file for cases with spaces in paths
      if(all(stringr::str_sub(ext_file, 1, 1) != "'",
             stringr::str_sub(ext_file, 1, 1) != "'",
             stringr::str_detect(string = ext_file,
                                 pattern = " "))){
        ext_file_cmd <- paste0("'", ext_file, "'")
      } else {
        ext_file_cmd <- ext_file
      }

    }

    ################
    #### Command line

    cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

    if(!is.null(parallel.java)){
      cmd <- c(cmd, "-proc ", as.character(parallel.java))
    }

    cmd <- c(cmd,
             "--project", proj_end_path_cmd,
             "--usehabitat", habitat,
             "--capa",  paste0("file=", ext_file_cmd),
             "id=Id", "capa=Capacity")

    ########################################################
  } else if(mode == "neigh"){


    # Check for not null parameters and return a message if not used
    if(!is.null(ext_file)){
      message("Argument 'ext_file' is not used when 'mode='neigh''.")
    } else if(!is.null(exp)){
      message("Argument 'exp' is not used when 'mode='neigh''.")
    }

    #######################
    # Check for thr
    if(is.null(thr)){
      stop("'thr' must be a specified numeric or integer value
           when 'mode='neigh''.")
    } else if (!(inherits(thr, c("integer", "numeric")))){
      stop("'thr' must be a specified numeric or integer value
           when 'mode='neigh''.")
    }

    #######################
    # Check for codes

    ### Source raster codes
    raster_codes <- get_graphab_raster_codes(proj_name = proj_name,
                                             mode = 'all',
                                             proj_path = proj_path)

    if(is.null(codes)){
      stop("'codes' must be integer values when 'mode='neigh''.")
    } else if (!(inherits(codes, c("integer", "numeric")))){
      stop("'codes' must be numeric or integer values.")
    } else if(!(all(codes %in% raster_codes))){
      stop("All 'codes' values must be values existing in the source raster.")
    }

    #######################
    # Check for linkset
    if(is.null(linkset)){
      stop("'linkset' must be a character string when 'mode='neigh''.")
    } else if(!inherits(linkset, "character")){
      stop("'linkset' must be a character string when 'mode='neigh''.")
    } else if (!(linkset %in% project_linkset[["Linkset names"]])){
      stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
    }

    ###### Print used costs and codes
    df_cost <- project_linkset[["Linkset cost parameters"]][[linkset]]
    print(paste0("The following cost parameters will be used to ",
                 "weight the distances to neighbouring patches when ",
                 "computing the new capacities."))
    print(df_cost)

    ###################
    #### Command line
    cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab)

    if(!is.null(parallel.java)){
      cmd <- c(cmd, "-proc ", as.character(parallel.java))
    }

    cmd <- c(cmd,
             "--project", proj_end_path_cmd,
             "--usehabitat", habitat,
             "--uselinkset", linkset,
             "--capa")

    #############
    # Add max_cost converting it or not
    if(cost_conv){
      cmd <- c(cmd, paste0("maxcost={", thr, "}"))
    } else {
      cmd <- c(cmd, paste0("maxcost=", thr))
    }

    #############
    # Add codes
    vec_codes <- paste0("codes=", paste(codes, collapse = ",", sep = ""))
    cmd <- c(cmd, vec_codes)

    ###########
    # Add weight if necessary
    if(weight){
      cmd <- c(cmd, "weight")
      message(paste0("Weighting parameter: ",
                     "p(", thr, ") = 0.05"))
    }

    #########################################################################

  } else {
    stop("'mode' must be a character string equal to either 'area',
         'ext_file' or 'neigh'.")

  }

  ##############################################################################
  ##############################################################################
  #########################################
  # Check for Graphab
  gr <- get_graphab(res = FALSE, return = TRUE)

  if(gr == 1){
    message("Graphab has been downloaded")
  }

  #########################################
  # Get java path
  java.path <- Sys.which("java")


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
      message(paste0("Patch capacities of '", habitat,
                     "' habitat have been updated. ",
                     "Use 'get_graphab_metric()' to get values"))
    }
  } else {
    message(paste0("Patch capacities of '", habitat,
                   "' habitat have been updated. ",
                   "Use 'get_graphab_metric()' to get values"))
  }

}
