#' Computes custom capacities of patches in the Graphab project
#'
#' @description The function computes custom capacities of patches
#' in the Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param mode A character string indicating the way capacities are
#' computed. It must be either:\itemize{
#' \item{\code{mode='area'}(default): The capacity of the patches is computed
#' as the area of each habitat patch. The argument \code{exp} makes it
#' possible to raise area to a power given by an exposant.}
#' \item{\code{mode='ext_file'}: The capacity of the patches is given by an
#' external .csv file. See argument \code{ext_file} below.}
#' \item{\code{mode='neigh'}: The capacity is computed depending on the
#' neighbouring raster cells from each habitat patch. The number of cells
#' with a value given by \code{codes} argument is summed up to the
#' distance \code{thr}. This number can be weighted according to the
#' \code{weight} argument.}
#' }
#' @param patch_codes (optional, default=NULL) An integer value or vector
#' specifying the codes corresponding to the habitat pixel whose corresponding
#' patches are included to compute the capacity as the area of the habitat
#' when \code{mode='area'}. Patches corresponding to other initial habitat
#' codes are weighted by 0.
#' @param exp An integer value specifying the power to which patch area are
#' raised when \code{mode='area'}. When not specified, \code{exp=1} by default.
#' @param ext_file A character string specifying the name of the .csv file in
#' which patch capacities are stored. It must be located either in the working
#' directory or in the directory defined by \code{proj_path}. It must have
#' as many rows as there are patches in the project and its column names
#' must include 'Id' and 'Capacity'. The 'Id' column must correspond to the
#' patch ID in the 'patches' layer (see \code{\link{get_graphab_metric}}).
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
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process. Increasing this
#' value can speed up the computations. Too large values may not be compatible
#' with your machine settings.
#' @details See more information in Graphab 2.8 manual:
#' \url{https://sourcesup.renater.fr/www/graphab/download/manual-2.8-en.pdf}
#' Be careful, when capacity has been changed. The last changes are taken into
#' account for subsequent calculations in a project.
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_capacity(proj_name = "grphb_ex",
#'                  mode = "area")
#' }





graphab_capacity <- function(proj_name,         # character
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


  # Distinguish the modes
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
    }


    # Get graphab path
    version <- "graphab-2.8.jar"
    path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)
    #### Command line
    cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
             "--project", proj_end_path,
             "--capa",  "area")


    if(!is.null(exp)){
      if(!inherits(exp, c("numeric", "integer"))){
        stop("'exp' argument must be a numeric or integer value")
      } else if (length(exp) > 1){
        stop("'exp' argument must be a numeric or integer value")
      } else {
        cmd <- c(cmd, paste0("exp=", exp))
      }
    }


    if(!is.null(patch_codes)){
      if(!inherits(patch_codes, c("numeric", "integer"))){
        stop("'patch_codes' argument must be a numeric or integer vector/value")
      } else {


        # Get habitat codes
        hab_codes <- get_graphab_raster_codes(proj_name = proj_name,
                                              mode = "habitat",
                                              proj_path = proj_path)

        # Check whether all habitat codes are in patch_codes
        if(!all(hab_codes %in% patch_codes)){
          # Distinguish the habitat codes
          hab_w <- patch_codes
          hab_now <- hab_codes[!(hab_codes %in% patch_codes)]

          cmd <- c(cmd,
                   paste0(paste(hab_w, collapse = ","), "=1"),
                   paste0(paste(hab_now, collapse = ","), "=0"))
        } else {
          cmd <- c(cmd, paste0(paste(patch_codes, collapse = ","), "=1"))
        }

      }
    }


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
    } else if(!is.null(patch_codes)){ # NEW
      message("Argument 'patch_codes' is not used when 'mode='ext_file''.")
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

      patches <- foreign::read.dbf(file = paste0(proj_path, "/",
                                                 proj_name, "/patches.dbf"))
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

    }

    ###############################################################
    # Get graphab path
    version <- "graphab-2.8.jar"
    path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)
    #### Command line
    cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
             "--project", proj_end_path,
             "--capa",  paste0("file=", ext_file),
             "id=Id", "capa=Capacity")



  } else if(mode == "neigh"){


    # Check for not null parameters and return a message if not used
    if(!is.null(ext_file)){
      message("Argument 'ext_file' is not used when 'mode='neigh''.")
    } else if(!is.null(exp)){
      message("Argument 'exp' is not used when 'mode='neigh''.")
    } else if(!is.null(patch_codes)){
      message("Argument 'patch_codes' is not used when 'mode='neigh''.")
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
    if(is.null(codes)){
      stop("'codes' must be integer values when 'mode='neigh''.")
    } else if (!(inherits(codes, c("integer", "numeric")))){
      stop("'codes' must be numeric or integer values.")
    } else {
      list_codes <- graph4lg::get_graphab_raster_codes(proj_name = proj_name,
                                                       mode = "all",
                                                       proj_path = proj_path)
      if(!(all(codes %in% list_codes))){
        stop("All 'codes' values must be values existing in the source raster.")
      }
    }


    #######################
    # Check for linkset
    if(is.null(linkset)){
      stop("'linkset' must be a character string when 'mode='neigh''.")
    } else if(!inherits(linkset, "character")){
      stop("'linkset' must be a character string when 'mode='neigh''.")
    } else if (!(paste0(linkset, "-links.csv") %in% list.files(path = paste0(proj_path, "/", proj_name)))){
      stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
    }

    #########################################
    # Check for cost_conv
    if(!is.logical(cost_conv)){
      stop("'cost_conv' must be a logical (TRUE or FALSE).")
    }

    #########################################
    # Check for weight
    if(!is.logical(weight)){
      stop("'weight' must be a logical (TRUE or FALSE).")
    }


    ###### Print used costs and codes
    df_cost <- graph4lg::get_graphab_linkset_cost(proj_name = proj_name,
                                                  linkset = linkset,
                                                  proj_path = proj_path)
    print(paste0("The following cost parameters will be used for ",
                 "computing capacities neighbouring the patches"))
    print(df_cost)


    ###############################################################
    # Get graphab path
    version <- "graphab-2.8.jar"
    path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/", version)
    #### Command line
    cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
             "--project", proj_end_path,
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
      message(paste0("Patch capacities have been updated. ",
                     "Use 'get_graphab_metric()' to get values"))
    }
  } else {
    message(paste0("Patch capacities have been updated. ",
                   "Use 'get_graphab_metric()' to get values"))
  }

}
