#' Get unique raster codes from a Graphab project
#'
#' @description The function extracts unique raster codes from a Graphab project
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml will be created.
#' @param mode A character string equal to either 'all' (default) or 'habitat'
#' indicating whether the returned codes are all the codes of the source raster
#' used for creating the project or only the code corresponding to the
#' habitat patches.
#' @param proj_path (optional) A character string indicating the path to the
#' directory that contains the project directory. It should be used when the
#' project directory is not in the current working directory. Default is NULL.
#' When 'proj_path = NULL', the project directory is equal to \code{getwd()}.
#' @return The function returns a vector of integer values corresponding to
#' the source raster codes (all the codes or only the one corresponding to
#' habitat patches).
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' proj_name <- "grphb_ex"
#' get_graphab_raster_codes(proj_name = proj_name,
#'                mode = "all")
#' }

get_graphab_raster_codes <- function(proj_name,
                                     mode = "all",
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


  #######################
  # Check for mode
  if(!inherits(mode, "character")){
    stop("'mode' must be a character string.")
  } else if (!(mode %in% c("all", "habitat"))){
    stop("'mode' must be equal to 'all' or 'habitat'.")
  }

  #########################################################
  # Copy the .xml file as a .txt file in temp files to open it

  xml <- tempfile(pattern = ".txt")
  file.copy(from = proj_end_path,
            to = xml)
  file_data <- utils::read.table(xml)


  #####################################
  # Get no data value
  na_code <- file_data[which(stringr::str_sub(file_data[, 1],
                                              1, 8) == "<noData>"), 1]
  na_code <- stringr::str_sub(na_code, 9, -10)

  if(na_code == "NaN"){
    na_pres <- FALSE
  } else {
    na_pres <- TRUE

    # Remove the .0 value
    if(stringr::str_sub(na_code, -2, -1) == ".0"){
      na_code <- stringr::str_sub(na_code, 1, -3)
    }

  }

  ###############################################
  # Get the code values

  if(mode == "all"){

    first_code <- min(which(file_data[, 1] == "<codes>")) + 1
    last_code <- min(which(file_data[, 1] == "</codes>")) - 1
    vec_codes <- file_data[first_code:last_code, 1]

  } else if(mode == "habitat"){

    first_code <- min(which(file_data[, 1] == "<patchCodes>")) + 1
    last_code <- min(which(file_data[, 1] == "</patchCodes>")) - 1
    vec_codes <- file_data[first_code:last_code, 1]

  }

  # Extract the codes
  vec_codes <- unlist(lapply(vec_codes,
                      FUN = function(x){stringr::str_sub(x, 6, -7)}))

  # Remove No data if present
  if(na_pres){
    if(na_code %in% vec_codes){
    vec_codes <- vec_codes[-which(vec_codes == na_code)]
    # Print a message
    message(paste0("No data value is equal to ", na_code))
    }
  }

  vec_codes <- as.numeric(vec_codes)

  return(vec_codes)


}
