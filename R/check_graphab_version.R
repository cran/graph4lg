#' Check whether the Graphab project has been created with graphab.jar >= 3.0
#'
#' @description The function checks whether the Graphab project has been
#' created with graphab.jar >= 3.0
#' @param proj_path The path to the project .xml file.
#' @return If the version is < graphab-3.0.jar, the function returns a
#' message specifying that the project has been created with an earlier
#' version and that this version of graph4lg should not be used, or that
#' another project should be created.
#' @export
#' @author P. Savary
#' @keywords internal
#' @examples
#' \dontrun{
#' proj_end_path <- "./graphab_example/graphab_example.xml"
#' check_graphab_version(proj_path = proj_end_path)
#' }


check_graphab_version <- function(proj_path){

  #########################################################
  # Copy the .xml file as a .txt file in temp files to open it

  xml <- tempfile(pattern = ".txt")
  file.copy(from = proj_path,
            to = xml)
  file_data <- utils::read.table(xml)

  #####################################
  # Get merge line

  if(stringr::str_sub(file_data[3, 1], 1, 10) != "<version>3"){
    message("This project has been created with Graphab < 3.0.")
    message("graph4lg <= 1.9 should be used for this project.")
    message("You can also create a new project with this version of graph4lg.")
  }
}
