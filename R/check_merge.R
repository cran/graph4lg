#' Check whether the option 'nomerge' was used when building the landscape
#' graph with Graphab
#'
#' @description The function checks whether the option 'nomerge' was used when
#' building the landscape graph with Graphab
#' @param proj_end_path The path to the project .xml file.
#' @return The function returns a logical indicating whether 'nomerge' was used.
#' If nomerge=TRUE, then it returns FALSE. If nomerge=FALSE, it returns TRUE.
#' @export
#' @author P. Savary
#' @keywords internal
#' @examples
#' \dontrun{
#' proj_name <- "grphb_ex"
#' check_merge(proj_name = proj_name)
#' }


check_merge <- function(proj_end_path){

  #########################################################
  # Copy the .xml file as a .txt file in temp files to open it

  xml <- tempfile(pattern = ".txt")
  file.copy(from = proj_end_path,
            to = xml)
  file_data <- utils::read.table(xml)

  #####################################
  # Get merge line
  merge_line <- file_data[which(stringr::str_sub(file_data[, 1],
                                              1, 7) == "<merge>"), 1]
  merge_param <- stringr::str_sub(merge_line, 8, 11)

  if(merge_param == "fals"){
    merge_res <- FALSE
  } else if (merge_param == "true"){
    merge_res <- TRUE
  }
  return(merge_res)
}
