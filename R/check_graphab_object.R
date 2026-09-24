#' Check whether a Graphab object exists in the Graphab project.
#'
#' @description The function checks whether a given object exists in the
#' Graphab project. This can be a habitat, a linkset or a graph.
#' @param proj_path The path to the project .xml file.
#' @param object_type A character string with the name of the object type whose
#' list will be checked. Can be either "habitat", "linkset", or "graph".
#' @param name A character string (or vector) with the name(s) of the
#' object(s) whose existence is checked.
#' @return A logical specifying whether the object exists.
#' @export
#' @author P. Savary
#' @keywords internal
#' @examples
#' \dontrun{
#' proj_end_path <- "./graphab_example/graphab_example.xml"
#' check_graphab_object(proj_path = proj_end_path,
#'                      object_type = "habitat",
#'                      name = "forest")
#' }

check_graphab_object <- function(proj_path,
                                 object_type,
                                 name){

  #########################################################
  # Copy the .xml file as a .txt file in temp files to open it
  xml <- tempfile(pattern = ".txt")
  file.copy(from = proj_path,
            to = xml)
  file_data <- utils::read.table(xml)

  if(object_type == "habitat"){
    type_begin <- "<habitats>"
    type_end <- "</habitats>"
    target <- paste0("<string>", name, "</string>")

    # Get object type lines
    if(length(which(file_data[, 1] == type_begin)) != 0){
      first_line <- min(which(file_data[, 1] == type_begin))
      last_line <- min(which(file_data[, 1] == type_end))
    } else {
      stop("No habitat has been created in this project.")
    }

  } else if(object_type == "linkset"){
    type_begin <- "<linksets>"
    type_end <- "</linksets>"
    target <- paste0("<string>", name, "</string>")

    # Get object type lines
    if(length(which(file_data[, 1] == type_begin)) != 0){
      first_line <- min(which(file_data[, 1] == type_begin))
      last_line <- min(which(file_data[, 1] == type_end))
    } else {
      stop("No link set has been created in this project.")
    }

  } else if(object_type == "graph"){
    type_begin <- "<graphs>"
    type_end <- "</graphs>"
    target <- paste0("<string>", name, "</string>")

    # Get object type lines - in two steps because of MultiGraph possibility

    if(length(which(file_data[, 1] == type_begin)) != 0){
      first_line <- min(which(file_data[, 1] == type_begin))
    } else {
      stop("No graph has been created in this project.")
    }

    last_lines <- which(file_data[, 1] == type_end)
    if(any(file_data[, 1] == "<MultiGraph>")){
      last_multigraph <- which(file_data[, 1] == "</MultiGraph>")
      last_lines <- last_lines[-which(last_lines == (last_multigraph - 1))]
      last_line <- min(last_lines)
    } else {
      last_line <- min(last_lines)
    }
  }

  # Search for the target name in these lines
  line_check <- any(file_data[first_line:last_line,
                              1] == target)

  # Return a logical
  return(line_check)
}
