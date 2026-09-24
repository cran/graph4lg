#' Show the objects created in the Graphab project
#'
#' @description The function lists the objects created in the project.
#' @param proj_path The path to the project .xml file.
#' @return The list of objects
#' @export
#' @author P. Savary
#' @keywords internal
#' @details Given that this function calls Java and is heavier than other
#' check-ups, functions reading the .xml file are usually preferred.
#' @return A list of the existing objects, with NULL values if none exists.
#' @examples
#' \dontrun{
#' proj_end_path <- "./graphab_example/graphab_example.xml"
#' graphab_show(proj_path = proj_end_path)
#' }


graphab_show <- function(proj_path){

  #########################################
  # Get java path
  java.path <- Sys.which("java")
  #########################################
  # Get graphab path
  version <- "graphab-3.0.jar"
  path_to_graphab <- paste0(rappdirs::user_data_dir(), "/graph4lg2_jar/", version)

  #######################################
  # Add '' to proj_path for cases with spaces in paths not transformed before
  if(all(stringr::str_sub(proj_path, 1, 1) != "'",
         stringr::str_sub(proj_path, 1, 1) != "'",
         stringr::str_detect(string = proj_path,
                             pattern = " "))){
    proj_path <- paste0("'", proj_path, "'")
  }


  #########################################
  # Command line
  cmd <- c("-Djava.awt.headless=true", "-jar", path_to_graphab,
           "--project", proj_path, "--show")

  #########################################
  # Run the command line
  rs <- system2(java.path, args = cmd, stdout = TRUE)

  hab_begin <- which(stringr::str_detect(string = rs,
                                         pattern = "= Habitats ="))

  link_begin <- which(stringr::str_detect(string = rs,
                                          pattern = "= Jeux de liens ="))

  ## Detect if system is French or English
  if(length(link_begin) == 0){

    link_begin <- which(stringr::str_detect(string = rs,
                                            pattern = "= Linksets ="))
    graph_begin <- which(stringr::str_detect(string = rs,
                                             pattern = "= Graphs ="))
    metric_begin <- which(stringr::str_detect(string = rs,
                                              pattern = "Metrics =="))

  } else {

    graph_begin <- which(stringr::str_detect(string = rs,
                                             pattern = "= Graphes ="))
    metric_begin <- which(stringr::str_detect(string = rs,
                                              pattern = "triques =="))
  }

  # Habitats
  if(link_begin - hab_begin == 2){
    exist_hab <- NULL
  } else {
    exist_hab <- rs[(hab_begin+1):(link_begin-2)]
  }
  # Linksets
  if(graph_begin - link_begin == 2){
    exist_link <- NULL
  } else {
    exist_link <- rs[(link_begin+1):(graph_begin-2)]
  }
  # Graphs
  if(metric_begin - graph_begin == 2){
    exist_graph <- NULL
  } else {
    exist_graph <- rs[(graph_begin+1):(metric_begin-2)]
  }
  # Metrics
  if(metric_begin == length(rs)){
    exist_metric <- NULL
  } else {
    exist_metric <- rs[(metric_begin+1):length(rs)]
  }

  # List results
  list_objects <- list(exist_hab, exist_link, exist_graph, exist_metric)
  names(list_objects) <- c("Habitats", "Linksets", "Graphs", "Metrics")

  # Return list
  return(list_objects)

}
