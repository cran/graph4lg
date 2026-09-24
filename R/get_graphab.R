#' Download Graphab if not present on the user's machine
#'
#' @description The function checks for the presence of Graphab (.jar) on the
#' user's machine and downloads it if absent. It also checks that users have
#' installed java on their machine.
#'
#' @param res Logical indicating whether a message says if Graphab has been
#' downloaded or not.
#' @param return Logical indicating whether the function returns a 1 or a 0
#' to indicate if Graphab has been downloaded or not.
#' @return If res = TRUE, the function displays a message indicating to users
#' what has been done.
#' If return = TRUE, it returns a 0 if Graphab is already on the machine and
#' a 1 if it has been downloaded.
#' @details If the download does not work, you can create a directory named
#' 'graph4lg_jar' in the directory \code{rappdirs::user_data_dir()} and copy
#' Graphab software downloaded from
#' \url{https://thema.umlp.fr/productions/download.php?name=graphab&version=3.0.0&username=Graph4lg&institution=R}
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' get_graphab()
#' }

get_graphab <- function(res = TRUE, return = FALSE){

  #########################################################################
  #########################################################################
  #########################################################################

  # Check if java is installed
  if(Sys.which("java") == ""){
    warning("Please install java if you want to use Graphab")
  }

  # Check for graphab.jar and download it if necessary

  data_dir <- rappdirs::user_data_dir()

  if("graphab-2.8.jar" %in% list.files(paste0(data_dir, "/graph4lg_jar"))) {

    message("Graphab 2.8 is on your machine, but Graphab 3.0 will be used.")
    message("To keep using Graphab 2.8, please use graph4lg < 2.0.")
    message("You can download it at:")
    message("remotes::install_gitlab('psavary3/graph4lg@archive-1-9')")

  }

  if(!("graphab-3.0.jar" %in% list.files(paste0(data_dir, "/graph4lg2_jar")))){

    if(!dir.exists(paths = paste0(data_dir, "/graph4lg2_jar"))){

      dir.create(path = paste0(data_dir, "/graph4lg2_jar"))

    }

    url <- "https://thema.umlp.fr/productions/download.php?name=graphab&version=3.0.0&username=Graph4lg&institution=R"
    #url <- "https://thema.umlp.fr/productions/software/graphab/download/graphab-3.0.0.jar"
    destfile <- "/graph4lg2_jar/graphab-3.0.jar"

    utils::download.file(url, paste0(data_dir, "/", destfile),
                         method = "auto",
                         mode = "wb")

    graphab <- 1

    if(res){
      message("Graphab has been downloaded")
    }

  } else {

    graphab <- 0

    if(res){
      message("Graphab 3.0 is already on your machine")
    }
  }

  if(return){
    return(graphab)
  }

}




