#' Compute cost distances between points on a raster
#'
#' @description The function computes cost-distances associated to least cost
#' paths between point pairs on a raster with specified cost values.
#'
#' @param raster A parameter indicating the raster file on which cost distances
#' are computed. It can be:\itemize{
#' \item{A character string indicating the path to a raster file in format
#' .tif or .asc.}
#' \item{A \code{RasterLayer} object already loaded in R environment}
#' }
#' All the raster cell values must be present in the column 'code' from
#' \code{cost} argument.
#' @param pts A parameter indicating the points between which cost distances
#' are computed. It can be either: \itemize{
#' \item{A character string indicating the path to a .csv file. It must have
#' three columns:\itemize{
#' \item{ID: The ID of the points.}
#' \item{x: A numeric or integer indicating the longitude of the points.}
#' \item{y: A numeric or integer indicating the latitude of the points.}
#' }}
#' \item{A \code{data.frame} with the spatial coordinates of the points.
#' It must have three columns:\itemize{
#' \item{ID: The ID of the points.}
#' \item{x: A numeric or integer indicating the longitude of the points.}
#' \item{y: A numeric or integer indicating the latitude of the points.}
#' }}
#' \item{A \code{SpatialPointsDataFrame} with at least an attribute column
#' named "ID" with the point IDs.}
#' }
#' The point coordinates must be in the same spatial coordinate reference system
#' as the raster file.
#' @param cost A \code{data.frame} indicating the cost values associated to each
#' raster value. It must have two columns:\itemize{
#' \item{'code': raster cell values}
#' \item{'cost': corresponding cost values}
#' }
#' @param method A character string indicating the method used to compute the
#' cost distances. It must be:\itemize{
#' \item{'gdistance': uses the functions from the package \pkg{gdistance}
#' assuming that movement is possible in 8 directions from each cell, that
#' a geo-correction is applied to correct for diagonal movement lengths and that
#' raster cell values correspond to resistance (and not conductance).}
#' \item{'java': uses a .jar file which is downloaded on the user's machine if
#' necessary and if java is installed. This option substantially reduces
#' computation times and makes possible the parallelisation.}
#' }
#' @param direction An integer (4, 8, 16) indicating the directions in which
#' movement can take place from a cell. Only used when \code{method="gdistance"}.
#' By default, \code{direction=8}.
#' @param parallel.java An integer indicating how many computer cores are used
#' to run the .jar file. By default, \code{parallel.java=1}.
#' @param alloc_ram (optional, default = NULL) Integer or numeric value
#' indicating RAM gigabytes allocated to the java process when used. Increasing
#' this value can speed up the computations. Too large values may not be
#' compatible with your machine settings.
#' @param return A character string indicating whether the returned object is a
#' \code{data.frame} (\code{return="df"}) or a pairwise
#' \code{matrix} (\code{return="mat"}).
#' @return The function returns:\itemize{
#' \item{If \code{return="mat"}, a pairwise \code{matrix} with cost-distance
#' values between points.}
#' \item{If \code{return="df"}, an object of type \code{data.frame} with three columns:
#' \itemize{
#' \item{from: A character string indicating the ID of the point of origin.}
#' \item{to: A character string indicating the ID of the point of destination.}
#' \item{cost_dist: A numeric indicating the accumulated cost-distance along
#' the least-cost path between point ID1 and point ID2.}
#' }
#' }
#' }
#' @export
#' @author P. Savary
#' @examples
#' x <- raster::raster(ncol=10, nrow=10, xmn=0, xmx=100, ymn=0, ymx=100)
#' raster::values(x) <- sample(c(1,2,3,4), size = 100, replace = TRUE)
#' pts <- data.frame(ID = 1:4,
#'                   x = c(10, 90, 10, 90),
#'                   y = c(90, 10, 10, 90))
#' cost <- data.frame(code = 1:4,
#'                    cost = c(1, 10, 100, 1000))
#' mat_cost_dist(raster = x,
#'               pts = pts, cost = cost,
#'               method = "gdistance")


mat_cost_dist <- function(raster,
                          pts,
                          cost,
                          method = "gdistance",
                          return = "mat",
                          direction = 8,
                          parallel.java = 1,
                          alloc_ram = NULL){


  # Check raster argument

  if(inherits(raster, "RasterLayer")){

    r_type <- "RasterLayer"

  } else if (file.exists(raster)){

    r_type <- "RasterFile"

    raster_path <- raster

  } else {

    stop("'raster' must be either a RasterLayer object or a valid path
           to a .tif or .asc raster layer file.")
  }

  # Check pts argument

  if(inherits(pts, "SpatialPointsDataFrame")){

    if(!("ID" %in% colnames(pts@data))){
      stop("'pts' must include an 'ID' column")
    }

    p_type <- "SpatialPointsDataFrame"

  } else if (inherits(pts, "data.frame")){

    if(all(c('ID', 'x', 'y') %in% colnames(pts))){

      p_type <- "df"

      pts <- pts[, c('ID', 'x', 'y')]

    } else {

      stop("The columns of pts must include 'ID', 'x' and 'y'")

    }

  } else if (inherits(pts, "character")){

    if(file.exists(pts)){

      pts_path <- pts

      pts <- utils::read.csv(file = pts_path)

      if(all(c('ID', 'x', 'y') %in% colnames(pts))){

        p_type <- "csv"

      } else {

        stop("The columns of pts must include 'ID', 'x' and 'y'")

      }
    }

  } else {

    stop("'pts' must be either a SpatialPointsDataFrame object or a data.frame,
          or a valid path to a .csv file")
  }

  idp <- as.character(pts$ID)

  # Check cost argument

  if(!inherits(cost, "data.frame")){
    stop("'cost' must be a data.frame object")
  } else {
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

    c_type = "data.frame"
  }


  # Check return argument
  if(!inherits(return, "character")){
    stop("'return' must be a character string.")
  } else if(!(return %in% c("df", "mat"))){
    stop("'return' must be either 'df' or 'mat'.")
  }


  # Method : gdistance

  if(method == "gdistance"){

    ############################
    # Create an ascii file
    if(r_type == "RasterLayer"){

      rast_val <- unique(raster::values(raster))

      if(any(is.na(rast_val))){
        rast_val <- rast_val[-which(is.na(rast_val))]
      }

      if(!all(rast_val %in% cost$code)){
        stop("Specify the cost value associated to every raster cell value")
      }

    } else if (r_type == "RasterFile"){

      raster <- raster::raster(x = raster_path)

      rast_val <- unique(raster::values(raster))

      if(any(is.na(rast_val))){
        rast_val <- rast_val[-which(is.na(rast_val))]
      }

      if(!all(rast_val %in% cost$code)){
        stop("Specify the cost value associated to every raster cell value")
      }


    }
    ############################

    if(p_type == "csv"){

      pts <- utils::read.csv(file = pts_path)

      if(!all(c('ID', 'x', 'y') %in% colnames(pts))){

        stop("The columns of pts must include 'ID', 'x' and 'y'")

      }

    }

    if(p_type %in% c("df", "csv")){
      pts <- suppressWarnings(sp::SpatialPointsDataFrame(coords = pts[, c('x', 'y')], data = pts))
    }

    # Raster reclass

    mat_rcl <- matrix(c(cost$code, cost$cost), nrow = nrow(cost), ncol = 2)
    raster <- raster::reclassify(raster, rcl = mat_rcl)

    reso <- raster::res(raster)[1]

    trans <- gdistance::transition(x = raster,
                                   transitionFunction = function(x) 1/mean(x),
                                   directions = direction)
    trans <- gdistance::geoCorrection(trans)
    cost_dist <- gdistance::costDistance(x = trans,
                                         fromCoords = pts,
                                         toCoords = pts)

    cost_dist <- cost_dist/reso
    row.names(cost_dist) <- colnames(cost_dist) <- idp

    cost_dist <- graph4lg::pw_mat_to_df(pw_mat = cost_dist)

    colnames(cost_dist) <- c("from", "to", "id_link", "cost_dist")

    cost_dist <- cost_dist[, c("from", "to", "cost_dist")]


    # Method : java

  } else if (method == "java"){

    if(Sys.which("java") == ""){
      stop("Please install java to use 'method == 'java''")
    }

    java.path <- Sys.which("java")

    # Check for costdist.jar and download it if necessary

    data_dir <- rappdirs::user_data_dir()

    if("costdist-0.3.jar" %in% list.files(paste0(data_dir, "/graph4lg_jar"))){

      message("costdist-0.3.jar will be used")

    } else {

      message("costdist-0.3.jar will be downloaded")

      if(!dir.exists(paths = paste0(data_dir, "/graph4lg_jar"))){

        dir.create(path = paste0(data_dir, "/graph4lg_jar"))

      }

      url <- "https://thema.univ-fcomte.fr/productions/download.php?name=graphab&prog=costdist&version=0.3&username=Graph4lg&institution=R"
      #url <- "https://sourcesup.renater.fr/www/graphab/download/costdist-0.3.jar"

      destfile <- "/graph4lg_jar/costdist-0.3.jar"

      utils::download.file(url, paste0(data_dir, "/", destfile),
                    method = "auto",
                    mode = "wb")

    }

    # Raster

    # Create an ascii file
    if(r_type == "RasterLayer"){

      file_rast <- tempfile(fileext = ".asc")
      raster::writeRaster(raster,
                          file = file_rast,
                          overwrite = TRUE)
      del_rast <- 1

    } else if (r_type == "RasterFile"){

      if(stringr::str_detect(raster_path, pattern = ".tif")){

        file_rast <- tempfile(fileext = ".asc")
        raster <- raster::raster(x = raster_path)
        raster::writeRaster(raster,
                            file = file_rast,
                            overwrite = TRUE)

        del_rast <- 1

      } else {
        file_rast <- raster_path

        del_rast <- 0
      }
    }

    # Point

    if(p_type == "SpatialPointsDataFrame"){

      pts <- cbind(pts@data[, 'ID'], pts@coords)
      colnames(pts) <- c('ID', 'x', 'y')

    } else if (p_type == "csv"){

      pts <- utils::read.csv(file = pts_path)
      pts <- pts[, c('ID', 'x', 'y')]

    }

    file_pts <- tempfile(fileext = ".csv")
    utils::write.csv(pts, file = file_pts, row.names = FALSE)


    # Cost values argument
    ncode <- nrow(cost)

    vec_cost <- c()
    for(i in 1:ncode){
      vec_cost <- c(vec_cost, paste0(cost[i, "code"], "=", cost[i, 'cost']))

    }

    file_res <- tempfile(fileext = ".txt")

    # Run java code

    cmd <- c("-Djava.awt.headless=true", "-jar",
             paste0(data_dir, "/graph4lg_jar/costdist-0.3.jar"),
             parallel.java, file_pts, file_rast, file_res, vec_cost)


    if(!is.null(alloc_ram)){
      if(inherits(alloc_ram, c("integer", "numeric"))){
        cmd <- c(paste0("-Xmx", alloc_ram, "g"), cmd)
      } else {
        stop("'alloc_ram' must be a numeric or an integer")
      }
    }


    system2(java.path, args = cmd)

    ###########################################################################################
    # Open res

    cost_dist <- utils::read.table(file = file_res, header = FALSE, sep = ",")

    colnames(cost_dist) <- c("from", "to", "cost_dist")

    file.remove(file_res)
    file.remove(file_pts)

    if(del_rast == 1){
      file.remove(file_rast)
    }
  } else {
    stop("'method' must be either 'gdistance' or 'java'")
  }

  cost_dist <- cost_dist[order(cost_dist$from, cost_dist$to), ]
  row.names(cost_dist) <- paste0(cost_dist$from, "_", cost_dist$to)

  if(return == "mat"){
    cost_dist <- df_to_pw_mat(data = cost_dist,
                              from = "from", to = "to",
                              value = "cost_dist")
  }

  return(cost_dist)

}

