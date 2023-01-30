#' Sample points or patches on a categorical raster layer
#'
#' @description The function samples points or patches on a categorical raster
#' layer.
#'
#' @param raster A RasterLayer object corresponding to a categorical raster layer
#' @param class An integer value or vector with the value(s) corresponding to
#' the code values of the raster layer within which points will be sampled.
#' @param nb_pts An integer value indicating the number of points to be sampled
#' @param dist_min An integer value indicating the minimum distance separating
#' the sampled points (default = 0).
#' @param edge_size An integer value indicating the width of the edge of the
#' raster layer which is ignored during the sampling (default = 0). It prevents
#' from sampling in the margins of the study area.
#' @param by_patch A logical value indicating whether contiguous patches with
#' cells having the same code value are delineated prior to
#' sampling (default = TRUE). It prevents from sampling several points in the
#' same contiguous patch.
#' @param neighborhood An integer value indicating which cells are considered
#' adjacent when contiguous patches are delineated (it should be 8
#' (default, Queen's case) or 4 (Rook's case)). This parameter is ignored
#' when \code{by_patch = FALSE}.
#' @param surf_min An integer value indicating the minimum surface of a patch
#' considered for the sampling in number of raster cells. This parameter is used
#' whatever the \code{by_patch} argument is. Default is 0.
#' @param prop_area A logical value indicating whether sampling in large patches
#' is more likely (default = TRUE). If \code{by_patch = FALSE}, this parameter
#' is ignored. When \code{prop_area = TRUE}, the probability to sample a given
#' patch is proportional to its area.
#' @param step_max An integer value indicating how many sampling steps are
#' performed to identify a point set satisfying all the conditions before
#' returning an error.
#' @param output A character string indicating the type of returned output:
#' \itemize{
#' \item{'data.frame': A \code{data.frame} with three/four columns:
#' \itemize{
#' \item{ID: The point or patch centroid ID}
#' \item{x: The point or patch centroid longitude}
#' \item{y: The point or patch centroid latitude}
#' \item{area: The area of the sampled patch (only if \code{by_patch = TRUE})}
#' }}
#' \item{'pts_layer': A \code{SpatialPointsDataFrame} layer corresponding
#' to the sampled point (points or patch centroids)}
#' \item{'poly_layer': A \code{SpatialPolygonsDataFrame} layer corresponding
#' to the sampled patch polygons}
#' }
#' @param desc A logical value indicating whether the result should be
#' described or not (default = FALSE). If \code{desc = TRUE}, then the Gini
#' coefficient of the distances between points and of the patch areas (if
#' \code{by_patch = TRUE}) is computed with the \code{\link{gini_coeff}}.
#' An histogram of the link weights is also described.
#' @return A list of object(s) with one or several elements according to the
#' \code{output} and \code{desc} arguments.
#' @export
#' @keywords internal
#' @author P. Savary


sample_raster <- function(raster,
                          class,
                          nb_pts,
                          dist_min = 0,
                          edge_size = 0,
                          by_patch = TRUE,
                          neighborhood = 8,
                          surf_min = 0,
                          prop_area = TRUE,
                          step_max = 1000,
                          output = "df",
                          desc = TRUE){


  # ____________________________
  # ____________________________
  # ------- Check the function arguments

  # raster
  if(!inherits(raster, c("RasterLayer"))){
    stop("'RasterLayer' must be either a 'RasterLayer' object.")
  }

  # nb_pts
  if(!inherits(nb_pts, c("integer", "numeric"))){
    stop("'nb_pts' must be either a 'numeric' or 'integer' value.")
  }
  nb_pts <- round(nb_pts, digits = 0)

  # dist_min
  if(!inherits(dist_min, c("integer", "numeric"))){
    stop("'dist_min' must be either a 'numeric' or 'integer' value.")
  }
  dist_min <- round(dist_min, digits = 0)

  # edge_size
  if(!inherits(edge_size, c("integer", "numeric"))){
    stop("'edge_size' must be either a 'numeric' or 'integer' value.")
  }
  edge_size <- round(edge_size, digits = 0)

  # by_patch
  if(!(by_patch %in% c(TRUE, FALSE))){
    # If by_patch is not TRUE nor FALSE, then return an error
    stop("'by_patch' must be either TRUE or FALSE.")
  }

  # neighborhood
  if(!inherits(neighborhood, c("integer", "numeric"))){
    stop("'neighborhood' must be either a 'numeric' or 'integer' value.")
  } else if(!(neighborhood %in% c(4, 8))){
    stop("'neighborhood' must be equal to either 4 or 8.")
  }

  # surf_min
  if(!inherits(surf_min, c("integer", "numeric"))){
    stop("'surf_min' must be either a 'numeric' or 'integer' value.")
  }
  surf_min <- round(surf_min, digits = 0)

  # prop_area
  if(!(prop_area %in% c(TRUE, FALSE))){
    # If prop_area is not TRUE nor FALSE, then return an error
    stop("'prop_area' must be either TRUE or FALSE.")
  }

  # step_max
  if(!inherits(step_max, c("integer", "numeric"))){
    stop("'step_max' must be either a 'numeric' or 'integer' value.")
  }
  step_max <- round(step_max, digits = 0)

  # output
  if(!inherits(output, c("character"))){
    stop("'output' must be either a 'character' string.")
  } else if(!(output %in% c("df", "pts_layer", "poly_layer"))){
    stop("'output' must be equal to either 'df', 'pts_layer' or 'poly_layer'.")
  }

  # desc
  if(!(desc %in% c(TRUE, FALSE))){
    # If desc is not TRUE nor FALSE, then return an error
    stop("'desc' must be either TRUE or FALSE.")
  }

  # ____________________________
  # ____________________________
  # ------ Edge definition

  # Get the extent of the raster
  extent_r <- raster::extent(raster)

  if(edge_size == 0){

    # If edge_size is 0, then sample in the whole raster
    polyg_sample <- methods::as(extent_r, "SpatialPolygons")

  } else {

    # If edge_size != 0, then sample in the raster without its edges
    extent_r@xmin <- extent_r@xmin + edge_size
    extent_r@ymin <- extent_r@ymin + edge_size

    extent_r@xmax <- extent_r@xmax - edge_size
    extent_r@ymax <- extent_r@ymax - edge_size

    polyg_sample <- methods::as(extent_r, "SpatialPolygons")

  }

  # Define the CRS of the polygon in which we sample
  raster::crs(polyg_sample) <- raster::crs(raster)

  # Crop the raster with the polygon
  rast_without_edge <- raster::crop(raster, polyg_sample)

  # ____________________________
  # ____________________________
  # ------ Code values selection

  # Copy the raster and remove values other than the class code
  # (allows for multiple values in class)
  r_class <- rast_without_edge
  r_class[which(!(raster::values(r_class) %in% class))] <- NA

  # Check whether there remains some non-NA values
  # otherwise return an error.
  if(length(unique(raster::values(r_class))) == 1){
    if(is.na(unique(raster::values(r_class)))){
      stop("The 'class' value must be a class code value
           from 'raster'")
    }
  }

  # ____________________________
  # ____________________________
  # ------ Clump to define habitat patches

  # Clump
  r_clump <- raster::clump(r_class, directions = neighborhood)

  # Check whether there is only one patch

  # If it is the case, then r_clump contains either only 1 value
  if(length(unique(raster::values(r_clump))) == 1){
    if(unique(raster::values(r_clump)) == 1){
      message("There is only one patch as all class code values were given
            as 'class' argument.")
      # Error if by_patch is TRUE, given that sampling cannot take place
      # if there is only one patch
      if(by_patch){
        stop("'by_patch' option is not possible if there is only one patch.")
      }
    }
    # or contains 1 and NA values
  } else if(length(unique(raster::values(r_clump))) == 2){

    if(unique(raster::values(r_clump))[2] == 1){
      message("There is only one patch as all class code values were given
            as 'class' argument.")
      # Error if by_patch is TRUE, given that sampling cannot take place
      # if there is only one patch
      if(by_patch){
        stop("'by_patch' option is not possible if there is only one patch.")
      }
    }
  }

  # Get clump values corresponding to the ID of each patch
  val_cl <- as.vector(r_clump)
  # Remove NA values outside patches
  val_cl <- val_cl[-which(is.na(val_cl))]
  # Summarise the data to get the number of pixel per patch
  val_tab <- data.frame(table(val_cl))

  # ____________________________
  # ____________________________
  # ------ Patch selection according to area

  # Remove patches with small surfaces when there is a minimal surface
  if(surf_min != 0){
    val_sup <- val_tab[which(val_tab$Freq < surf_min), ]
    val_sup <- as.numeric(as.character(val_sup$val_cl))
    r_clump[which(raster::values(r_clump) %in% val_sup)] <- NA
  }

  # ____________________________
  # ____________________________
  # ____________________________
  # ____________________________
  # ------ Sampling


  # ____________________________
  # ____________________________
  # ---------- Distinguish cases in which we need only one point per patch (sample polygon)
  # and cases in which it does not matter (sample raster)

  if(by_patch){

    # ____________________________
    # ------ Sample in polygons

    # Convert raster to polygon, dissolved neighboring same values
    r_poly <- raster::rasterToPolygons(r_clump, dissolve = TRUE)

    # We copy it for further use.
    # Row numbers of the polygons in initial r_poly will be polygon IDs later on
    r_poly_copy <- r_poly

    # Convert r_poly to sf multipolygon simple features
    r_poly <- sf::st_as_sf(r_poly)

    # Get a df with the data (ID), the patch centroid coordinates and patch areas
    df_poly <- data.frame(cbind(c(1:nrow(r_poly)),
                                suppressWarnings(sf::st_coordinates(sf::st_centroid(r_poly))),
                                sf::st_area(r_poly)))
    colnames(df_poly) <- c("ID", "x", "y", "area")

  }

  # ____________________________
  # ____________________________
  # ____________________________
  # Sample with the same code whatever 'by_patch' and 'prop_area'

  # ____________________________
  # INITIAL SAMPLING: sample nb_pts/2 patches within df_poly or r_clump

  if(by_patch){
    # ____________________________
    # Sample polygons

    if(prop_area){
      # ____________________________
      # Sampling proportional to area
      prob <- c(df_poly$area/ sum(df_poly$area))

    } else {
      # ____________________________
      # Sampling independent from area
      prob <- c(rep(1, nrow(df_poly))/ nrow(df_poly))
    }

    samp_init <- sample(df_poly$ID,
                        size = round(nb_pts/2, digits = 0),
                        prob = prob,
                        replace = FALSE)
    df_init_samp <- df_poly[samp_init, ]

  } else {

    # ____________________________
    # Sample points
    samp_init <- raster::sampleRandom(r_clump, nb_pts/2, sp = TRUE)
    samp_init <- cbind(data.frame(ID = 1:nrow(samp_init)),
                       data.frame(sp::coordinates(samp_init)))
    colnames(samp_init) <- c("ID", "x", "y")
    df_init_samp <- samp_init
  }

  # ____________________________
  # We compute the pairwise distances between the selected points or centroids
  mat_rast <- suppressMessages(graph4lg::mat_geo_dist(df_init_samp,
                                                      ID = "ID", x = "x", y = "y"))

  # ____________________________
  # Diagonal values are 0 but become dist_min + 1 in order to
  # select different points separated by less than dist_min.
  diag(mat_rast) <- dist_min + 1

  # ____________________________
  # Points too close to each other are removed
  if( length(which(apply(mat_rast, 1, min) < dist_min)) != 0){

    # Get the values of mat_rast
    val_mat <- as.vector(mat_rast)
    # Order by increasing value
    val_mat <- val_mat[order(val_mat)]
    # Get the unique values below dist_min
    val_min <- unique(val_mat[which(val_mat < dist_min)])

    id_min <- c()
    for(n in 1:length(val_min)){
      pn <- data.frame(which(mat_rast == val_min[n], arr.ind = TRUE))[1, ]
      id_min <- c(id_min, row.names(pn))
    }

    df_init_samp <- df_init_samp[-which(df_init_samp$ID %in% id_min), ]
  }

  # ____________________________
  # t is the number of points initially selected
  t <- nrow(df_init_samp)

  # ____________________________
  # Number the loops, initialise to 1
  loop <- 1

  # ____________________________
  # ____________________________
  # WHILE LOOP to add new points while satisfying conditions
  # t will rise to nb_pts * 2 because some points will again be removed after the loop

  # We initialise df_samp
  df_samp <- df_init_samp

  while (t < (round(nb_pts * 1.5, digits = 0)) ){

    if(by_patch){

      # ____________________________
      # Sample polygons

      # We create a data.frame without patches previously sampled
      df_to_samp_step <- df_poly[-which(df_poly$ID %in% df_samp$ID), ]

      if(prop_area){
        # ____________________________
        # Sampling proportional to area
        prob <- c(df_to_samp_step$area/sum(df_to_samp_step$area))
      } else {
        # ____________________________
        # Sampling independent from area
        prob <- c(rep(1, nrow(df_to_samp_step))/nrow(df_to_samp_step))
      }

      # We sample nb_pts/2 points randomly
      samp_step <- sample(df_to_samp_step$ID,
                          size = round(nb_pts/2, digits = 0),
                          prob = prob,
                          replace = FALSE)

      df_step_samp <- df_poly[which(df_poly$ID %in% samp_step), ]

    } else {

      # We sample nb_pts/2 points randomly
      df_step_samp <- raster::sampleRandom(r_clump, nb_pts/2, sp = TRUE)
      df_step_samp <- cbind(data.frame(ID = 1:nrow(df_step_samp)),
                            data.frame(sp::coordinates(df_step_samp)))
      colnames(df_step_samp) <- c("ID", "x", "y")

    }

    # We compute the pairwise distances between the points
    mat_rast <- suppressMessages(graph4lg::mat_geo_dist(df_step_samp,
                                                        ID = "ID", x = "x", y = "y"))

    # Diagonal values are 0 but become dist + 1 in order to
    # select different points separated by less than dist.
    diag(mat_rast) <- dist_min + 1

    # Points too close from each other are removed
    if( length(which(apply(mat_rast, 1, min) < dist_min)) != 0){
      val_mat <- as.vector(mat_rast)
      val_mat <- val_mat[order(val_mat)]
      val_min <- unique(val_mat[which(val_mat < dist_min)])

      id_min <- c()
      for(n in 1:length(val_min)){
        pn <- data.frame(which(mat_rast == val_min[n], arr.ind = TRUE))[1, ]
        id_min <- c(id_min, row.names(pn))
      }
      df_step_samp <- df_step_samp[-which(df_step_samp$ID %in% id_min), ]
    }

    # ____________________________
    # Append the final data.frame
    df_samp <- rbind(df_samp, df_step_samp)
    df_samp$ID <- 1:nrow(df_samp)

    # We compute the pairwise distances between the points
    mat_rast <- suppressMessages(graph4lg::mat_geo_dist(df_samp,
                                                        ID = "ID", x = "x", y = "y"))

    # Diagonal values are 0 but become dist + 1 in order to
    # select different points separated by less than dist.
    diag(mat_rast) <- dist_min + 1

    # Points too close from each other are removed
    if( length(which(apply(mat_rast, 1, min) < dist_min)) != 0){
      val_mat <- as.vector(mat_rast)
      val_mat <- val_mat[order(val_mat)]
      val_min <- unique(val_mat[which(val_mat < dist_min)])

      if(0 %in% val_min){
        df_samp <- df_samp[-which(duplicated(df_samp$x, df_samp$y)), ]
        val_min <- val_min[-which(val_min == 0)]
      }

      if(length(val_min) > 0){
        id_min <- c()
        for(n in 1:length(val_min)){
          pn <- data.frame(which(mat_rast == val_min[n], arr.ind = TRUE))[1, ]
          id_min <- c(id_min, row.names(pn))
        }

        df_samp <- df_samp[-which(df_samp$ID %in% id_min), ]
      }
    }

    # t is the number of points selected
    t <- nrow(df_samp)

    # Add one to loop
    loop <- loop + 1

    # Print the loop ID
    #print(paste0("Loop ", loop))


    if(loop >= (1 + step_max)){
      stop(paste0("No solution was obtained over ", step_max, " steps."))
    }

  }

  #_____________________________________
  # Sample to select only nb_pts points
  df_samp <- df_samp[sample(1:nrow(df_samp), nb_pts), ]
  df_samp$ID <- 1:nrow(df_samp)

  #_____________________________________
  #_____________________________________
  #_____________________________________
  # Descriptions

  if(desc){

    if(by_patch){
      gini_area <- gini_coeff(df_samp$area)
    }

    mat_dist <- suppressMessages(graph4lg::mat_geo_dist(data = df_samp,
                                                        ID = "ID",
                                                        x = "x", y = "y"))

    gini_dist <- gini_coeff(ecodist::lower(mat_dist))

    g <- graph4lg::gen_graph_topo(mat_w = mat_dist, topo = "comp")

    hist <- graph4lg::plot_w_hist(g) +
      ggplot2::labs(x = "Euclidean distance",
                    y = "Number of point pairs")
    #print(hist)

    if(by_patch){
      description <- list(gini_area, gini_dist, hist)
      names(description) <- c("Gini coefficient for patch areas",
                              "Gini coefficient for Euclidean distances",
                              "Histogram of Euclidean distances between points")
    } else {
      description <- list(gini_dist, hist)
      names(description) <- c("Gini coefficient for Euclidean distances",
                              "Histogram of Euclidean distances between points")
    }


  }

  #_____________________________________
  #_____________________________________
  #_____________________________________
  # Outputs


  if (output == "df"){
    #_____________________________________
    # If output is a data.frame, then it is df_samp
    out <- df_samp

  } else if(output == "pts_layer"){

    #_____________________________________
    # If output is a points layer, then we create it

    xy <- df_samp[,c('x','y')]
    mxy <- as.matrix(xy)

    list_pts <- list()
    for(i in 1:nrow(xy)){
      list_pts[[i]] <- sf::st_point(mxy[i, ])
    }
    pts_lay <- sf::st_sfc(list_pts)

    pts_lay <- sf::st_sf(pts_lay,
                         df_samp[, c("ID")])

    pts_lay_spat <- sf::as_Spatial(pts_lay, IDs = "ID")

    # Add df_samp as attribute table
    pts_lay_spat@data <- df_samp

    # Define CRS
    raster::crs(pts_lay_spat) <- raster::crs(raster)

    out <- pts_lay_spat

  } else if(output == "polygon_layer"){

    #_____________________________________
    # If output is a polygon layer, we create it from r_poly_copy saved earlier
    poly_lay <- r_poly_copy[df_samp$ID, ]

    # We add df_samp as attribute table
    poly_lay@data <- df_samp

    out <- poly_lay

  }

  if(desc){
    res <- list(out, description)
  } else {
    res <- out
  }

  return(res)

}

