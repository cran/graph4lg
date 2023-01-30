#' Compare two link sets created in a Graphab project
#'
#' @description The function compares two link sets created in a Graphab project
#' both quantitatively and spatially.
#'
#' @param proj_name A character string indicating the Graphab project name.
#' The project name is also the name of the project directory in which the
#' file proj_name.xml is. It can be created with \code{\link{graphab_project}}
#' @param linkset1 A character string indicating the name of the first link set
#' involved in the comparison. The link set has to be present in the project
#' and can be created with \code{\link{graphab_link}}.
#' @param linkset2 A character string indicating the name of the second link set
#' involved in the comparison. The link set has to be present in the project
#' and can be created with \code{\link{graphab_link}}.
#' @param buffer_width (default=200) An integer or numeric value
#' indicating the width of the buffer created in each side of the links prior
#' to spatial intersection. It is expressed in meters.
#' @param min_length (optional, default=NULL) An integer or numeric value
#' indicating the minimum length in meters of the links to be compared. Links
#' whose length is larger than \code{min_length} will be ignored in
#' the comparison.
#' @param proj_path (optional, default=NULL) A character string indicating the
#' path to the directory that contains the project directory. It should be used
#' when the project directory is not in the current working directory.
#' Default is NULL. When 'proj_path = NULL', the project directory is equal
#' to \code{getwd()}.
#' @details The function compares two link sets linking the same habitat patches
#' of the Graphab project but computed using different cost scenarios. It
#' creates a buffer in each side of every link and then overlaps every link
#' in linkset1 with the same link in linkset2. It returns the area of both
#' buffered links and the area of their intersection. It also computes the
#' Mantel correlation coefficient between the cost distances associated to the
#' same links in both linksets.
#' @export
#' @keywords internal
#' @author P. Savary
#' @examples
#' \dontrun{
#' link_compar(proj_name = "grphb_ex",
#'               linkset1 = "lcp1",
#'               linkset2 = "lcp2"
#'               buffer_width = 200)
#' }


link_compar <- function(proj_name,         # character
                        linkset1, # character
                        linkset2, # character
                        buffer_width = 200,
                        min_length = NULL,
                        proj_path = NULL){ # if NULL getwd() otherwise a character path

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

  #########################################
  # Check for linkset1 class
  if(!inherits(linkset1, "character")){
    stop("'linkset1' must be a character string specifying the name of the
         first link set involved in the comparison.")
  } else if (!(paste0(linkset1, "-links.csv") %in%
               list.files(path = paste0(proj_path, "/", proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  #########################################
  # Check for linkset2 class
  if(!inherits(linkset2, "character")){
    stop("'linkset2' must be a character string specifying the name of the
         second link set involved in the comparison.")
  } else if (!(paste0(linkset2, "-links.csv") %in%
               list.files(path = paste0(proj_path, "/", proj_name)))){
    stop("The linkset you refer to does not exist.
           Please use graphab_link() before.")
  }

  #########################################
  # Check for buffer_width
  if(!inherits(buffer_width, c("numeric", "integer"))){
    stop("'buffer_width' must be a numeric or integer value")
  }


  #########################################
  # Check for min_length
  if(!is.null(min_length)){
    if(!inherits(min_length, c("numeric", "integer"))){
      stop("'min_length' must be a numeric or an integer threshold value.")
    }
  }


  ###########################################
  # Open layers corresponding to linksets

  ls1 <- suppressWarnings(sf::as_Spatial(sf::st_read(dsn = paste0(proj_path,
                                                                  "/", proj_name),
                                                     layer = paste0(linkset1,
                                                                    "-links"))))

  ls2 <- suppressWarnings(sf::as_Spatial(sf::st_read(dsn = paste0(proj_path,
                                                                  "/", proj_name),
                                                     layer = paste0(linkset2,
                                                                    "-links"))))


  ###########################################
  # Filter to limit the number of links by removing short links

  if(!is.null(min_length)){
    ls1 <- ls1[which(ls1$DistM >= min_length), ]
    ls2 <- ls2[which(ls2$DistM >= min_length), ]
  }

  ###########################################
  # Check that they share the same links

  if(nrow(ls1) != nrow(ls2)){
    stop("'linkset1' and 'linkset2' must have the same number of links.")
  } else if(!(all(ls1$Id %in% ls2$Id))){
    stop("'linkset1' and 'linkset2' must have the same link IDs.")
  }


  ###########################################
  # Add buffer
  print(paste0("The buffer width on each side of the links has been set to ",
               buffer_width, " m."))

  ls1_b <- raster::buffer(ls1, width = buffer_width, dissolve = FALSE)
  ls2_b <- raster::buffer(ls2, width = buffer_width, dissolve = FALSE)


  data1 <- ls1_b@data
  data2 <- ls2_b@data

  #########################################################################
  #########################################################################

  # ids <- ls1_b$Id
  #
  # spat_over <- function(id){
  #   ls1_bi <- ls1_b[which(data1$Id == id), ]
  #   ls2_bi <- ls2_b[which(data2$Id == id), ]
  #
  #   area_1 <- sf::st_area(sf::st_as_sf(ls1_bi))
  #   area_2 <- sf::st_area(sf::st_as_sf(ls2_bi))
  #
  #   inter_ls <- suppressWarnings(sf::st_intersection(sf::st_as_sf(ls1_bi),
  #                                                    sf::st_as_sf(ls2_bi)))
  #   inter_area <- sf::st_area(inter_ls)
  #
  #   return(list(area_1, area_2,
  #               inter_area,
  #               ls1_bi$Dist, ls2_bi$Dist,
  #               ls1_bi$DistM, ls2_bi$DistM))
  #
  # }
  # res_lap <- lapply(ids, FUN = spat_over)
  #
  # df_res <- data.frame(id_link = ids,
  #                      area_1 = unlist(lapply(res_lap, "[", 1)),
  #                      area_2 = unlist(lapply(res_lap, "[", 2)),
  #                      cost_dist_1 = unlist(lapply(res_lap, "[", 4)),
  #                      cost_dist_2 = unlist(lapply(res_lap, "[", 5)),
  #                      euc_dist_1 = unlist(lapply(res_lap, "[", 6)),
  #                      euc_dist_2 = unlist(lapply(res_lap, "[", 7)),
  #                      area_overlap = unlist(lapply(res_lap, "[", 3)))

  #########################################################################
  #########################################################################


  df_res <- data.frame(id_link = NA,
                       area_1 = NA,
                       area_2 = NA,
                       cost_dist_1 = NA,
                       cost_dist_2 = NA,
                       euc_dist_1 = NA,
                       euc_dist_2 = NA,
                       area_overlap = NA)
  df_res <- df_res[-1, ]


  for(i in 1:nrow(ls1_b)){

    id <- data1[i, 'Id']

    ls1_bi <- ls1_b[which(data1$Id == id), ]
    ls2_bi <- ls2_b[which(data2$Id == id), ]

    inter_ls <- suppressWarnings(sf::st_intersection(sf::st_as_sf(ls1_bi),
                                                     sf::st_as_sf(ls2_bi)))
    if(nrow(inter_ls) != 0){
      inter_area <- sf::st_area(inter_ls)
    } else {
      inter_area <- 0
    }

    ls1_area <- sf::st_area(sf::st_as_sf(ls1_bi))
    ls2_area <- sf::st_area(sf::st_as_sf(ls2_bi))

    df_res <- rbind(df_res,
                    data.frame(id_link = id,
                               area_1 = as.numeric(ls1_area),
                               area_2 = as.numeric(ls2_area),
                               cost_dist_1 = ls1_bi$Dist,
                               cost_dist_2 = ls2_bi$Dist,
                               euc_dist_1 = ls1_bi$DistM,
                               euc_dist_2 = ls2_bi$DistM,
                               area_overlap = as.numeric(inter_area)))

  }

  correl <- stats::cor(df_res$cost_dist_1, df_res$cost_dist_2)

  res <- list(df_res, correl)
  names(res) <- c("Spatial overlap table",
                  "Correlation coefficient between cost distances")

  return(res)

}

