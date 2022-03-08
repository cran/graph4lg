#' Describe the objects of a Graphab project
#'
#' @description The function describes the objects of a Graphab project
#'
#' @inheritParams get_graphab_linkset
#' @param mode A character string indicating the objects of the project that
#' are described. It must be either:\itemize{
#' \item{\code{mode='patches'}(default): The habitat patches are described
#' with synthetic descriptors (code, number, mean capacity, median capacity,
#' capacity harmonic mean, capacity Gini coefficient) and a histogram of
#' capacity distribution.}
#' \item{\code{mode='linkset'}: The links of a link set are described
#' with synthetic descriptors (codes, costs, number, mean cost distance,
#' median cost distance, cost distance harmonic mean, cost distance Gini
#' coefficient) and a histogram of cost distance distribution.}
#' \item{\code{mode='both'}: Both the patches and links of a linkset are
#' described}
#' }
#' @param fig Logical (default = FALSE) indicating whether to plot a figure of
#' the resulting spatial graph. The figure is plotted using function
#' \code{\link{plot_graph_lg}}. The plotting can be long if the graph has many
#' nodes and links.
#' @param return_val Logical (default = TRUE) indicating whether the project
#' features are returned as a list (TRUE) or only displayed in the
#' R console (FALSE).
#' @import ggplot2
#' @export
#' @author P. Savary
#' @examples
#' \dontrun{
#' graphab_project_desc(proj_name = "grphb_ex",
#'                      mode = "patches",
#'                      fig = FALSE)
#' }


graphab_project_desc <- function(proj_name,
                                 mode = "patches",
                                 linkset = NULL,
                                 proj_path = NULL,
                                 fig = FALSE,
                                 return_val = TRUE){


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

  ##########################################
  # Check for mode
  if(!inherits(mode, "character")){
    stop("'mode' must be a character string")
  } else if(!(mode %in% c("patches", "linkset", "both"))){
    stop("'mode' must be equal to either 'patches', 'linkset' or 'both'.")
  }


  ##########################################
  # Check for return_val
  if(!inherits(return_val, "logical")){
    stop("'return_val' must be a logical")
  }


  #####
  # Describe patches
  if(mode %in% c("patches", "both")){

    # Get all codes
    all_codes <- get_graphab_raster_codes(proj_name = proj_name,
                                          mode = 'all',
                                          proj_path = proj_path)
    # Get habitat code
    hab_code <- get_graphab_raster_codes(proj_name = proj_name,
                                         mode = 'habitat',
                                         proj_path = proj_path)
    # Get patches information
    patches <- get_graphab_metric(proj_name = proj_name,
                                  proj_path = proj_path)

    all_codes_p <- paste0("Raster source layer codes: ",
                          paste(all_codes, collapse = ", "))
    hab_codes_p <- paste0("Habitat patch codes: ",
                          paste(hab_code, collapse = ", "))
    nb_p <- paste0("Number of patches: ", nrow(patches))
    cap <- "Patch capacities:"
    total_cap <- paste0("  Total: ", sum(patches$Capacity))
    min_cap <- paste0("  Minimum: ", min(patches$Capacity))
    max_cap <- paste0("  Maximum: ", max(patches$Capacity))
    mean_cap <- paste0("  Mean: ", mean(patches$Capacity))
    median_cap <- paste0("  Median: ", stats::median(patches$Capacity))
    sd_cap <- paste0("  Standard deviation: ", stats::sd(patches$Capacity))
    hm_cap <- paste0("  Harmonic mean: ", harm_mean(patches$Capacity))
    gini_cap <- paste0("  Gini index: ", gini_coeff(patches$Capacity))

    cat(c(all_codes_p, "\n",
          hab_codes_p, "\n",
          nb_p, "\n",
          cap, "\n",
          total_cap, "\n",
          min_cap, "\n",
          max_cap, "\n",
          mean_cap, "\n",
          median_cap, "\n",
          sd_cap, "\n",
          hm_cap, "\n",
          gini_cap, "\n"),
        "\n", sep = "")

    if(return_val){
      res_p <- list(all_codes, hab_code,
                    nrow(patches), sum(patches$Capacity),
                    min(patches$Capacity), max(patches$Capacity),
                    mean(patches$Capacity), stats::median(patches$Capacity),
                    stats::sd(patches$Capacity), harm_mean(patches$Capacity),
                    gini_coeff(patches$Capacity))
      names(res_p) <- c("Raster source layer codes",
                        "Habitat patch codes",
                        "Number of patches",
                        "Total patch capacities",
                        "Min. patch capacity",
                        "Max. patch capacity",
                        "Mean patch capacity",
                        "Median patch capacity",
                        "Std. deviation patch capacity",
                        "Harmonic mean patch capacity",
                        "Gini coeff. patch capacity")
    }


    if(fig){
      range <- max(patches$Capacity) - min(patches$Capacity)

      if(range > 0){
        b_w <- range/80

        fig_patch <- ggplot(data = patches,
                            aes(x = .data$Capacity)) +
          geom_histogram(binwidth = b_w,
                         fill = "#396D35",
                         color = "#776F62", size = .2) +
          labs(x = "Patch capacities",
               y = "Frequencies")

        print(fig_patch)

      } else {
        message(paste0("The range of patch capacities is equal to 0. ",
                       "Plotting an histogram does not make sense."))
        fig_patch <- NULL
      }

    } else {

      fig_patch <- NULL
    }

  }


  if(mode %in% c("linkset", "both")){

    # Get linkset
    linkset_desc <- get_graphab_linkset_cost(proj_name = proj_name,
                                             linkset = linkset,
                                             proj_path = proj_path)
    # Get linkset values
    all_links <- get_graphab_linkset(proj_name = proj_name,
                                     linkset = linkset,
                                     proj_path = proj_path)

    all_codes_l <- paste0("Raster source layer codes considered for the costs: ",
                          paste(linkset_desc$code, collapse = ", "))
    all_cost_l <- paste0("Corresponding costs: ",
                         paste(linkset_desc$cost, collapse = ", "))
    nb_l <- paste0("Number of links: ", nrow(all_links))
    cd <- "Cost distances:"
    min_cd <- paste0("  Minimum: ", min(all_links$Dist))
    max_cd <- paste0("  Maximum: ", max(all_links$Dist))
    mean_cd <- paste0("  Mean: ", mean(all_links$Dist))
    median_cd <- paste0("  Median: ", stats::median(all_links$Dist))
    sd_cd <- paste0("  Standard deviation: ", stats::sd(all_links$Dist))
    hm_cd <- paste0("  Harmonic mean: ", harm_mean(all_links$Dist))
    gini_cd <- paste0("  Gini index: ", gini_coeff(all_links$Dist))

    cat(c(all_codes_l, "\n",
          all_cost_l, "\n",
          cd, "\n",
          nb_l, "\n",
          min_cd, "\n",
          max_cd, "\n",
          mean_cd, "\n",
          median_cd, "\n",
          sd_cd, "\n",
          hm_cd, "\n",
          gini_cd, "\n"),
        "\n", sep = "")


    if(return_val){
      res_l <- list(linkset_desc$code, linkset_desc$cost,
                    nrow(all_links),
                    min(all_links$Dist), max(all_links$Dist),
                    mean(all_links$Dist), stats::median(all_links$Dist),
                    stats::sd(all_links$Dist), harm_mean(all_links$Dist),
                    gini_coeff(all_links$Dist))
      names(res_l) <- c("Raster source layer codes considered for the costs",
                        "Corresponding costs",
                        "Number of links",
                        "Min. cost distance",
                        "Max. cost distance",
                        "Mean cost distance",
                        "Median cost distance",
                        "Std. deviation cost distance",
                        "Harmonic mean cost distance",
                        "Gini coeff. cost distance")
    }


    if(fig){
      range <- max(all_links$Dist) - min(all_links$Dist)

      if(range > 0){
        b_w <- range/80

        fig_cd <- ggplot(data = all_links,
                         aes(x = .data$Dist)) +
          geom_histogram(binwidth = b_w,
                         fill = "#396D35",
                         color = "#776F62", size = .2) +
          labs(x = "Cost distances",
               y = "Frequencies")

        print(fig_cd)

      } else {
        message(paste0("The range of cost distances is equal to 0. ",
                       "Plotting an histogram does not make sense."))
        fig_cd <- NULL
      }

    } else {

      fig_cd <- NULL
    }

  }


  if(return_val){
    if(mode == "both"){
      res <- list(res_p, res_l)
    } else if(mode == "patches"){
      res <- res_p
    } else if(mode == "linkset"){
      res <- res_l
    }
    return(res)
  }

}



