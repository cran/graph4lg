#' Plot histograms of link weights
#'
#' @description The function enables to plot histogram to visualize the
#' distribution of the link weights
#'
#' @param graph A graph object of class \code{igraph} whose links are weighted
#' @param fill A character string indicating the color used to fill
#' the bars (default: "#396D35"). It must be a hexadecimal color code or
#' a color used by default in R.
#' @param class_width (default values: NULL) A numeric or an integer specifying
#' the width of the classes displayed on the histogram. When it is not
#' specified, the width is equal to the difference between the minimum and
#' maximum values divided by 80.
#' @return A ggplot2 object to plot
#' @import ggplot2
#' @export
#' @author P. Savary
#' @examples
#' data(data_ex_genind)
#' mat_w <- mat_gen_dist(data_ex_genind, dist = "DPS")
#' gp <- gen_graph_topo(mat_w = mat_w, topo = "gabriel")
#' hist <- plot_w_hist(gp)


plot_w_hist <- function(graph,
                        fill = "#396D35",
                        class_width = NULL){

  # Check whether 'graph' has weighted links
  if(is.null(igraph::E(graph)$weight)){
    stop("Your graph must have weighted links.")
  }

  # Create a data.frame from the graph links
  graph_df <- data.frame(igraph::as_edgelist(graph))

  # Add the link weights in the data.frame
  graph_df$weight <- igraph::E(graph)$weight
  # Give columns' names to 'graph_df'
  names(graph_df) <- c("from", "to", "weight")

  # Set the class width

  # If not specified, then it depends on the range of values
  if(is.null(class_width)){

    # Get minimum and maximum weights values
    min_w <- min(graph_df$weight)
    max_w <- max(graph_df$weight)

    # Calculate the interval to use in order to have 80 classes
    # between the min and max values
    b_w <- (max_w - min_w)/80

    # Else, it is the specified value, provided it is a numeric or an integer
  } else {

    if(!inherits(class_width, c("numeric", "integer"))){
      stop("'class_width' must be a numeric or an integer")
    } else {
      b_w <- class_width
    }
  }

  # Plot the histogram
  hist <- ggplot(data = graph_df,
                 aes_string(x = 'weight')) +
    geom_histogram(binwidth = b_w, fill = fill, color = "#776F62", size = .2) +
    labs(x="Link weight",
         y="Number of node pairs") +
    theme_bw()

  #print(hist)

  return(hist)

}

