#' Prune a graph using the 'percolation threshold' method
#'
#' @description The function allows to prune a graph by removing
#' the links with the largest weights until the graph breaks into
#' two components. The returned graph is the last graph with only one
#' component.
#'
#' @param x A symmetric \code{matrix} of pairwise distances between nodes
#' @param val_step The number of classes to create to search for the
#' threshold value without testing all the possibilities. By default,
#' 'val_step = 20'.
#' @return A graph object of type \code{igraph}
#' @export
#' @author P. Savary
#' @examples
#' data(data_simul_genind)
#' suppressWarnings(mat_w <- graph4lg::mat_geo_dist(data=pts_pop_simul,
#'                             ID = "ID",
#'                             x = "x",
#'                             y = "y"))
#' g_percol(x = mat_w)

g_percol <- function(x, val_step = 20){

  # Check whether x is a symmetric matrix
  if(class(x) != "matrix"){
    stop("'x' must be a matrix")
  } else if(!isSymmetric(x)){
    stop("The matrix 'x' must be symmetric")
  } else {
    # Creation of the complete initial graph
    g1 <- igraph::graph.adjacency(x,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
  }

  # Edge list of the complete graph
  g1_df <- data.frame(igraph::as_edgelist(g1))
  # We add the weight of each link in the edge list data frame
  g1_df$w <- igraph::E(g1)$weight

  # We order the df with the larger weights first
  g1_df <- g1_df[order(g1_df$w, decreasing = TRUE), ]

  # We calculate the limits of the classes (id_sup)
  val_part <- round(nrow(g1_df)/val_step, digits = 0)
  # We create id_sup. Its first element is 1
  id_sup <- c(1)
  for (i in 1:(val_step-1)){
    id_sup[i+1] <- val_part*i
  }
  id_sup[val_step] <- nrow(g1_df)-1
  # The first value of id_sup is 1, the last is the number of edges - 1

  # We calculate the weight of the link at the class limits (threshold)
  threshold <- g1_df[id_sup, 'w']

  # We will find in which class, weights become so small that the graph breaks into
  # two components

  # If we remove the links with a big weight, the graph stays connected logically
  # Once it breaks into two components, it means that the threshold is
  # in the last class before the value

  comp <- 1
  i = 1
  while(comp < 2){
    mat_g2 <- x
    mat_g2[mat_g2 > threshold[i]] <- 0
    g2 <- igraph::graph.adjacency(mat_g2,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
    comp <- igraph::components(g2)$no[1]
    i <- i + 1
  }
  t1 <- i - 1

  # t1 is the first iteration for which we get two components
  # So threshold[t1] is lower than the real threshold
  # The real threshold is between the id_sup[t1-1]-th and the id_sup[t1]-th
  # rows of g1_df

  # We reiterate the operation in the identified class
  val_step2 <- val_step/2
  val_part2 <- round(val_part/val_step2, digits = 0)
  id_sup2 <- c(id_sup[t1 - 1])
  for (i in 1:(val_step2-1)){
    id_sup2[i+1] <- id_sup[t1-1]+val_part2*i
  }
  id_sup2[val_step2] <- id_sup[t1]

  # We calculate the weight of the link at the class limits (threshold)
  threshold2 <- g1_df[id_sup2, 'w']

  # We will find in which class, weights become so small that the graph breaks into
  # two components

  # If we remove the links with a big weight, the graph stays connected logically
  # Once it breaks into two components, it means that the threshold is
  # in the last class before the value

  comp <- 1
  i = 1
  while(comp < 2){
    mat_g2 <- x
    mat_g2[mat_g2 > threshold2[i]] <- 0
    g2 <- igraph::graph.adjacency(mat_g2,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
    comp <- igraph::components(g2)$no[1]
    i <- i + 1
  }
  t2 <- i - 1

  # The threshold is the weight of one of the link
  # between the rows id_sup2[t2 - 1] and id_sup2[t2]
  # of g1_df

  id_sup3 <- c(id_sup2[t2-1]:id_sup2[t2])
  threshold3 <- g1_df[id_sup3, 'w']

  comp <- 1
  i = 1
  while(comp < 2){
    mat_g2 <- x
    mat_g2[mat_g2 > threshold3[i]] <- 0
    g2 <- igraph::graph.adjacency(mat_g2,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
    comp <- igraph::components(g2)$no[1]
    i <- i + 1
  }
  t3 <- i - 1

  thr2 <- g1_df[id_sup3[t3], 'w']
  mat_g2 <- x
  mat_g2[mat_g2 > thr2] <- 0
  g2 <- igraph::graph.adjacency(mat_g2,
                                mode = "undirected",
                                weighted = TRUE,
                                diag = FALSE)
  comp2 <- igraph::components(g2)$no[1]

  thr <- g1_df[id_sup3[t3]-1, 'w']
  mat_g_end <- x
  mat_g_end[mat_g_end > thr] <- 0
  g_end <- igraph::graph.adjacency(mat_g_end,
                                   mode = "undirected",
                                   weighted = TRUE,
                                   diag = FALSE)
  comp1 <- igraph::components(g_end)$no[1]

  if(sum(comp1, comp2) == 3){

    message(paste("Number of conserved links :",
                  length(igraph::E(g_end))), sep = "")
    message(paste("Maximum weight of the conserved links :",
                  max(igraph::E(g_end)$weight)), sep = "")

  } else {
    stop("Error: there are probably equal link weights.")
  }
  return(g_end)

}

