#' Compute an index comparing graphs' topologies
#'
#' @description The function computes several indices in order to compare two
#' graphs' topologies. One of the graph has the "true" topology the other is
#' supposed to reproduce. The indices are then a way to assess the reliability
#' of the latter graph.
#' Both graphs must have the same number of nodes, but not necessarily the
#' same number of links. They must also have the same nodes' names and in
#' the same order.
#'
#' @details The indices are calculated from a confusion matrix counting the number
#' of links that are in the "observed" graph ("true") and also in the "predicted"
#' graph (true positives : TP), that are in the "observed" graph but not in the
#' "predicted" graph (false negatives : FN), that are not in the "observed"
#' graph but in the "predicted" graph (false positives : FP) and that are not
#' in the "observed" graph and not in the "predicted" graph neither (true
#' negatives: TN). K is the total number of links in the graphs.
#' K is equal to \eqn{n\times(n-1)} if the graphs are directed and to
#' \eqn{\frac{n\times(n-1)}{2}} if they are not directed, with n the number of nodes.
#' OP = TP + FN, ON = TN + FP, PP = TP + FP and PN = FN + TN.
#'
#' The Matthews Correlation Coefficient (MCC) is computed as follows:
#' \eqn{MCC = \frac{TP\times TN-FP\times FN}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}}
#'
#' The Kappa index is computed as follows:
#' \eqn{Kappa = \frac{K\times (TP + TN) - (ON \times PN) - (OP \times PP)}{K^{2} - (ON \times PN) - (OP \times PP)}}
#'
#' The False Discovery Rate (FDR) is calculated as follows:
#' \eqn{FDR = \frac{FP}{TP+FP}}
#'
#' The Accuracy is calculated as follows:
#' \eqn{Acc = \frac{TP + TN}{K}}
#'
#' The Sensitivity is calculated as follows:
#' \eqn{Sens = \frac{TP}{TP+FN}}
#'
#' The Specificity is calculated as follows:
#' \eqn{Spec = \frac{TN}{TN+FP}}
#'
#' The Precision is calculated as follows:
#' \eqn{Prec = \frac{TP}{TP+FP}}
#'
#'Self loops are not taken into account.
#' @param obs_graph A graph object of class \code{igraph} with n nodes.
#' It is the observed graph that \code{pred_graph} is supposed to approach.
#' @param pred_graph A graph object of class \code{igraph} with n nodes.
#' It is the predicted graph that is supposed to be akin to \code{obs_graph}.
#' @param mode A character string specifying which index to compute in order
#' to compare the topologies of the graphs.
#' \itemize{
#' \item{If 'mode = 'mcc'' (default), the Matthews Correlation Coefficient (MCC) is computed.}
#' \item{If 'mode = 'kappa'', the Kappa index is computed.}
#' \item{If 'mode = 'fdr'', the False Discovery Rate (FDR) is computed.}
#' \item{If 'mode = 'acc'', the Accuracy is computed.}
#' \item{If 'mode = 'sens'', the Sensitivity is computed.}
#' \item{If 'mode = 'spec'', the Specificity is computed.}
#' \item{If 'mode = 'prec'', the Precision is computed.}
#' }
#' @param directed Logical (TRUE or FALSE) specifying whether both graphs
#' are directed or not.
#' @return The value of the index computed
#' @export
#' @author P. Savary
#' @references \insertRef{dyer2004population}{graph4lg}
#' \insertRef{baldi2000assessing}{graph4lg}
#' \insertRef{matthews1975comparison}{graph4lg}
#' @examples
#' data(data_simul_genind)
#' data(pts_pop_simul)
#' mat_dist <- suppressWarnings(graph4lg::mat_geo_dist(data=pts_pop_simul,
#'       ID = "ID",
#'       x = "x",
#'       y = "y"))
#' mat_dist <- mat_dist[order(as.character(row.names(mat_dist))),
#'                       order(as.character(colnames(mat_dist)))]
#' graph_obs <- gen_graph_thr(mat_w = mat_dist, mat_thr = mat_dist,
#'                             thr = 15000, mode = "larger")
#' mat_gen <- mat_gen_dist(x = data_simul_genind, dist = "DPS")
#' graph_pred <- gen_graph_topo(mat_w = mat_gen, mat_topo = mat_dist,
#'                             topo = "gabriel")
#' graph_topo_compar(obs_graph = graph_obs,
#'                   pred_graph = graph_pred,
#'                   mode = "mcc",
#'                   directed = FALSE)


graph_topo_compar <- function(obs_graph, pred_graph, mode = "mcc", directed = FALSE){

  # Check whether obs_graph and pred_graph are graphs
  if(class(obs_graph) != "igraph"){
    stop("'obs_graph' must be a graph object of class 'igraph'.")
  } else if (class(pred_graph) != "igraph"){
    stop("'pred_graph' must be a graph object of class 'igraph'.")
  }

  # Check whether they have the same nodes' number
  if(length(igraph::V(obs_graph)) != length(igraph::V(pred_graph))){
    stop("Both graphs must have the same nodes' number.")
  }


  # Check whether the graphs' nodes have names
  if(is.null(igraph::V(obs_graph)$name)){
    stop("The nodes of 'obs_graph' must have names.")
  } else if(is.null(igraph::V(pred_graph)$name)){
    stop("The nodes of 'pred_graph' must have names.")
  }

  # Check whether the graphs have the same names and in the same order
  if(!all(igraph::V(obs_graph)$name == igraph::V(pred_graph)$name)){
    stop("Both graphs must have the same nodes' names and the nodes ranked in the same order.")
  }


  nb_nodes <- length(igraph::V(obs_graph)$name)

  # If the graph has directed links
  if(directed == TRUE){

    # Number of directed links
    K <- (nb_nodes - 1) * nb_nodes

    # Get the adjacency matrices of both graphs
    adj1 <- igraph::as_adjacency_matrix(obs_graph, names = TRUE , sparse = FALSE)
    adj2 <- igraph::as_adjacency_matrix(pred_graph, names = TRUE , sparse = FALSE)

    diag(adj1) <- diag(adj2) <- 0

    # 1obs 1pred (TP)
    tp <- length(which((adj1 & adj2)))
    # 1 values in both adjacency matrices

    # 0obs 0pred (TN)
    tn <- length(which(((1 - adj1) & (1 - adj2)))) - nb_nodes
    # 0 values become 1 in both adjacency matrices when using 1 - adj1 and 1 - adj2
    # We remove the number of nodes as both diagonal values
    # of 1 - adj1 and 1 - adj2 are ones.

    # 0obs 1pred (FP)
    fp <- length(which(((1 - adj1) & adj2)))
    # We do not remove the number of nodes as only the
    # diagonal values of 1 - adj1 are ones.

    # 1obs 0pred (FN)
    fn <- length(which((adj1 & (1 - adj2))))
    # We do not remove the number of nodes as only the
    # diagonal values of 1 - adj2 are ones.

  # If the graph has non-directed links
  } else if (directed == FALSE){

    K <- (nb_nodes - 1) * nb_nodes/2
    # We consider half the number of links in that case
    # given the graphs are not directed.

    adj1 <- igraph::as_adjacency_matrix(obs_graph, type = "both", names = TRUE , sparse = FALSE)
    adj2 <- igraph::as_adjacency_matrix(pred_graph, type = "both", names = TRUE , sparse = FALSE)

    diag(adj1) <- diag(adj2) <- 0

    # 1obs 1pred (TP)
    tp <- length(which((adj1 & adj2)))/2
    # We divide by 2 in each case as the adjacency matrices are symmetric

    # 0obs 0pred (TN)
    tn <- length(which(((1 - adj1) & (1 - adj2))))/2 - (nb_nodes / 2)
    # We divide by 2 and then we remove half the number of nodes as both diagonal values
    # of 1 - adj1 and 1 - adj2 are ones.

    # 0obs 1pred (FP)
    fp <- length(which(((1 - adj1) & adj2)))/2

    # 1obs 0pred (FN)
    fn <- length(which((adj1 & (1 - adj2))))/2

  } else {
    stop("'directed' must be TRUE or FALSE.")
  }

  # Predicted positive = true positive + false positive
  pp <- tp + fp
  # Predicted negative = false negative + true negative
  pn <- fn + tn

  # Observed positive = true positive + false negative
  op <- tp + fn
  # Observed negative = true negative + false positive
  on <- tn + fp

  # We compute products because otherwise we reach the limit
  # of the numbers R can handle in calculations.
  pp_op <- pp * op
  pn_on <- pn * on

  mcc <- (tp * tn - fp * fn) / ( sqrt(pp_op) * sqrt(pn_on) )
  kappa <- ( K * (tp + tn) - (on * pn) - (op * pp)) / (K^2 - (on * pn) - (op * pp))
  fdr <- fp / (tp + fp)
  acc <- (tp + tn) / K
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  prec <- tp / (tp + fp)

  if(mode == "mcc"){
    message(paste("Matthews Correlation Coefficient : ", mcc, sep = ""))
    return(mcc)
  } else if(mode == "kappa"){
    message(paste("Kappa Index : ", kappa, sep = ""))
    return(kappa)
  } else if(mode == "fdr"){
    message(paste("False Discovery Rate : ", fdr, sep = ""))
    return(fdr)
  } else if(mode == "acc"){
    message(paste("Accuracy : ", acc, sep = ""))
    return(acc)
  } else if(mode == "sens"){
    message(paste("Sensitivity : ", sens, sep = ""))
    return(sens)
  } else if(mode == "spec"){
    message(paste("Specificity : ", spec, sep = ""))
    return(spec)
  } else if(mode == "prec"){
    message(paste("Precision : ", prec, sep = ""))
    return(prec)
  }

}










