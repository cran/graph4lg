#' Create an independence graph of genetic differentiation
#' from genetic data of class genind
#'
#' @description The function allows to create genetic graphs from genetic data
#' by applying the conditional independence principle. Populations whose allelic
#' frequencies covary significantly once the covariance with the other
#' populations has been taken into account are linked on the graphs.
#'
#' @details The function allows to vary many parameters such as the genetic
#' distance used, the formula used to compute the covariance, the statistical
#' tolerance threshold, the p-values adjustment, among others.
#'
#' @param x An object of class \code{genind} that contains the multilocus
#' genotype (format 'locus') of the individuals as well as their population
#' and their geographical coordinates.
#' @param dist A character string indicating the method used to compute the
#' multilocus genetic distance between populations
#' \itemize{
#' \item{If 'dist = 'basic'' (default), then the multilocus genetic distance is
#' computed using a Euclidean genetic distance formula (Excoffier et al., 1992)}
#' \item{If 'dist = 'weight'', then the multilocus genetic distance is computed
#' as in Fortuna et al. (2009). It is a Euclidean genetic distance giving more
#' weight to rare alleles}
#' \item{If 'dist = 'PG'', then the multilocus genetic distance is computed as
#' in popgraph::popgraph function, following several steps of PCA and SVD
#' (Dyer et Nason, 2004).}
#' \item{If 'dist = 'PCA'', then the genetic distance is computed following a
#' PCA of the matrix of allelic frequencies by population. It is a Euclidean
#' genetic distance between populations in the multidimensional space defined
#' by all the independent principal components.}
#' }
#' @param cov A character string indicating the formula used to compute the
#' covariance matrix from the distance matrix
#' \itemize{
#' \item{If 'cov = 'sq'' (default), then the covariance matrix is calculated
#' from the matrix of squared distances as in Everitt et Hothorn (2011)}
#' \item{If 'cov = 'dist'', then the covariance matrix is calculated from the
#' matrix of distances as in Dyer et Nason (2004) and popgraph function}
#' }
#' @param pcor A character string indicating the way the partial correlation
#' matrix is computed from the covariance matrix.
#' \itemize{
#' \item{If 'pcor = 'magwene'', the steps followed are the same as in
#' Magwene (2001) and in popgraph::popgraph function. It is the recommended
#' option as it meets mathematical requirements.}
#' \item{If 'pcor = 'other'', the steps followed are the same as used
#' by Fortuna et al. (2009). They are not consistent with the approach
#' of Magwene (2001).}
#' }
#' @param test A character string indicating the method used to test the
#' significance of the partial correlation coefficients.
#' \itemize{
#' \item{If 'test = 'EED'' (default), then the Edge Exclusion Deviance
#' criterion is used (Whittaker, 2009). Although other methods exist, this is
#' the most common and thus the only one implemented here.}
#' }
#' @param alpha A numeric value corresponding to the statistical tolerance
#' threshold used to test the difference from 0 of the partial correlation
#' coefficients. By default, 'alpha=0.05'.
#' @param adj A character string indicating the way of adjusting p-values to
#' assess the significance of the p-values
#' \itemize{
#' \item{If 'adj = 'none'' (default), there is no p-value adjustment correction}
#' \item{If 'adj = 'holm'', p-values are adjusted using the sequential
#' Bonferroni correction (Holm, 1979)}
#' \item{If 'adj = 'bonferroni'', p-values are adjusted using the classic
#' Bonferroni correction}
#' \item{If 'adj = 'BH'', p-values are adjusted using Benjamini et Hochberg
#' (1995) correction controlling false discovery rate}
#' }
#' @param output A character string indicating the matrices included in
#' the output list.
#' \itemize{
#' \item{If 'output = 'all'' (default), then D (distance matrix),
#' C (covariance matrix), Rho (partial correlation matrix),
#' M (graph incidence matrix) and S (strength matrix) are included}
#' \item{If 'output = 'dist_graph'', then the distance matrix D is returned
#' only with the values corresponding to the graph edges}
#' \item{If 'output = 'str_graph'', then the strength values matrix S is
#' returned only with the values corresponding to the graph edges}
#' \item{If 'output = 'inc'', then the binary adjacency matrix M is returned}
#' \item{If 'output = 'igraph'', then a graph of class \code{igraph}
#' is returned}
#' }
#' @return A \code{list} of objects of class \code{matrix}, an object of
#' class \code{matrix} or a graph object of class \code{igraph}
#' @export
#' @author P. Savary
#' @examples
#' data(data_pc_genind)
#' dist_graph_test <- gen_graph_indep(x = data_pc_genind, dist = "basic",
#'                              cov = "sq", pcor = "magwene",
#'                              alpha = 0.05, test = "EED",
#'                              adj = "none", output = "igraph")
#' @references \insertRef{dyer2004population}{graph4lg}
#' \insertRef{benjamini1995controlling}{graph4lg}
#' \insertRef{bowcock1994high}{graph4lg}
#' \insertRef{everitt2011introduction}{graph4lg}
#' \insertRef{excoffier1992analysis}{graph4lg}
#' \insertRef{fortuna2009networks}{graph4lg}
#' \insertRef{holm1979simple}{graph4lg}
#' \insertRef{magwene2001new}{graph4lg}
#' \insertRef{wermuth1977algorithm}{graph4lg}
#' \insertRef{whittaker2009graphical}{graph4lg}


gen_graph_indep <- function(x, dist = "basic", cov = "sq",
                            pcor = "magwene", alpha = 0.05,
                            test = "EED", adj = "none",
                            output = "igraph") {

  ############################################

  ##### Data preparation and allelic frequencies computation


  ############################################

  # tab is a matrix with as many rows as there are individuals and
  # as many columns as there are loci.
  tab <- x@tab
  # group_init is the vector with the populations of every individual
  group_init <- x@pop
  group_init <- factor(as.character(group_init))
  # n is the number of individuals
  n <- nrow(tab)
  # rows is the number of populations
  rows <- length(levels(x@pop))

  # Creation of a vector with the number of different alleles at each locus.
  # It has the length of the total number of alleles in the data
  df_freq_all <- data.frame(Locus = x@loc.fac,
                            Allele = as.character(unlist(x@all.names)))
  df_freq_all$id <- paste(df_freq_all$Locus, ".",
                          df_freq_all$Allele, sep = "")

  # x_w is a large matrix with the frequency of each allele in each individual
  #sink("aux")
  x_w <- adegenet::makefreq(x, missing = 0)
  #sink(NULL)

  # freq_all is a vector with the mean allelic frequencies over all individuals
  freq_all <- apply(x_w, 2, mean)

  # Attribute its frequency to each locus-allele combination
  # after checking the order
  if(unique(df_freq_all$id == names(freq_all)) == TRUE){
    df_freq_all$Frequency <- freq_all
  } else {
    stop("Locus-Allele combinations are not in the same order.")
  }

  # Number of loci
  n.loci <- length(unique(x@loc.fac))
  # Number of alleles
  n.all <- length(x@loc.fac)

  # Creation of a matrix of size number of populations * number of alleles
  # with allelic frequencies at each locus in each population
  #sink("aux")
  A <- adegenet::makefreq(adegenet::genind2genpop(x))
  #sink(NULL)
  # Number of columns = number of alleles
  columns <- ncol(A)
  # Groups = names of the populations
  groups <- row.names(A)

  # D is an empty symmetric matrix with as many rows as there are populations.
  # Rows' names and columns' names are populations' names
  D <- matrix(0, nrow = rows, ncol = rows)
  row.names(D) <- colnames(D) <- groups

  #####################################

  ### Distance computation #####

  #####################################

  # Fortuna et al. 2009 method : allele weighting
  if (dist == "weight") {
    # Creation of a vector with the number of alleles at each locus
    # repeated as many times as there are alleles (for each allele, it
    # gives the number of alleles at the same locus)
    vec.n.all <- as.vector(rep(x@loc.n.all, times = x@loc.n.all))

    # Creation of a vector with the frequency of the alleles at each locus
    # over all the data
    Freq <- df_freq_all$Frequency

    # The loop will run along all the populations and all the alleles
    for (i in 1:rows) {
      for (j in 1:rows) {
        # If non-diagonal element
        if (i != j) {
          # Initialization of suma_loci --> 0
          suma_loci <- 0
          # The loop will run along all the columns of A
          for (m in 1:columns) {
            # k is the number of alleles at the locus of the allele considered
            k <- vec.n.all[m]
            # p_k is the frequency of the allele considered
            p_k <- Freq[m]
            # w is the weight given to the allele in the computation of
            # the element of the distance between i and j corresponding
            # to this allele
            w <- 1/(k * p_k)
            # We add the element corresponding to this allele to suma_loci
            # This way, we compute a Euclidean distance with a particular
            # weighting of each element of the sum.
            suma_loci <- suma_loci + (w * ((A[i, m] - A[j, m])^2))
          }
          # d_ij is the distance between i and j : square root of (1/2 x sum of
          # squared weighted elements)
          d_ij <- sqrt(0.5 * suma_loci)
          D[i, j] <- d_ij
          # If diagonal element --> 0
        } else {
          D[i, j] <- 0
        }
      }
    }

    # PCA-derived distance
  } else if (dist == "PCA") {

    # PCA of the total table of allele frequencies by individuals
    pcfit <- stats::prcomp(x_w, retx = TRUE)
    # We project every individual along the principal components
    # We keep only the number of principal components equal to the total
    # number of alleles minus the total number of loci
    x_w2 <- pcfit$x[, 1:(n.all - n.loci)]
    # We update columns, the number of "genetic variables" considered
    columns <- ncol(x_w2)
    # We have ncol(x_w2) = n.all - n.loci

    # Mean coordinates per population and principal components
    x_w2_mean <- tapply(x_w2, list(rep(group_init, columns), col(x_w2)), mean)

    # Double-centering by columns and by rows so that covariance calculation
    # is correct
    x_w2_mean <- scale(x_w2_mean, scale = FALSE, center = TRUE)
    x_w2_mean <- t(scale(t(x_w2_mean), scale=FALSE, center = TRUE))

    # Euclidean distance between population in the new space
    D <- as.matrix(stats::dist(x_w2_mean, diag = TRUE, upper = TRUE,
                               method = "euclidean"))

    # Basic method : no weighting, no transformation (loci with many alleles
    # have more weight)
  } else if (dist == "basic") {

    # No PCA
    # No double-centering here because row sums are constant
    # Centering by rows does not affect distance (be careful if missing data)
    # Centering by columns never affects distance by definition
    # rowSums(A)
    D <- as.matrix(stats::dist(A, diag = TRUE, upper = TRUE,
                               method = "euclidean"))

    # Popgraph method : many transformations to reduce dimensionality
    # and colinearity before distance calculation
  } else if (dist == "PG") {
    # This part of the code is directly inspired from popgraph package
    tol <- 1e-04
    critVal <- stats::qchisq(1 - alpha, 1)

    # throw warning of some groups are small
    t <- table(group_init)
    if (any(t < 4)) {
      popnames <- paste(names(which(t < 4)), collapse = ", ")
      message(paste("You have strata (", popnames, ") that have fewer
                    than 4 individuals. This analysis needs to have a good
                    estimate of within stratum variance."))
    }

    # PCA of all the data rotate all data and keep stuff that is not invariant.

    pcfit <- stats::prcomp(x_w, retx = TRUE)
    mv <- pcfit$x[, pcfit$sdev > tol]

    P <- ncol(mv)
    Pop.priors <- as.numeric(table(group_init)/n)
    pop.means <- tapply(mv, list(rep(group_init, P), col(mv)), mean)
    sigma.w <- sqrt(diag(stats::var(mv - pop.means[group_init, ])))

    if (any(sigma.w < tol))
      message(paste("Dropped rotated genetic variable '",
                    paste((as.character(1:P)[sigma.w < tol]), collapse = ","),
                    "' from the analysis due to constancy within groups",
                    sep = ""))

    scaling <- diag(1/sigma.w, , P)
    fac <- 1/(n - rows)
    X <- sqrt(fac) * (mv - pop.means[group_init, ]) %*% scaling
    X.s <- svd(X, nu = 0)
    rank <- sum(X.s$d > tol)
    # if(rank < P) warning( paste( (P-rank), ' variables are collinear and
    # being dropped from the discriminant rotation.', sep=''))
    scaling <- scaling %*% X.s$v[, 1:rank] %*% diag(1/X.s$d[1:rank], , rank)
    mu <- colSums(Pop.priors %*% pop.means)
    X <- sqrt((n * Pop.priors)/(rows - 1)) * scale(pop.means,
                                                   center = mu,
                                                   scale = FALSE) %*% scaling
    X.s <- svd(X, nu = 0)
    rank <- sum(X.s$d > tol * X.s$d[1L])
    scaling <- scaling %*% X.s$v[, 1L:rank]

    means <- colMeans(pop.means)
    LDValues <- scale(mv, center = means, scale = FALSE) %*% scaling

    allLD <- tapply(LDValues, list(rep(group_init, rank), col(LDValues)), mean)
    #allLD <- centroid_distance(LDValues, groups)

    allSD <- rep(NA, rows)
    names(allSD) <- groups

    for( i in 1:rows ){
      xp <- LDValues[ group_init == groups[i],]
      allSD[i] <- sum(diag(cov(xp)))
    }

    #allSD <- centroid_variance(LDValues, group_init)

    D <- as.matrix(stats::dist(allLD, diag = TRUE, upper = TRUE,
                               method = "euclidean"))

    # There are differences between two ways of computing euclidean distances
    # only at the 14th digits when rounding.

    #for (i in seq(1, rows)) for (j in seq(i, rows)) {
    # if (i != j) {
    #  p1 <- unlist(allLD[i, ])
    #  p2 <- unlist(allLD[j, ])
    #  D[i, j] <- D[j, i] <- sqrt(sum((p1 - p2)^2))
    #}
    #}

    ## Be careful to the names order in the resulting matrix.
    ## It should different from the order in the matrices obtained
    ## with other methods


  #} else if (dist == "DPS") {

   # D <- mat_gen_dist(x, dist = "DPS")
   # # We make the distance Euclidean with the Cailliez transformation
   # D <- as.matrix(ade4::cailliez(stats::as.dist(D)))

  #} else if (dist == "FST_lin") {

   # D <- mat_gen_dist(x, dist = "FST_lin", null_val = TRUE)
   # # We make the distance Euclidean
   # D <- as.matrix(ade4::cailliez(stats::as.dist(D)))

  #} else if (dist == "FST") {

   # D <- mat_gen_dist(x, dist = "FST", null_val = TRUE)
   # # We make the distance Euclidean
   # D <- as.matrix(ade4::cailliez(stats::as.dist(D)))

  } else {
    stop("You must specify a correct 'dist' option.")
  }

  #########################################################

  ##### Covariance matrix

  #########################################################

  # Creation of an empty covariance matrix
  C <- matrix(0, nrow = rows, ncol = rows)
  row.names(C) <- colnames(C) <- row.names(D)

  # Covariance matrix computed from matrix of squared distances
  if (cov == "sq") {
    dist_sq <- D * D
    mean_rows <- rowMeans(dist_sq)
    mean_columns <- colMeans(dist_sq)
    mean_total <- mean(dist_sq)

    for (i in 1:rows) {
      for (j in 1:rows) {
        d_ij <- dist_sq[i, j]
        di <- mean_rows[i]
        dj <- mean_columns[j]
        C[i, j] <- -0.5 * (d_ij - di - dj + mean_total)
      }
    }

    # Covariance matrix computed from matrix of distances
  } else if (cov == "dist") {
    mean_rows <- rowMeans(D)
    mean_columns <- colMeans(D)
    mean_total <- mean(D)

    for (i in 1:rows) {
      for (j in 1:rows) {
        d_ij <- D[i, j]
        di <- mean_rows[i]
        dj <- mean_columns[j]
        C[i, j] <- -0.5 * (d_ij - di - dj + mean_total)
      }
    }

  } else {
    stop("You must specify a correct 'cov' option.")
  }

  #####################################################################

  ##### Partial correlation calculation

  #####################################################################

  # Computation steps directly inspired from Magwene et al. 2001
  if (pcor == "magwene") {

    # Correlation matrix
    R <- as.matrix(stats::cov2cor(C))

    # Correlation matrix inversion MASS::ginv(x) to calculate Moore-Penrose
    # generalized inverse matrix
    W <- MASS::ginv(R)
    row.names(W) <- colnames(W) <- row.names(R)

    # Partial correlation matrix
    Rho <- as.matrix(stats::cov2cor(W))
    Rho <- -Rho
    diag(Rho) <- 1

    #Rho <- matrix(0, nrow = rows, ncol = rows)
    #row.names(Rho) <- colnames(Rho) <- row.names(W)
    #for (i in 1:rows) {
    # for (j in 1:rows) {
    #  if (i != j) {
    #   Rho[i, j] <- -W[i, j]/((W[i, i] * W[j, j])^0.5)
    #  } else {
    #   Rho[i, j] <- 1
    # }
    #}
    #}
    # Differences after 13 digits when rounding

    # Computation step as used by Fortuna et al. (2009)
  } else if (pcor == "other") {

    # Precision matrix
    P <- MASS::ginv(C)
    row.names(P) <- colnames(P) <- row.names(C)

    # Partial correlation matrix
    Rho <- as.matrix(stats::cov2cor(P))
    Rho <- -Rho
    diag(Rho) <- 1

    #Rho <- matrix(0, nrow = rows, ncol = rows)
    #row.names(Rho) <- colnames(Rho) <- row.names(P)
    #for (i in 1:rows) {
    # for (j in 1:rows) {
    #  if (i == j) {
    #    Rho[i, j] <- P[i, j]/((P[i, i] * P[j, j])^0.5)
    #  } else {
    #    Rho[i, j] <- -P[i, j]/((P[i, i] * P[j, j])^0.5)
    #  }
    #}
    #}
    # Differences after 13 digits when rounding

  } else {
    stop("You must specify a correct 'pcor' option.")
  }

  ##############################################################

  #### Test of the significance of the partial correlation coefficient values.

  ##############################################################

  if (test == "EED") {

    # Without considering differently the negative partial correlation
    # coefficient values
    if (dist %in% c("PG", "weight")) {

      # Deviance matrix
      Dev <- matrix(0, nrow = rows, ncol = rows)
      row.names(Dev) <- colnames(Dev) <- row.names(Rho)

      Rho2 <- 1 - Rho^2
      Rho2[Rho2 < 0] <- 0
      Dev <- -n * log(Rho2)

      # Considering the negative partial correlation
      # coefficient values as zeros
    } else {

      # Deviance matrix
      Dev <- -n * log(1 - (Rho^2))
      Dev[Rho < 0] <- 0
      Dev[Rho == 1] <- 10
      diag(Dev) <- 0

    }


    ###### Test with or without adjustment option #####

    if (adj == "none") {

      # Significance threshold at the level alpha and 1 degree of freedom
      # for the chi-square distribution
      critVal <- stats::qchisq(1 - alpha, 1)

      # Incidence matrix M based on deviance values significance
      M <- matrix(0, nrow = rows, ncol = rows)
      row.names(M) <- colnames(M) <- row.names(Dev)
      for (i in 1:rows) {
        for (j in 1:rows) {
          # If Dev < critVal, non-significant partial correlation
          # and no link
          if (Dev[i, j] < critVal) {
            M[i, j] <- 0
          # If Dev >= critVal, significant partial correlation --> link
          } else {
            M[i, j] <- 1
          }
        }
      }

    } else if (adj %in% c("holm","bonferroni","BH")) {

      # pval associated to each EED value
      pval1 <- 1 - stats::pchisq(Dev, df = 1)

      # Data frame with the p-values (NA initially), the names of the pop
      # correponding to every partial correlation value.
      pval1_df <- data.frame(pop1 = rep(row.names(pval1), rows),
                             pop2 = rep(row.names(pval1), each = rows),
                             p = rep(NA, rows * rows))
      # Attribute the correspong p-value to every row
      for (i in 1:rows) {
        for (j in 1:i) {
          pop1 <- row.names(pval1)[i]
          pop2 <- colnames(pval1)[j]
          p <- pval1[i, j]
          pval1_df[which(pval1_df$pop1 == pop1 & pval1_df$pop2 == pop2), 3] <- p
        }
      }
      # Remove the diagonal elements of the matrix in the data.frame
      pval1_df <- pval1_df[-which(pval1_df$pop1 == pval1_df$pop2), ]
      pval1_df <- pval1_df[-which(is.na(pval1_df$p)), ]
      # Computation of the adjusted p-value
      pval1_df$p_cor <- stats::p.adjust(as.vector(pval1_df$p), method = adj)

      # Incidence matrix based on deviance values significance
      M <- matrix(0, nrow = rows, ncol = rows)
      row.names(M) <- colnames(M) <- row.names(pval1)
      for (i in 1:rows) {
        for (j in 1:i) {
          pop1 <- row.names(M)[i]
          pop2 <- colnames(M)[j]
          if(pop1 != pop2){
            M[i, j] <- M[j, i] <- ifelse(pval1_df[which(pval1_df$pop1 == pop1 &
                                                        pval1_df$pop2 == pop2),
                                                  "p_cor"] < alpha, 1, 0)
          } else {
            M[i, j] <- M[j, i] <- 0
          }
        }
      }


    } else {
      stop("You must specify a correct 'adj' option.")
    }


  } else {
    stop("You must specify a correct 'test' option.")
  }

  ########################################################################

  ###### Output definition

  ########################################################################



  if (output == "all") {
    res_list <- list("D", D, "C", C, "Rho", Rho, "M", M)
    return(res_list)
  } else if (output == "dist_graph") {
    D[M == 0] <- 0
    return(D)
  } else if (output == "inc") {
    return(M)
  } else if (output == "igraph") {
    D[M == 0] <- 0
    graph <- igraph::graph.adjacency(D, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)
    class(graph) <- c("igraph")
    return(graph)
  } else {
    stop("You must specify a correct 'output' option.")
  }

}


# Not included finally

# \item{If 'dist = 'DPS'', then the genetic distance used is equal to
# 1 - the proportion of shared alleles (Bowcock, 1994)}
# \item{If 'dist = 'FST'', then the genetic distance used is the pairwise
# FST (Weir et Cockerham, 1984)}
# \item{If 'dist = 'FST_lin'', then the genetic distance used is the linearised
# pairwise FST (Weir et Cockerham, 1984)}

