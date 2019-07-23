#' Compute a pairwise matrix of genetic distances between populations
#'
#' @description The function computes a pairwise matrix of genetic distances
#' between populations and allows to implement several formula.
#'
#' @param x An object of class \code{genind} that contains the multilocus
#' genotypes (format 'locus') of the individuals as well as their populations.
#' @param dist A character string indicating the method used to compute the
#' multilocus genetic distance between populations
#' \itemize{
#' \item{If 'dist = 'basic'' (default), then the multilocus genetic distance is
#' computed using a formula of Euclidean genetic distance (Excoffier et al., 1992)}
#' \item{If 'dist = 'weight'', then the multilocus genetic distance is computed
#' as in Fortuna et al. (2009). It is a Euclidean genetic distance giving more
#' weight to rare alleles}
#' \item{If 'dist = 'PG'', then the multilocus genetic distance is computed as
#' in popgraph::popgraph function, following several steps of PCA and SVD
#' (Dyer et Nason, 2004).}
#' \item{If 'dist = 'DPS'', then the genetic distance used is equal to
#' 1 - the proportion of shared alleles (Bowcock, 1994)}
#' \item{If 'dist = 'FST'', then the genetic distance used is the pairwise
#' FST (Weir et Cockerham, 1984)}
#' \item{If 'dist = 'FST_lin'', then the genetic distance used is the linearised
#' pairwise FST (Weir et Cockerham, 1984)(FST_lin = FST/(1-FST))}
#' \item{If 'dist = 'PCA'', then the genetic distance is computed following a
#' PCA of the matrix of allelic frequencies by population. It is a Euclidean genetic
#' distance between populations in the multidimensional space defined by all
#' the independent principal components.}
#' \item{If 'dist = 'GST'', then the genetic distance used is the G'ST (Hedrick, 2005)}
#' \item{If 'dist = 'D'', then the genetic distance used is Jost's D (Jost, 2008)}
#' }
#' @param null_val (optional) Logical. Should negative and null FST, FST_lin,
#' GST or D values be replaced by half the minimum positive value?
#' This option allows to compute Gabriel graphs from these "distances".
#' Default is null_val = FALSE.
#' This option only works if 'dist = 'FST'' or 'FST_lin' or 'GST' or 'D'
#' @return An object of class \code{matrix}
#' @export
#' @details Negative values are converted into 0.
#' Euclidean genetic distance \eqn{d_{ij}} between population i and j is computed as
#' follows: \deqn{d_{ij}^{2} = \sum_{k=1}^{n} (x_{ki} - x_{kj})^{2} } where
#' \eqn{x_{ki}} is the allelic frequency of allele k in population i and n is
#' the total number of alleles. Note that when 'dist = 'weight'', the formula
#' becomes \deqn{d_{ij}^{2} = \sum_{k=1}^{n} (1/(K*p_{k}))(x_{ki} - x_{kj})^{2}}
#' where K is the number of alleles at the locus of the allele k and \eqn{p_{k}}
#' is the frequency of the allele k in all populations.
#' Note that when 'dist = 'PCA'', n is the number of conserved independent
#' principal components and \eqn{x_{ki}} is the value taken by the principal
#' component k in population i.
#' @author P. Savary
#' @references \insertRef{bowcock1994high}{graph4lg}
#' \insertRef{excoffier1992analysis}{graph4lg}
#' \insertRef{dyer2004population}{graph4lg}
#' \insertRef{fortuna2009networks}{graph4lg}
#' \insertRef{weir1984estimating}{graph4lg}
#' \insertRef{hedrick2005standardized}{graph4lg}
#' \insertRef{jost2008gst}{graph4lg}
#' @examples
#' data(data_simul_genind)
#' x <- data_simul_genind
#' D <- mat_gen_dist(x = x, dist = "basic")

mat_gen_dist <- function(x,
                         dist = "basic",
                         null_val = FALSE){

  # Check whether 'x' is of class 'genind'
  if(class(x) != "genind"){
    stop("x must be an object of type genind")
  }

  # Check whether 'null_vall' is specified only with the right type of distance
  if(!any(c(dist == "FST", dist == "FST_lin", dist == "GST", dist =="D"))){
    if(null_val){
      stop("This option only works if dist = 'FST', 'FST_lin', 'GST' or 'D'.")
    }
  }

  # Get the allelic data tab
  tab <- x@tab
  # Get the groups names of each individual
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

  # Fortuna et al. 2009 method : allele weighting
  if (dist == "weight"){
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
            # w is the weight given to the allele in the computation of the
            # element of the distance between i and j corresponding to this allele
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
    pcfit <- stats::prcomp(x_w, retx = T)
    # We project every individual along the principal components
    # We keep only the number of principal components equal to the total
    # number of alleles minus the total number of loci
    x_w2 <- pcfit$x[, 1:(n.all - n.loci)]
    # We update columns, the number of "genetic variables" considered
    columns <- ncol(x_w2)
    # We have ncol(x_w2) = n.all - n.loci

    # Mean coordinates per population and principal components
    x_w2_mean <- tapply(x_w2, list(rep(group_init, columns), col(x_w2)), mean)

    # Double-centering by columns and by rows so that covariance calculation is correct
    x_w2_mean <- scale(x_w2_mean, scale = FALSE, center = TRUE)
    x_w2_mean <- t(scale(t(x_w2_mean), scale=FALSE, center = TRUE))

    # Euclidean distance between population in the new space
    D <- as.matrix(stats::dist(x_w2_mean, diag = TRUE, upper = TRUE, method = "euclidean"))

    # Basic method : no weighting, no transformation (loci with many alleles have more weight)
  } else if (dist == "basic") {

    # No PCA
    # No double-centering here because row sums are constant
    # Centering by rows does not affect distance (be careful if missing data)
    # Centering by columns never affects distance by definition
    # rowSums(A)
    D <- as.matrix(stats::dist(A, diag = TRUE, upper = TRUE, method = "euclidean"))


  # Popgraph method : many transformations to reduce dimensionality
  # and colinearity before distance calculation
  } else if (dist == "PG"){
    # This part of the code is directly inspired from popgraph package
    tol <- 1e-4
    pcfit <- stats::prcomp(x_w, retx = T)
    mv <- pcfit$x[, pcfit$sdev > tol]

    P <- ncol(mv)
    Pop.priors <- as.numeric(table(group_init)/n)
    pop.means <- tapply(mv, list(rep(group_init, P), col(mv)), mean)
    sigma.w <- sqrt(diag(stats::var(mv - pop.means[group_init, ])))

    if (any(sigma.w < tol))
      warning(paste("Dropped rotated genetic variable '", paste((as.character(1:P)[sigma.w < tol]), collapse = ","), "' from the analysis due to constancy within groups",
                    sep = ""))

    scaling <- diag(1/sigma.w, , P)
    fac <- 1/(n - rows)
    X <- sqrt(fac) * (mv - pop.means[group_init, ]) %*% scaling
    X.s <- svd(X, nu = 0)
    rank <- sum(X.s$d > tol)
    # if(rank < P) warning( paste( (P-rank), ' variables are collinear and being dropped from the discriminant rotation.', sep=''))
    scaling <- scaling %*% X.s$v[, 1:rank] %*% diag(1/X.s$d[1:rank], , rank)
    mu <- colSums(Pop.priors %*% pop.means)
    X <- sqrt((n * Pop.priors)/(rows - 1)) * scale(pop.means, center = mu, scale = FALSE) %*% scaling
    X.s <- svd(X, nu = 0)
    rank <- sum(X.s$d > tol * X.s$d[1L])
    scaling <- scaling %*% X.s$v[, 1L:rank]

    means <- colMeans(pop.means)
    LDValues <- scale(mv, center = means, scale = FALSE) %*% scaling

    allLD <- tapply(LDValues, list(rep(group_init, rank), col(LDValues)), mean)

    D <- as.matrix(stats::dist(allLD, diag = TRUE, upper = TRUE, method = "euclidean"))

  # DPS : based on the proportion of shared alleles
  } else if (dist == "DPS"){

    # Use of the function 'mat_pw_dps()'
    D <- mat_pw_dps(x)

  # FST
  } else if (dist == "FST"){

    x@tab <- x@tab[order(as.character(x@pop)), ]
    x@pop <- x@pop[order(as.character(x@pop))]

    # Use of the function 'mat_pw_fst()'
    D <- mat_pw_fst(x = x)
    # Negative values are converted into 0
    D[D < 0] <- 0

  } else if (dist == "GST"){

    x@tab <- x@tab[order(as.character(x@pop)), ]
    x@pop <- x@pop[order(as.character(x@pop))]

    # Use of the function 'mat_pw_gst()'
    D <- mat_pw_gst(x = x)
    # Negative values are converted into 0
    D[D < 0] <- 0

  } else if (dist == "D"){

    x@tab <- x@tab[order(as.character(x@pop)), ]
    x@pop <- x@pop[order(as.character(x@pop))]

    # Use of the function 'mat_pw_d_j()'
    D <- mat_pw_d_j(x = x)
    # Negative values are converted into 0
    D[D < 0] <- 0

  } else if (dist == "FST_lin"){

    x@tab <- x@tab[order(as.character(x@pop)), ]
    x@pop <- x@pop[order(as.character(x@pop))]

    # Use of the function 'mat_pw_fst()'
    mat <- mat_pw_fst(x = x)
    # Negative values are converted into 0
    mat[mat < 0] <- 0

    # FST_lin = FST/(1 - FST)
    D <- mat
    for (i in 1:rows){
      for (j in 1:rows){
        D[i, j] <- mat[i, j]/( 1 - mat[i, j])
      }
    }

  } else  {
    stop("You must specify a correct 'dist' option.")
  }

  # Null_vall option implementation
  if(any(c(dist, dist, dist, dist) == c("FST","FST_lin","GST","D"))){
    if(null_val){
      mat2 <- D
      diag(mat2) <- 1
      mat2[mat2 == 0] <- 4
      val_rep <- 0.5 * min(mat2)
      D[which(mat2 == 4)] <- val_rep
    }
  }

  return(D)

}


