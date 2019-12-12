#' Compute a pairwise genetic distance matrix between populations
#' using Bowcock et al. (1994) formula
#'
#' @description The function computes the pairwise DPS, a genetic distance
#' based on the proportion of shared alleles.
#'
#' @details The formula used is inspired from MSA software :
#' \deqn{D_{PS}=1-\frac{\sum_{d}^{D}\sum_{k}^{K}\min (f_{a_{kd}i},f_{a_{kd}j})}{D} }
#' such as \eqn{a_{kd}} is the allele \eqn{k} at locus \eqn{d}
#' \eqn{D} is the total number of loci
#' \eqn{K} is the allele number at each locus
#' \eqn{\gamma_{a_{kd^{ij}}}=0} if individuals \eqn{i} and \eqn{j}
#' do not share allele \eqn{a_{kd}}
#' \eqn{\gamma_{a_{kd^{ij}}}=1} if one of individuals \eqn{i} and \eqn{j}
#' has a copy of \eqn{a_{kd}}
#' \eqn{\gamma_{a_{kd^{ij}}}=2} if both individuals have 2 copies
#' of \eqn{a_{kd}} (homozygotes)
#' \eqn{f_{a_{kd}i}} is allele \eqn{a_{kd}} frequency in
#' individual \eqn{i} (0, 0.5 or 1).
#' More information in :
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/7510853}{Bowcock et al., 1994}
#' and \href{http://i122server.vu-wien.ac.at/MSA/info.html/Distances.html#DPS}{MSA manual}.
#' This function uses functions from \pkg{adegenet} package
#' Note that in the paper of Bowcock et al. (1994), the denominator is 2D.
#' But, in MSA software manual, the denominator is D.
#'
#' @param x An object of class \code{genind}
#' @return A pairwise matrix of genetic distances between populations
#' @export
#' @author P. Savary
#' @examples
#' data("data_simul_genind")
#' dist_bowcock <- mat_pw_dps(data_simul_genind)
#' @references \insertRef{bowcock1994high}{graph4lg}

################################################################################

mat_pw_dps <- function(x) {

  # Check whether 'x' is a 'genind' object
  if(!inherits(x, "genind")){
    stop("Input 'x' must be an object of class 'genind'.")
  }

  #sink("aux")
  # f is a large matrix with the frequency of each allele in each population
  f <- adegenet::makefreq(adegenet::genind2genpop(x), missing = 0)
  #sink()

  # n.loci = number of loci
  n.loci <- length(unique(x@loc.fac))

  # K : number of populations
  K <- nrow(f)
  # Empty symmetric matrix of dimension K*K
  ret <- matrix(0, K, K)
  rownames(ret) <- colnames(ret) <- rownames(f)
  # Fill the matrix by applying the formula
  for (i in 1:K) {
    for (j in 1:i) {
      if (i != j) {
        ret[i, j] <- ret[j, i] <- 1 - (sum(apply(f[c(i, j), ], 2, min))/n.loci)
      }
    }
  }
  return(ret)
}

