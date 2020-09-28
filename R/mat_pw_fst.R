#' Compute a pairwise FST matrix between populations
#'
#' @description The function computes the pairwise FST matrix between
#' populations from an object of class \code{genind} or
#' directly from a GENEPOP file.
#'
#' @param x An object of class \code{genind}, or the
#' character string indicating the path of the GENEPOP file.
#' @param pop_names (optional) A vector of class \code{character} of the same
#' length as the number of populations (row and column number in the returned
#' matrix). It contains the name of the populations.
#' @return A pairwise \code{matrix} of FST with as many rows and columns as
#' there are populations in the input data.
#' @details The formula used is that of Weir et Cockerham (1984).
#' This functions uses directly the function \code{diffCalc}
#' from \pkg{diveRsity}.
#' See \url{https://genepop.curtin.edu.au:443/help_input.html} for details on the
#' GENEPOP file format and see Raymond (1995) for detail about GENEPOP software.
#' @section Warnings:
#' The order of populations matters :
#' \itemize{
#' \item If \code{x} is an object of class \code{genind}, individuals are
#' re-ordered by populations and populations are ordered in alphabetic order.
#' \item If \code{x} is the path to a GENEPOP file, population order
#' in \code{pop_names} must be the same as in the GENEPOP file.
#' }
#' Negative values are converted into 0
#' @keywords internal
#' @export
#' @author P. Savary
#' @references \insertRef{weir1984estimating}{graph4lg}
#' \insertRef{raymond1995genepop}{graph4lg}
#' @examples
#' data("data_ex_genind")
#' mat_d_j <- mat_pw_d_j(data_ex_genind)
#' path_in <- system.file('extdata', 'gpop_simul_10_g100_04_20.txt',
#'                        package = 'graph4lg')
#' file_n <- file.path(tempdir(), "gpop_simul_10_g100_04_20.txt")
#' file.copy(path_in, file_n, overwrite = TRUE)
#' mat_pw_fst(x = file_n, pop_names = as.character(order(as.character(1:10))))
#' file.remove(file_n)

##################################

mat_pw_fst <- function(x, pop_names = NULL){

  # If 'x' is a 'genind' object
  if (inherits(x, "genind")){
    # Create a temporary text file name
    tmp <- tempfile(fileext = ".txt")
    # Convert the 'genind' object into a GENEPOP formatted text file
    genind_to_genepop(x, output = tmp)

    # Compute the distance matrix between populations using diffCalc from
    # diveRsity package
    mat_fst <- diveRsity::diffCalc(infile = tmp,
                                   outfile = NULL,
                                   fst = TRUE,
                                   pairwise = TRUE)

    # Get the FST matrix
    mat_fst <- mat_fst$pairwise$Fst
    # Add rows and columns names
    pop_names <- x@pop[order(as.character(x@pop))]
    pop_names <- as.character(pop_names[-which(duplicated(pop_names))])
    row.names(mat_fst) <- colnames(mat_fst) <- pop_names

  # If 'x' is the path to a GENEPOP formatted text file
  } else if (class(x) == "character"){

    # Compute the distance matrix between populations using diffCalc from
    # diveRsity package, directly from the GENEPOP file
    mat_fst<-diveRsity::diffCalc(infile = x,
                                 outfile = NULL,
                                 fst = TRUE,
                                 pairwise = TRUE)

    # Get the FST matrix
    mat_fst<-mat_fst$pairwise$Fst

    # Add rows and columns names after having checked that pop_names
    # has as many elements as there are rows and columns in the distance matrix
    if (is.vector(pop_names)){
      if(length(pop_names) != nrow(mat_fst)){
        stop(paste("'pop_names' must have ",
                   length(row.names(mat_fst)), " elements, it has ",
                   length(pop_names),sep=""))
      }

      # Also check that no population is duplicated in 'pop_names'
      if(any(duplicated(pop_names))){
        stop("At least one population appears twice in 'pop_names'.")
      }
    }
    row.names(mat_fst) <- colnames(mat_fst) <- pop_names

  } else {
    stop("Input value x must be either a 'genind' object or a
         character string.")
  }

  # Create a function to make a matrix symmetric
  makeSymm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }

  # Create the full distance matrix with the function makeSymm
  if(all(row.names(mat_fst) == colnames(mat_fst))){
    mat_fst <- makeSymm(mat_fst)
  } else {
    stop("An error occured while making the matrix symmetric.")
  }
  # Diagonal elements are 0
  diag(mat_fst) <- rep(0, length(diag(mat_fst)))
  # Negative values are replaced by 0
  mat_fst[mat_fst < 0] <- 0

  return(mat_fst)

}




