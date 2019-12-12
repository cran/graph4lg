#' Convert a loci object into a genind object
#'
#' @description This function is exactly the same as \code{loci2genind}
#' from \pkg{pegas} package
#'
#' @param x An object of class \code{loci} to convert
#' @param ploidy An integer indicating the ploidy level
#' (by default, 'ploidy = 2')
#' @param na.alleles  A character vector indicating the coding of the alleles
#' to be treated as missing data (by default, 'na.alleles = c("NA")')
#' @return An object of class \code{genind}
#' @export
#' @author P. Savary
#' @examples
#' data("data_pc_loci")
#' genind <- loci_to_genind(data_pc_loci, ploidy = 2, na.alleles = "NA")

loci_to_genind <- function(x,
                           ploidy = 2,
                           na.alleles = c("NA")){

  # Check whether 'x' is a 'loci' object
  if(!inherits(x, "loci")){
    stop("Input 'x' must be a 'loci' object.")
  }

  # Reorder individuals if necessary
  if(!all(x$population == x$population[order(x$population)])){
    message("Individuals in 'x' were not ordered, they have
            been ordered by populations and populations ordered in alphabetic
            order for the convertion.")

    x <- x[order(x$population), ]
  }

  # Convert 'x' into a 'genind' object
  data_genind <- pegas::loci2genind(x,
                                    ploidy,
                                    na.alleles)

  return(data_genind)

}





