#' Compute population-level genetic indices
#'
#' @description The function computes population-level genetic indices from an
#' object of class \code{genind}.
#'
#' @param x An object of class \code{genind}
#' from package \pkg{adegenet}.
#' @param pop_names (optional) A character vector indicating populations' names.
#' It is of the same length as the number of populations. Without this
#' argument, populations are given the names they have initially in the
#' 'genind' object (which is sometimes only a number). The order of the
#' populations' names must match with their order in the 'genind' object.
#' The function does not reorder them. Users must be careful.
#' @param indices (optional) A character vector indicating the population-level
#' indices to compute. These indices can be:
#' \itemize{
#' \item{Mean allelic richness by locus by population (\code{indices = c("A", ...)})}
#' \item{Mean expected heterozygosity by locus by population (\code{indices = c("He",...)})}
#' \item{Mean observed heterozygosity by locus by population (\code{indices = c("Ho",...)})}
#' \item{Number of individuals by population (\code{indices = c("Nb_ind", ...)})}
#' \item{Total allelic richness by population (\code{indices = c("A_tot",...)})}
#' }
#' By default, \code{indices = c("Nb_ind", "A", "He", "Ho")}.
#' @return An object of class \code{data.frame} whose rows correspond to populations
#' and columns to populations' attributes (ID, size, genetic indices). By default,
#' the first column corresponds to the populations' names (ID). The order of the
#' columns depends on the vector 'indices'.
#' @export
#' @author P. Savary
#' @examples
#' data(data_simul_genind)
#' x <- data_simul_genind
#' pop_names <- levels(x@pop)
#' df_pop_indices <- pop_gen_index(x = x,
#'                    pop_names = pop_names,
#'                    indices = c("Nb_ind", "A", "He", "Ho"))

pop_gen_index <- function(x,
                          pop_names = NULL,
                          indices = c("Nb_ind", "A", "He", "Ho")){

  # Check whether 'x' is a 'genind' object
  if(class(x) != "genind"){
    stop("x must be a 'genind' object.")
  }

  # Get the names of the populations in the 'genind' object 'x'
  pop_x <- levels(x@pop)
  # Get the number of populations
  nb_pop <- length(pop_x)

  # Check whether 'pop_names' is of class "character"
  if(class(pop_names) != "character"){
    stop("'pop_names' must be a character vector.")
  }

  # Check whether 'pop_names' has as many elements as there are populations in 'x'
  if(length(pop_names) != nb_pop){
    stop("'pop_names' must have as many as elements as there are populations in 'x'.")
  }

  # If 'pop_names' is NULL (default), then populations' names are the same as in
  # the 'genind' object
  # ELSE they are derived from 'pop_names' vector when provided.
  if(is.null(pop_names)){
    pop_names <- pop_x
  } else {
    pop_names <- pop_names
  }

  # Check whether indices consist of valid character strings
  if(!all(indices %in% c("A", "He", "Ho", "Nb_ind", "A_tot"))){
    stop("indices must consist of valid character strings
         ('A', 'He' , 'Ho', 'Nb_ind', 'A_tot')")
  }

  # Summary of the data
  summ_x <- adegenet::summary(x)

  # Number of individuals by population
  Nb_ind <- summ_x[[2]]

  # Total allelic richness by population
  A_tot <- summ_x[[4]]
  # Number of loci
  n_loci <- length(summ_x[[3]])
  # Mean allelic richness by population and by locus
  A <- A_tot/n_loci

  # We split the data into the different populations
  x.pop <- adegenet::seppop(x)

  # We create a list of summaries of the data by population
  summ_x_pop <- lapply(x.pop, summary)

  # We get the population-level heterozygosities (observed and expected)
  Ho <- rep(NA, length(summ_x_pop))
  He <- rep(NA, length(summ_x_pop))


  for (j in 1:length(summ_x_pop)){
    summ_x_pop_j <- summ_x_pop[[j]]
    Ho[j] <- summ_x_pop_j['Hobs']
    He[j] <- summ_x_pop_j['Hexp']
  }

  # We gather the data into a data.frame
  df_pop_indices <- data.frame(ID = pop_names,
                               Nb_ind = Nb_ind,
                               A = A,
                               Ho = Ho,
                               He = He,
                               A_tot = A_tot)

  # We keep only the columns indicated by the vector 'indices'
  df_pop_indices <- df_pop_indices[, c(1, which(colnames(df_pop_indices) %in% indices))]

  # The function returns the data.frame 'df_pop_indices'
  return(df_pop_indices)

}







