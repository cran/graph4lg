#' Convert a file in STRUCTURE format into a genind object
#'
#' @description The function converts a text file in STRUCTURE format into
#' a genind object to use in R
#'
#' @param path A character string indicating the path to the STRUCTURE file in
#' format .txt, or alternatively the name of the file in the working directory.
#' The STRUCTURE file must only have :
#' \itemize{
#' \item A first column with the IDs of the individuals (can be a simple number)
#' \item A second column with the IDs of the populations (can be a simple number)
#' \item Some loci columns : as many columns as loci in the data
#' }
#' The row for loci's names is optional but recommended.
#' Each individual is displayed on 2 rows.
#'
#' @param loci_names A character vector with the names of the loci if not
#' specified in the file's first row. This argument is mandatory if the
#' STRUCTURE file does not include the names of the loci in the first row.
#' In other cases, the names of the loci is extracted from the file's first row
#' @param pop_names (optional) A character vector indicating the populations'
#' names in the same order as in the STRUCTURE file. It is of the same length
#' as the number of populations. Without this argument, populations are
#' numbered from 1 to the total number of individuals.
#' @param ind_names (optional) A character vector indicating the individuals'
#' names in the same order as in the STRUCTURE file. It is of the same length
#' as the number of individuals. Without this argument, individuals are
#' numbered from 1 to the total number of individuals.
#' @return An object of type \code{genind}.
#' @details The columns' order of the resulting object can be different from
#' that of objects returned by \code{\link{gstud_to_genind}}
#' and \code{\link{genepop_to_genind}}, depending on alleles' and loci's coding
#' This function uses functions from \pkg{pegas} package.
#' For details about STRUCTURE file format : \href{http://www.ccg.unam.mx/~vinuesa/tlem09/docs/structure_doc.pdf}{STRUCTURE user manual}
#' @export
#' @author P. Savary
#' @examples
#' data("data_pc_genind")
#' loci_names <- levels(data_pc_genind@loc.fac)
#' pop_names <- levels(data_pc_genind@pop)
#' ind_names <- row.names(data_pc_genind@tab)
#' path_in <- system.file('extdata', 'tab_gstud_structure.txt',
#'                        package = 'graph4lg')
#' file_n <- file.path(tempdir(), "tab_gstud_structure.txt")
#' file.copy(path_in, file_n, overwrite = TRUE)
#' str <- structure_to_genind(path = file_n, loci_names = loci_names,
#'                            pop_names = pop_names, ind_names = ind_names)
#' file.remove(file_n)

##################################

structure_to_genind <- function(path,
                                pop_names = NULL,
                                loci_names = NULL,
                                ind_names = NULL){

  # Check whether 'path' is a character string
  if(!is.character(path)){
    stop("You must specify the file's path as a character string input 'path'.")
  }

  # Check whether individuals are regrouped in the data
  if (any(duplicated(pop_names))){
    warning("At least 1 population appears twice in 'pop_names' vector.
            You should regroup individuals from the same population.")
  }

  # Read the text file
  data_str <- utils::read.table(file = path,
                                sep = "\t")

  # If the first line does not correspond to the first individual
  if(stringr::str_sub(data_str[1, ], 1, 2) != "1 "){
    # Get the loci's names from the first line
    loci_names <- unlist(strsplit(as.character(data_str[1, ]), split = " "))

    # Replace '.' by '_' in the loci's names
    loci_names <- gsub(loci_names,
                       pattern = "\\.",
                       replacement = "_")

    # Re-open the data skipping the first line and with another separator
    data_str <- utils::read.table(file = path, sep = "", skip = 1)

  # If the first line corresponds to the first individual
  } else {

    # Check whether loci's names are given
    if (!is.vector(loci_names)){
      stop("You must specify the names of the loci in a character vector
           because input file does not include them.")
    }
    # Re-open the data with another separator
    data_str <- utils::read.table(file = path, sep = "")

    # Check whether there are as many loci's names in 'loci_names' than
    # in the data
    if (length(loci_names)!= (ncol(data_str)-2)){
      stop("The length of 'loci_names' is not equal to the number of loci
            in the input file.")
    }

    # Replace '.' by '_' in the loci's names
    loci_names <- gsub(loci_names,
                       pattern = "\\.",
                       replacement = "_")
  }

  # Give columns' names to the data
  col_names <- unlist(c("ind", "pop", loci_names))
  colnames(data_str) <- col_names

  # Check whether there are as many populatons in 'pop_names' than in the data
  if (is.vector(pop_names)){
    ind.pop <- as.numeric(table(data_str$pop))

    if(length(pop_names)!=length(ind.pop)){
      stop("The length of 'pop_names' is not equal to populations' number
            in the input file.")
    }

    # Give populations' names to the individuals
    data_str$pop <- rep(pop_names, times = ind.pop)
  }

  # Check whether there are as many individuals in 'ind_names' than in the data
  if (is.vector(ind_names)){

    if(length(ind_names) != (nrow(data_str)/2)){
      stop("The length of 'ind_names' is not equal to individuals' number
           in the input file.")
    }

    # Give individuals' names to the individuals
    data_str$ind <- rep(ind_names, each = 2)
  }

  # Check whether there are at least 5 individuals in every population
  if(length(unique(as.numeric(table(data_str[, 'pop'])) < 11)) == 2){
    warning("There are populations with less than 5 individuals")
  }

  # Split the data in two: even rows and uneven rows
  # (every individual genotype is displayed on two consecutive rows)
  data_uneven <- data_str[seq(1, nrow(data_str), 2), ]
  data_even <- data_str[-seq(1, nrow(data_str), 2), ]

  # Merge the data by separating them with "/"
  data_union <- data_even
  for (i in 3:ncol(data_union)){
    data_union[, i] <- paste(data_even[, i], "/",
                             data_uneven[, i],
                             sep = "")
  }

  # Replace '-9/-9' by 'NA/NA'
  data_union[, 3:ncol(data_union)] <- lapply(data_union[, 3:ncol(data_union)],
                                             gsub,
                                             pattern = '-9/-9',
                                             replacement = 'NA/NA')

  # Add rows' names
  row.names(data_union) <- data_union$ind
  # Remove the 'ind' column
  data_union <- data_union[, -which(colnames(data_union) == "ind")]

  # Convert data_union into a 'loci' object and then into a 'genind' object
  data_genind <- pegas::loci2genind(pegas::as.loci(data_union,
                                                   col.pop = 1,
                                                   allele.sep = "/"),
                                    ploidy = 2,
                                    na.alleles = "NA")

  return(data_genind)

}



