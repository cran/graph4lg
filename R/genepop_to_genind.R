#' Convert a GENEPOP file into a genind object
#'
#' @description The function converts a text file in the format used by GENEPOP
#' software into a genind object
#'
#' @param path A character string with the path leading to the GENEPOP file
#' in format .txt, or alternatively the name of this file in the working
#' directory.
#' @param n.loci The number of loci in the GENEPOP file (integer or numeric).
#' @param pop_names (optional) Populations' names in the same order
#' as in the GENEPOP file.
#' Vector object (class character) of the same length as the number
#' of populations.
#' Without this parameter, populations are numbered from 1 to the number
#' of populations.
#' @param allele.digit.coding Number indicating whether alleles are coded
#' with 3 (default) or 2 digits.
#' @return An object of type \code{genind}.
#' @export
#' @details This function uses functions from \pkg{pegas} package.
#' GENEPOP file should can include microsatellites loci or SNPs with allele names
#' of length 2 or 3 (noted as 01, 02, 03 or 04 for SNPs).
#' The loci line(s) must not start with a spacing.
#' @seealso For more details about GENEPOP file formatting :
#' \url{https://genepop.curtin.edu.au:443/help_input.html}
#' For the opposite conversion, see \code{\link{genind_to_genepop}}.
#' The output file can be used to compute pairwise FST matrix
#' with \code{\link{mat_pw_fst}}
#' @author P. Savary
#' @examples
#' path_in <- system.file('extdata', 'gpop_simul_10_g100_04_20.txt',
#'                        package = 'graph4lg')
#' file_n <- file.path(tempdir(), "gpop_simul_10_g100_04_20.txt")
#' file.copy(path_in, file_n, overwrite = TRUE)
#' genepop_to_genind(path = file_n, n.loci = 20,
#'                   pop_names = as.character(order(as.character(1:10))))
#' file.remove(file_n)
#' @references \insertRef{raymond1995genepop}{graph4lg}

##################################


genepop_to_genind <- function(path,
                              n.loci,
                              pop_names = NULL,
                              allele.digit.coding = 3){

  # Check whether 'path' is a character string
  if(!is.character(path)){
    stop("You must specify the file path as a character string input 'path'.")
  }

  # Check whether 'n.loci' is an integer
  if(!inherits(n.loci, c("numeric", "integer"))){
    stop("You must specify the number of loci as an integer input 'n.loci'.")
  }

  # Check whether populations' names are not duplicated
  if (any(duplicated(pop_names))){
    message("At least 1 population appears twice in 'pop_names' vector.
            You should regroup individuals from the same population.")
  }

  # Check whether allele.digit.coding is 3 or 2
  if(allele.digit.coding != 3){
    if(allele.digit.coding == 2){
      message("2 allele.digit.coding was used.")
    }
  }

  # Open the text file
  data_gpop <- utils::read.table(file = path,
                                 sep = "\t")
  # Separate the columns of the text file
  sep_col <- suppressWarnings(tidyr::separate(data_gpop, col = 1,
                                              sep = " ",
                                              remove = TRUE,
                                              into = as.character(c(1:(1 + n.loci)))))

  # Sometimes, the loci names are the first lines of the file
  # Alternatively, they can all be in the first line
  # The following if conditions distinguish these different cases

  # Loci on the first lineS
  if (sep_col[(2 + n.loci), 1] == "Pop" || sep_col[(2 + n.loci), 1] == "POP"){
    # Get the loci names which correspond to the first lines of the GENEPOP file
    names.loci <- sep_col[2:(1 + n.loci), 1]
    names.loci <- gsub(names.loci,
                       pattern = "\\.",
                       replacement = "_")

    # Remove the first rows which correspond to the loci names
    sep_col <- sep_col[-(1:(1 + n.loci)), ]

    # Loci on the second line
  } else if (sep_col[3, 1] == "Pop" || sep_col[3, 1] == "POP"){

    # Check whether the GENEPOP is formatted correctly
    if(sep_col[2, 1] == ""){
      stop("The GENEPOP file is badly formatted. Delete spacings before
           loci names")
    }

    # In this case, the loci names are all in the second line
    names.loci <- sep_col[2, (1:n.loci)]
    names.loci <- gsub(names.loci, pattern = "\\.", replacement = "_")
    names.loci[1:(n.loci-1)] <- stringr::str_sub(names.loci[1:(n.loci - 1)],
                                                 1, -2)

    # Remove the two first lines (name of the file, name of the loci)
    sep_col <- sep_col[-(1:2), ]

  }

  # Get the number of populations (corresponds to the number
  # of times a row of sep_col starts with "POP")
  rpop <- which(sep_col[, 1] %in% c("Pop", "POP"))
  # Get the number of individuals by pop (number of rows between 2 "POP")
  n.ind.pop <- c()
  for (i in 1:length(rpop)){
    n.ind.pop[i] <- rpop[i + 1] - rpop[i] - 1
  }
  # Last one : number of rows between the last "POP" and the end
  n.ind.pop[length(rpop)] <- nrow(sep_col) - rpop[length(rpop)]


  # Warning message if there are less than 5 individuals
  if(length(unique(n.ind.pop < 5)) == 2){
    warning("There are populations with less than 5 individuals")
  }

  # Remove the rows with only "POP" at the beginning
  sep_col <- sep_col[-which(sep_col[, 1] %in% c("Pop", "POP")), ]
  # Columns names are "POP" and then the name of every locus
  colnames(sep_col) <- c("POP", as.character(names.loci))
  # Create the column "POP"
  sep_col$POP <- rep(1:length(rpop), times = n.ind.pop)

  # Fill sep_col with the codes of the msats separated by "/"
  all_sep <- lapply(sep_col[, 2:ncol(sep_col)],
                    function(x){
                      paste(stringr::str_sub(x,
                                             1,
                                             allele.digit.coding),
                            "/",
                            stringr::str_sub(x,
                                             allele.digit.coding+1,
                                             2*allele.digit.coding),
                            sep = "")})

  if(ncol(sep_col) == 2){
    sep_col[, 2:ncol(sep_col)] <- unlist(all_sep)
  } else {
    sep_col[, 2:ncol(sep_col)] <- all_sep
  }


  # Replace missing values by NA
  na_rep1 <- lapply(sep_col[, 2:ncol(sep_col)],
                    gsub,
                    pattern = paste0(paste(rep("0", allele.digit.coding),
                                           collapse = ""),
                                     "/",
                                     paste(rep("0", allele.digit.coding),
                                           collapse = "")),
                    replacement = 'NA/NA')

  if(ncol(sep_col) == 2){
    sep_col[, 2:ncol(sep_col)] <- unlist(na_rep1)
  } else {
    sep_col[, 2:ncol(sep_col)] <- na_rep1
  }


  na_rep2 <- lapply(sep_col[, 2:ncol(sep_col)],
                    gsub,
                    pattern = paste0("/",
                                     paste(rep("0", allele.digit.coding),
                                           collapse = "")),
                    replacement = '/NA')

  if(ncol(sep_col) == 2){
    sep_col[, 2:ncol(sep_col)] <- unlist(na_rep2)
  } else {
    sep_col[, 2:ncol(sep_col)] <- na_rep2
  }

  na_rep3 <- lapply(sep_col[, 2:ncol(sep_col)],
                    gsub,
                    pattern = paste0(paste(rep("0", allele.digit.coding),
                                           collapse = ""),
                                     "/"),
                    replacement = 'NA/')


  if(ncol(sep_col) == 2){
    sep_col[, 2:ncol(sep_col)] <- unlist(na_rep3)
  } else {
    sep_col[, 2:ncol(sep_col)] <- na_rep3
  }

  if (is.vector(pop_names)){
    if(length(pop_names) != length(table(sep_col$POP))){
      stop("The length of 'pop_names' is not equal to populations' number
            in the input file.")
    }

    sep_col$POP <- rep(pop_names, times = n.ind.pop)
  } else {
    sep_col$POP <- as.character(sep_col$POP)
  }

  # Convert sep_col to a loci object of pegas and then to a genind object
  # of adegent
  data_genind <- pegas::loci2genind(pegas::as.loci(sep_col, col.pop = 1,
                                                   allele.sep = "/"),
                                    ploidy = 2, na.alleles = "NA")

  return(data_genind)

}








