#' Convert a genind object into a GENEPOP file
#'
#' @description The function converts an object of class \code{genind} into
#' a GENEPOP file.
#' It then allows to use the functionalities of the GENEPOP software and
#' its derived package \pkg{GENEPOP} on R, as well as some functions
#' from other packages (differentiation test, F-stats calculations,
#' HWE test,...).
#' It is designed to be used with diploid microsatellite data with
#' alleles coded with 2 or 3 digits or SNPs genind objects.
#'
#' @param x An object of class \code{genind}
#' from package \pkg{adegenet}.
#' @param output A character string indicating the option used to select what
#' the function will return:
#' \itemize{
#' \item{If \code{output = "data.frame"}(default), then the function will
#' return an object 'x' of class \code{data.frame} ready to be saved as a
#' text file with the following command:
#' \code{write.table(x, file = "file_name.txt", quote=FALSE,
#' row.names=FALSE, col.names=FALSE)}}
#' \item{If \code{output = "path_to_file/file_name.txt"}, then the function
#' will write a text file named 'file_name.txt' in the directory corresponding
#' to 'path_to_file'. Without 'path_to_file', the text file is written in the
#' current working directory. The text file has the format required by GENEPOP
#' software.}
#' }
#' @return An object of type \code{data.frame} if \code{ouput = "data.frame"}.
#' If \code{output} is the path and/or the file name of a text file, then
#' nothing is returned in R environment but a text file is created with the
#' specified file name, either in the current working directory or in the
#' specified folder.
#' @section Warning:
#' \subsection{Confusion}{
#' Do not confound this function with \code{\link[adegenet]{genind2genpop}}
#' from \pkg{adegenet}. The latter converts an object of class \code{genind}
#' into an object of class \code{genpop}, whereas \code{genind_to_genepop}
#' converts an object of class \code{genind} into a text file compatible with
#' GENEPOP software (Rousset, 2008).
#' }
#' \subsection{Allele coding}{
#' This function can handle genetic data with different allele coding: 2 or 3
#' digit coding for microsatellite data or 2 digit coding for SNPs (A,C,T,G
#' become respectively 01, 02, 03, 04).
#' }
#' \subsection{Individuals order}{
#' When individuals in input data are not ordered by populations, individuals
#' from the same population can be separated by individuals from other
#' populations. It can be problematic when calculating then pairwise distance
#' matrices. Therefore, in such a case, individuals are ordered by populations
#' and populations ordered in alphabetic order.
#' }
#' @seealso For more details about GENEPOP file formatting :
#' \url{https://genepop.curtin.edu.au:443/help_input.html}.
#' For the opposite conversion, see \code{\link{genepop_to_genind}}.
#' The output file can be used to compute pairwise FST matrix
#' with \code{\link{mat_pw_fst}}
#' @export
#' @author P. Savary
#' @examples
#' data(data_ex_genind)
#' x <- data_ex_genind
#' df_genepop <- suppressWarnings(genind_to_genepop(x,
#'                                                  output = "data.frame"))
#' @references \insertRef{raymond1995genepop}{graph4lg}

##################################

genind_to_genepop <- function(x, output = "data.frame"){

  # Check whether 'x' is a genind object
  if(!inherits(x, "genind")){
    stop("Input 'x' must be an object of class 'genind'.")
  }

  # Get genetic data from x@tab
  data <- x@tab
  # Get pop_names
  pop_names <- x@pop

  ### Return a message if there is only one population
  if(length(unique(pop_names)) == 1){
    message("There is only one population in your dataset")
    unique_pop <- TRUE
  } else {
    unique_pop <- FALSE
  }

  # If the individuals are not ordered by populations, they are reordered
  # with their populations in alphabetic order.
  if(!all(pop_names == rep(pop_names[-which(duplicated(pop_names))],
                           times = table(pop_names)[unique(pop_names)]))){
    message("Individuals in the input data were not ordered, they have
          been ordered by populations and populations in alphabetic order
          for the conversion.")
    pop_names <- as.character(x@pop)[order(as.character(x@pop))]
    data <- x@tab[order(as.character(x@pop)),]
  } else {

    pop_names <- as.character(x@pop)[order(as.character(x@pop))]
    data <- x@tab[order(as.character(x@pop)),]

  }

  ### Identify the type of markers and modify SNPs data for GENEPOP usage
  if(all(unlist(unique(x@all.names)) %in% c("A", "T", "C", "G"))){
    m_type <- "snp"
    message("Your dataset is treated as a SNP dataset.
            Alleles initially coded A, T, C, G were respectively coded
            01, 02, 03 and 04")

    colnames(data) <- gsub(colnames(data),pattern=".A",replacement=".01")
    colnames(data) <- gsub(colnames(data),pattern=".T",replacement=".02")
    colnames(data) <- gsub(colnames(data),pattern=".C",replacement=".03")
    colnames(data) <- gsub(colnames(data),pattern=".G",replacement=".04")
    #####
  } else {
    m_type <- "msat"
  }

  # Get loci names
  loci_names_l <- x@loc.fac
  # Create an empty data.frame
  loc_all <- data.frame(col = colnames(data))
  # Check whether the names of the locus.allele combination are well formatted
  # with a '.' between the locus and the allele names.
  if(all(stringr::str_count(colnames(data), "\\.") == 1) != TRUE){
    stop("The columns' names of x@tab must be of form 'locus.allele' with only 1
         '.' between locus and allele")
  }
  # Separate the locus and the allele
  loc_all <- tidyr::separate(loc_all, col = 1, sep = "\\.",
                             into = c("locus", "allele"))

  # Avoid problems when allele names have not the same number of characters
  max_all_chr <- max(nchar(loc_all$allele))
  for(i in 1:nrow(loc_all)){
    name_all_chr <- nchar(loc_all[i, "allele"])
    loc_all[i, "allele"] <- ifelse(name_all_chr < max_all_chr,
                                   paste0(rep("0", times = max_all_chr - name_all_chr),
                                          loc_all[i, "allele"]),
                                   loc_all[i, "allele"])
  }


  loci_names <- as.character(loci_names_l[-which(duplicated(loci_names_l))])
  n.loci <- length(loci_names_l[-which(duplicated(loci_names_l))]  )

  # Create the data.frame data_gpop which will contain the genetic data
  # The first column contains the ID of each individual
  data_gpop <- data.frame(id = paste(pop_names, "_", row.names(data), ",",
                                     sep = ""))

  data_gpop[] <- lapply(data_gpop, gsub, pattern = " ", replacement = "")

  # The loop will run over all the loci
  for (i in 1:n.loci){
    # loc is the locus considered
    loc <- loci_names[i]
    # a is an empty vector which will contain the two alleles coding of all
    # individuals over all the loci
    a <- c()
    # The loop will run over all the individuals
    for (j in 1:nrow(data)){
      # Get the column corresponding to the locus
      col_loc <- which(loc_all[, 'locus'] == loc)
      # Get the locus.allele at which individual j is homozygote or heterozygote
      hom <- which(data[j,col_loc] == 2)
      het <- which(data[j,col_loc] == 1)

      # If homozygote, then copy the code of the allele twice
      if(length(hom) != 0){
        a[j] <- paste(loc_all[col_loc[hom], 'allele'],
                      loc_all[col_loc[hom], 'allele'],
                      sep = "")
        # If heterozygote, then copy the code of the two alleles of ind j
        # the lower code number is the first
      } else if (length(het) != 0){

        if (as.character(loc_all[col_loc[het[1]], 'allele']) <
            as.character(loc_all[col_loc[het[2]], 'allele'])){

          a[j] <- paste(loc_all[col_loc[het[1]], 'allele'],
                        loc_all[col_loc[het[2]], 'allele'],
                        sep = "")
        } else {
          a[j] <- paste(loc_all[col_loc[het[2]], 'allele'],
                        loc_all[col_loc[het[1]], 'allele'],
                        sep = "")
        }

      } else {
        # If missing data, set 000000 or 0000
        if(nchar(loc_all[1, 'allele'] == 6)){
          a[j] <- "000000"
        } else {
          a[j] <- "0000"
        }
        ###
      }

    }

    # For each locus, add the next column to data_gpop
    data_gpop <- cbind(data_gpop, a)
  }

  # Then, set the colnames
  colnames(data_gpop) <- c("ID", as.character(loci_names))

  # Convert all the columns as characters
  data_gpop[ , ] <- lapply(data_gpop[ , ], as.character)

  # Add the pop_names
  # Distinguish the case where there is only one pop
  if(unique_pop){

    # Create data_gpop2 whose first line only is "POP"
    data_gpop2 <- rbind(c("Pop",rep("", n.loci)), data_gpop)


  } else {


    # Get the rows' numbers of the end of the populations
    num_end_pop <- c()
    for  (i in (1:(length(pop_names) - 1))){
      if (pop_names[i] != pop_names[i + 1]){
        num_end_pop[i] <- i
      } else {
        num_end_pop[i] <- 0
      }
    }
    num_end_pop <- which(num_end_pop != 0)


    # Create data_gpop2 whose first line only is "POP"
    data_gpop2 <- rbind(c("Pop",rep("", n.loci)), data_gpop[1:num_end_pop[1], ])

    # Then, add the rows of each population, one line with "POP", etc..
    for (i in 1:(length(num_end_pop)-1)){
      data_gpop2 <- rbind(data_gpop2, c("Pop", rep("", n.loci)),
                          data_gpop[(num_end_pop[i]+1):num_end_pop[i+1], ])
    }

    # Add the last population
    data_gpop2 <- rbind(data_gpop2, c("Pop", rep("", n.loci)),
                        data_gpop[(num_end_pop[length(num_end_pop)]+1):nrow(data_gpop), ])

  }

  # The two first lines are the heading and the names of the loci
  if(n.loci > 1){
    data_gpop2 <- rbind(c("Conversion of an object of class 'genind' into a GENEPOP file with the package 'graph4lg'",
                          rep("", n.loci)),
                        c(paste(loci_names[1:(n.loci-1)],",", sep = ""),
                          loci_names[n.loci], ""),
                        data_gpop2)
    colnames(data_gpop2)[2:(length(colnames(data_gpop2))-1)] <-
      paste(colnames(data_gpop2)[2:(length(colnames(data_gpop2))-1)], ",", sep = "")

  } else {
    data_gpop2 <- rbind(c("Conversion of an object of class 'genind' into a GENEPOP file with the package 'graph4lg'",
                          rep("", n.loci)),
                        c(loci_names[n.loci], ""),
                        data_gpop2)
  }

  # Output
  if(output == "data.frame"){
    return(data_gpop2)
  } else if (stringr::str_sub(output,-4,-1) == ".txt") {
    utils::write.table(data_gpop2,file = output, quote=FALSE,
                       row.names=FALSE, col.names=FALSE)
  } else {
    stop("'output' parameter must be 'data.frame' (then the function returns an
object of class data.frame) or a path and file name ending in '.txt'. In the
          latter case, the functions creates a text file in the directory
         specified by the path or in the current working directory when a
         path is not specified.")
  }

}



