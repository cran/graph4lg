#' Convert a file from \pkg{gstudio} or \pkg{popgraph} into a genind object
#'
#' @description The function converts a file formatted to use \pkg{gstudio} or
#' \pkg{popgraph} package into a genind object (\pkg{adegenet} package)
#'
#' @param x An object of class \code{data.frame} with loci columns in
#' format \code{locus} (defined in package \pkg{gstudio}) with as many rows as
#' individuals and as many columns in format \code{locus} as there are loci and
#' additional columns
#' @param pop_col A character string indicating the name of the column with
#' populations' names in \code{x}
#' @param ind_col (optional) A character string indicating the name of the
#' column with individuals' ID in \code{x}
#' @return An object of class \code{genind}.
#' @export
#' @details This function uses functions from \pkg{pegas} package.
#' It can handle genetic data where alleles codings do not have same length,
#' (99:101, for example).
#' If the names of the loci include '.' characters, they will
#' be replaced by '_'.
#' @author P. Savary
#' @examples
#' data("data_ex_gstud")
#' x <- data_ex_gstud
#' pop_col <- "POP"
#' ind_col <- "ID"
#' data_genind <- gstud_to_genind(x, pop_col, ind_col)


gstud_to_genind <- function(x,
                            pop_col,
                            ind_col = NULL){

  # Check whether 'x' is a 'data.frame'
  if(!inherits(x, "data.frame")){
    stop("'x' must be an object of class 'data.frame' with columns
          corresponding to loci of class 'locus'.")
  }

  # Check whether 'pop_col' is a character string
  if(!is.character(pop_col)){
    stop("You must specify the column name with the population names
         as a character string input 'pop_col'.")
  }

  # Check whether at least one column of 'x' is of class 'locus'
  if(!any(unlist(lapply(x, function(x){inherits(x, "locus")})))){
    stop("Columns corresponding to loci must be of class 'locus'.")
  }

  # Get the number of the useful columns for this operation
  # (pop names and locus)
  col <- c()
  for (i in 1:ncol(x)){
    if(inherits(x[, i], "locus")){
      col[i] <- i
    } else {
      col[i] <- 0
    }
  }
  col[ncol(x)+1]<- c(which(colnames(x) == pop_col))
  col <- col[-which(col == 0)]
  col <- col[c(length(col), 1:(length(col)-1))]

  # Reorder individuals if necessary
  if(!all(x[, pop_col] == x[order(x[, pop_col]), pop_col])){
    message("Individuals in 'x' were not ordered, they have
            been ordered by populations and populations ordered in alphabetic
            order for the conversion.")

    x <- x[order(x[, pop_col]), ]
  }

  # Create data, a data.frame with columns from x with locus and pop names
  data <- x[, col]

  # The first column is the population names, following ones are loci
  # Loci columns become characters
  data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], as.character)

  # Allele separators become '/'
  data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)],
                                 gsub,
                                 pattern = ':',
                                 replacement = '/')

  # Missing values become 'NA/NA'
  data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)],
                                function(x) ifelse(x == "", "NA/NA", x))

  # Row names become 'ind_col' column when specified
  if(is.vector(ind_col)){
    row.names(data) <- x[,ind_col]
  }

  # '.' are replaced in columns' names by '_'
  colnames(data) <- gsub(colnames(data),
                         pattern = "\\.",
                         replacement = "_")

  # data is converted into a loci object and then into a genind object
  data_genind <- pegas::loci2genind(pegas::as.loci(data, col.pop = 1),
                                    ploidy = 2,
                                    na.alleles = "NA")

  # Display a warning message if there are populations with less
  # than 5 individuals
  if(length(unique(as.numeric(table(x[, pop_col])) < 5)) == 2){
    warning("There are populations with less than 5 individuals")
  }

  return(data_genind)
}
