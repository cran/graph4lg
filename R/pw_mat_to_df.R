#' Convert a pairwise matrix into an edge-list data.frame
#'
#' @description The function converts a pairwise matrix into an edge-list
#' data.frame
#'
#' @param pw_mat A pairwise matrix which can be:\itemize{
#' \item{An object of class \code{matrix}. It must have
#' the same row names and column names. If values represent distances,
#' diagonal elements should be equal to 0.}
#' \item{An object of class \code{dist}. In that, its column numbers are
#' used to create IDs in the resulting data.frame.}
#' }
#' @return An object of class \code{data.frame}
#' @export
#' @author P. Savary
#' @examples
#' data(data_tuto)
#' pw_mat <- data_tuto[[1]]
#' df <- pw_mat_to_df(pw_mat)


pw_mat_to_df <- function(pw_mat){

  if(!inherits(pw_mat, c("matrix", "dist"))){
    stop("'pw_mat' must be a matrix or dist object")
  } else if (inherits(pw_mat, "matrix")){

    if(is.null(row.names(pw_mat))){
      stop("'pw_mat' must have row names")
    } else if(is.null(colnames(pw_mat))){
      stop("'pw_mat' must have column names")
    }

    if(!isSymmetric(pw_mat)){
      stop("'pw_mat' must be a symmetric matrix")
    }

    if(!all(row.names(pw_mat) == colnames(pw_mat))){
      stop("Row names and column names of 'pw_mat' must be equal")
    }

    if(!all(diag(pw_mat) == 0)){
      warning("All diagonal elements of pw_mat were not equal to 0.")
    }

  } else if (inherits(pw_mat, "dist")){
    pw_mat <- as.matrix(pw_mat)
  }

  vec_name <- row.names(pw_mat)

  # Df with all the potential pairs
  df_pw <- data.frame(expand.grid(vec_name, vec_name))
  df_pw[, 1:2] <- lapply(df_pw[, 1:2], function(x){as.character(x)})

  # Remove the diagonal elements which corresponds to 0 values
  df_pw <- df_pw[-which(df_pw$Var1 == df_pw$Var2), ]

  # Remove values where Var1 is higher than Var2
  # (avoid duplicates)
  # If numeric conversion returns NA, then compare character
  if(suppressWarnings(any(is.na(as.numeric(df_pw$Var1))))){

    df_pw <- df_pw[-which(as.character(df_pw$Var1) > as.character(df_pw$Var2)), ]

    # If numeric conversion does not return NA, then compare numeric orders
  } else {

    df_pw <- df_pw[-which(as.numeric(df_pw$Var1) > as.numeric(df_pw$Var2)), ]

  }


  # Create a unique ID
  df_pw$id_link <- paste0(df_pw$Var1, "_", df_pw$Var2)

  # Set column names
  colnames(df_pw) <- c("id_1", "id_2", "id_link")

  # Set row names
  row.names(df_pw) <- df_pw$id_link

  # Create a column 'value' with NA at first
  df_pw$value <- NA

  # Fill 'value' column
  for(row in 1:nrow(df_pw)){

    id_1 <- df_pw[row, 'id_1']
    id_2 <- df_pw[row, 'id_2']

    df_pw[row, 'value'] <- pw_mat[id_1, id_2]

  }

  return(df_pw)
}
