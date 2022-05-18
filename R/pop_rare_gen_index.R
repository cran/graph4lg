#' Compute population-level rarefied genetic indices with ADZE software
#'
#' @description The function computes population-level rarefied genetic indices
#' from an object of class \code{genind} with the ADZE software.
#'
#' @param x An object of class \code{genind}
#' from package \pkg{adegenet}.
#' @param max_g (optional default = NULL) The maximum standardized sample size
#' used by ADZE software (MAX_G) in ADZE manual. It is equal to twice the
#' minimum number of individuals considered for the rarefaction analysis. By
#' default, it is equal to twice the number of individuals in the smallest
#' population. Ohterwise, it must be either a numeric or integer value.
#' @param pop_names (optional) A character vector indicating population names.
#' It is of the same length as the number of populations. Without this
#' argument, populations are given the names they have initially in the
#' 'genind' object (which is sometimes only a number). The order of the
#' population names must match with their order in the 'genind' object.
#' The function does not reorder them. Users must be careful.
#' @param OS A character string indicating whether you use a Linux ('linux')
#' or Windows ('win') operating system.
#' @return An object of class \code{data.frame} whose rows
#' correspond to populations and columns to population attributes
#' (ID, size, genetic indices). By default, the first column corresponds to
#' the population names (ID). The order of the columns depends on the
#' vector 'indices'.
#' @export
#' @keywords internal
#' @author P. Savary

pop_rare_gen_index <- function(x,
                               max_g = NULL,
                               pop_names = NULL,
                               OS = "linux"){


  ### Check for ADZE and param file
  data_dir <- rappdirs::user_data_dir()

  # Check for the directory
  if(!dir.exists(paths = paste0(data_dir, "/graph4lg_jar"))){
    dir.create(path = paste0(data_dir, "/graph4lg_jar"))
    stop(paste0("The directory ", paste0(data_dir, "/graph4lg_jar"),
                " has been created. You must now copy the adze-1.0 software ",
                "in this directory, as well as the 'paramfile_adze' file."))
  }

  # Check for the software
  if(!any(c("adze-1.0",
            "adze-1.0.exe") %in% list.files(paste0(data_dir,
                                                   "/graph4lg_jar")))){
    stop(paste0("You must copy the adze-1.0 software in the directory ",
                paste0(data_dir, "/graph4lg_jar"),
                " as well as the 'paramfile_adze' file."))
  }

  # Check for the param file
  if(!("paramfile_adze" %in% list.files(paste0(data_dir,
                                                   "/graph4lg_jar")))){
    stop(paste0("You must copy the paramfile_adze file in the directory ",
                paste0(data_dir, "/graph4lg_jar")))
  }

  # Check whether 'x' is a 'genind' object
  if(!inherits(x, "genind")){
    stop("x must be a 'genind' object.")
  }

  # Check whether 'pop_names' is either NULL or numeric
  if(!is.null(max_g)){
    if(!inherits(max_g, c("numeric", "integer"))){
      stop("max_g must be either NULL or a numeric or integer value.")
    }
  }


  # Check whether 'x' is a 'genind' object
  if(!inherits(OS, "character")){
    stop("'OS' must be a character string.")
  } else if(!(OS %in% c("linux", "win"))){
    stop("'OS' must be equal to either 'linux' or 'win'.")
  }

  # Get the names of the populations in the 'genind' object 'x'
  pop_x <- levels(x@pop)
  # Get the number of populations
  nb_pop <- length(pop_x)

  if(nb_pop == 1){
    message("There is only one population in your dataset")
  }

  # If 'pop_names' is NULL (default), then population names are the same as in
  # the 'genind' object
  # ELSE they are derived from 'pop_names' vector when provided.
  if(is.null(pop_names)){
    pop_names <- pop_x
  } else {
    pop_names <- pop_names
  }

  # Check whether 'pop_names' is of class "character"
  if(!inherits(pop_names, "character")){
    stop("'pop_names' must be a character vector.")
  }

  # Check whether 'pop_names' has as many elements as there are
  # populations in 'x'
  if(length(pop_names) != nb_pop){
    stop("'pop_names' must have as many as elements as there are
         populations in 'x'.")
  }

  # Get nb ind. and set data_lines
  nb_ind <- nrow(x@tab)
  data_lines <- 2 * nb_ind

  # nb_loci
  nb_loci <- length(unique(x@loc.fac))

  # max_g
  if(is.null(max_g)){
    max_g <- 2 * min(table(x@pop))
  }

  # Path to adze and to param_file
  path_to_adze <- paste0(rappdirs::user_data_dir(), "/graph4lg_jar/adze-1.0")

  param_file <- paste0(rappdirs::user_data_dir(),
                       "/graph4lg_jar/paramfile_adze")

  # data_file
  data_file <- tempfile(fileext = ".txt")
  # output_file
  r_out <- tempfile()
  p_out <- tempfile()
  c_out <- tempfile()
  k_range <- nb_pop - 1

  # generate data file in STRUCTURE format
  genind_to_structure(x = x,
                      output = data_file)


  if(OS == "linux"){
    # Run the command for Linux
    system(paste(path_to_adze,
                 param_file,
                 "-g", max_g,
                 "-d", data_lines,
                 "-l", nb_loci,
                 "-nr", 1,
                 "-nc", 2,
                 "-s", 2,
                 "-f", data_file,
                 "-r", r_out,
                 "-p", p_out,
                 "-k", k_range,
                 "-o", c_out))
  } else {
    # Run the command for Windows
    shell(paste(path_to_adze,
                 param_file,
                 "-g", max_g,
                 "-d", data_lines,
                 "-l", nb_loci,
                 "-nr", 1,
                 "-nc", 2,
                 "-s", 2,
                 "-f", data_file,
                 "-r", r_out,
                 "-p", p_out,
                 "-k", k_range,
                 "-o", c_out))
  }



  ## Analyze the results
  ### Open the files

  # Allelic richness
  all_rich <- utils::read.table(file = r_out)
  colnames(all_rich) <- c("pop", "g", "locus", "mean_ar", "var_ar", "std_ar")
  all_rich <- all_rich[which(all_rich$g == max_g), ]

  # Private allelic richness
  priv_rich <- utils::read.table(file = p_out)
  colnames(priv_rich) <- c("pop", "g", "locus", "mean_priv", "var_priv", "std_priv")
  priv_rich <- priv_rich[which(priv_rich$g == max_g), ]

  # Missing allelic richness
  miss_rich <- utils::read.table(file = paste0(c_out, "_", k_range))
  colnames(miss_rich) <- c(paste0("pop", as.character(1:k_range)), "g", "locus", "mean_miss", "var_miss", "std_miss")
  miss_rich <- miss_rich[which(miss_rich$g == max_g), ]
  miss_rich$pop <- as.character(c(miss_rich[1, 1:k_range],
                                  miss_rich[nrow(miss_rich), k_range]))
  miss_rich <- miss_rich[, c("pop", "g", "locus", "mean_miss", "var_miss", "std_miss")]

  df_res <- merge(all_rich, priv_rich, by = "pop")
  df_res <- merge(df_res, miss_rich, by = "pop")

  df_res <- df_res[, c("pop",
                         "mean_ar", "var_ar", "std_ar",
                         "mean_priv", "var_priv", "std_priv",
                         "mean_miss", "var_miss", "std_miss")]
  return(df_res)
}






