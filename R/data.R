################################################################################

##### LARGE SIMULATED DATA ###########################################################

################################################################################

#' data_simul_genind genetic dataset
#'
#' Genetic dataset from genetic simulation on CDPOP
#' 1500 individuals, 50 populations
#' 20 microsatellite loci (3 digits coding)
#' 50 generations simulated
#'
#' @format An object of type 'genind'
#' @references \insertRef{landguth2010cdpop}{graph4lg}
#' @details The simulation was made with CDPOP during 50 generations.
#' Dispersal was possible between the 50 populations. Its probability depended
#' on the cost distance between populations, calculated on a simulated
#' resistance surface (raster). Mutations were not possible. There
#' were initially 600 alleles in total (many disappeared because of drift).
#' Population stayed constant
#' with a sex-ratio of 1. Generations did not overlap.
#' This simulation includes a part of stochasticity and these data result
#' from only 1 simulation run.
#' @examples
#' data("data_simul_genind")
#' length(unique(data_simul_genind@pop))
"data_simul_genind"

############################################################################

#' pts_pop_simul : details on simulated populations
#'
#' Simulation dataset
#' 50 populations located on a simulated landscape
#'
#' @format An object of class 'data.frame' with the following columns :
#' \describe{
#' \item{ID}{Population ID of the 50 populations}
#' \item{x}{Site longitude (RGF93)}
#' \item{y}{Site latitude (RGF93)}
#' }
#' @references \insertRef{landguth2010cdpop}{graph4lg}
#' There are as many rows as there are sampled populations.
#' @examples
#' data("pts_pop_simul")
#' str(pts_pop_simul)
"pts_pop_simul"



############################################################################

#' data_tuto : data used to generate the vignette
#'
#' Data used to generate the vignette
#'
#' @format Several outputs or inputs to show how the package works in a list
#' \describe{
#' \item{mat_dps}{Genetic distance matrix example}
#' \item{mat_pg}{Second genetic distance matrix example}
#' \item{graph_ci}{Genetic independence graph example}
#' \item{dmc}{Output of the function 'dist_max_corr'}
#' \item{land_graph}{Landscape graph example}
#' \item{mat_ld}{Landscape distance matrix example}
#' }
#' @examples
#' data("data_tuto")
#' mat_dps <- data_tuto[[1]]
#' str(mat_dps)
"data_tuto"


################################################################################

##### SMALL SIMULATED DATA ###########################################################

################################################################################

#' data_ex_genind genetic dataset
#'
#' Genetic dataset from genetic simulation on CDPOP
#' 200 individuals, 10 populations
#' 20 microsatellite loci (3 digits coding)
#' 100 generations simulated
#'
#' @format An object of type 'genind'
#' @references \insertRef{landguth2010cdpop}{graph4lg}
#' @details The simulation was made with CDPOP during 100 generations.
#' Dispersal was possible between the 10 populations. Its probability depended
#' on the cost distance between populations, calculated on a simulated
#' resistance surface (raster). Mutations were not possible. There
#' were initially 600 alleles in total (many disappeared because of drift).
#' Population stayed constant
#' with a sex-ratio of 1. Generations did not overlap.
#' This simulation includes a part of stochasticity and these data result
#' from only 1 simulation run.
#' @examples
#' data("data_ex_genind")
#' length(unique(data_ex_genind@pop))
"data_ex_genind"


#' data_ex_gstud genetic dataset
#'
#' Genetic dataset from genetic simulation on CDPOP
#' 200 individuals, 10 populations
#' 20 microsatellite loci (3 digits coding)
#' 100 generations simulated
#'
#' @format A 'data.frame' with columns:
#' \describe{
#' \item{ID}{Individual ID}
#' \item{POP}{Population name}
#' \item{LOCI-1 to LOCI-20}{20 loci columns with microsatellite data with
#' 3 digits coding, alleles separated by ":", and blank missing data
#' (class 'locus' from \pkg{gstudio})}
#' }
#' @examples
#' data("data_ex_gstud")
#' str(data_ex_gstud)
#' length(unique(data_ex_gstud$POP))
"data_ex_gstud"


################################################################################

#' data_ex_loci genetic dataset
#'
#' Genetic dataset from genetic simulation on CDPOP
#' 200 individuals, 10 populations
#' 20 microsatellite loci (3 digits coding)
#' 100 generations simulated
#'
#' @format An object of class 'loci' and 'data.frame' with the columns :
#' \describe{
#' \item{population}{Population name}
#' \item{Other columns}{20 loci columns with microsatellite data with
#' 3 digits coding, alleles separated by "/", and missing data noted "NA/NA"}
#' }
#' Row names correspond to individuals' ID
#' @examples
#' data("data_ex_loci")
#' length(unique(data_ex_loci$population))
"data_ex_loci"



############################################################################

#' pts_pop_ex : details on simulated populations
#'
#' Simulation dataset
#' 10 populations located on a simulated landscape
#'
#' @format An object of class 'data.frame' with the following columns :
#' \describe{
#' \item{ID}{Population ID of the 10 populations}
#' \item{x}{Site longitude (RGF93)}
#' \item{y}{Site latitude (RGF93)}
#' }
#' @references \insertRef{landguth2010cdpop}{graph4lg}
#' There are as many rows as there are sampled populations.
#' @examples
#' data("pts_pop_ex")
#' str(pts_pop_ex)
"pts_pop_ex"



############################################################################

#' data_tuto : data used to generate the vignette
#'
#' Data used to generate the vignette
#'
#' @format Several outputs or inputs to show how the package works in a list
#' \describe{
#' \item{dmc}{Output of the function 'dist_max_corr'}
#' \item{graph_ci}{Genetic independence graph example}
#' \item{mat_dps}{Genetic distance matrix example}
#' \item{mat_pg}{Second genetic distance matrix example}
#' }
#' @examples
#' data("data_tuto")
#' mat_dps <- data_tuto[[1]]
#' str(mat_dps)
"data_tuto"


