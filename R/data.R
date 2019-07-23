################################################################################
################################################################################

##### SETOPHAGA PLUMBEA ########################################################

################################################################################

#' data_pc_genind genetic dataset
#'
#' Setophaga plumbea genetic dataset
#' Data collected in French Guadeloupe by A. Khimoun and S. Garnier
#' 432 individuals, 20 populations
#' 12 microsatellite loci (3 digits coding)
#' 169 different alleles
#'
#' @format An object of type 'genind'
#' @references \insertRef{khimoun2017landscape}{graph4lg}
#' @source \url{https://datadryad.org/resource/doi:10.5061/dryad.v8324}
#' @examples
#' data("data_pc_genind")
#' length(unique(data_pc_genind@pop))
"data_pc_genind"

############################################################################


#' data_pc_gstud genetic dataset
#'
#' Setophaga plumbea genetic dataset
#' Data collected in French Guadeloupe by A. Khimoun and S. Garnier
#' 432 individuals, 20 populations
#' 12 microsatellite loci (3 digits coding)
#' 169 different alleles
#'
#' @format A 'data.frame' with columns:
#' \describe{
#' \item{Species}{Species identifier: in this case 'PC'}
#' \item{Cluster}{Population name}
#' \item{ID}{Individual ID}
#' \item{id}{Population ID used during the sampling stage}
#' \item{Latitude}{Sampling site's (population) latitude (WGS 84 - UTM 20 N)}
#' \item{Longitude}{Sampling site's (population) longitude (WGS 84 - UTM 20 N)}
#' \item{DkiD102 to DkiD12}{12 loci's columns with microsatellites' data with
#' 3 digits coding, alleles separated by ":", and blank missing data
#' (class 'locus' from \pkg{gstudio})}
#' }
#' @references \insertRef{khimoun2017landscape}{graph4lg}
#' @source \url{https://datadryad.org/resource/doi:10.5061/dryad.v8324}
#' @examples
#' data("data_pc_gstud")
#' str(data_pc_gstud)
#' length(unique(data_pc_gstud$Cluster))
"data_pc_gstud"


################################################################################

#' data_pc_gpop genetic dataset
#'
#' Setophaga plumbea genetic dataset
#' Data collected in French Guadeloupe by A. Khimoun and S. Garnier
#' 432 individuals, 20 populations
#' 12 microsatellite loci (3 digits coding)
#' 169 different alleles
#'
#' @format A 'data.frame' with the format needed to save as a GENEPOP software
#' file :
#' \describe{
#' \item{1st column}{Individual identifier beginning with population name and
#' ending with a comma}
#' \item{Other columns}{12 loci's columns with microsatellites' data with
#' 3 digits coding, alleles not separated, and missing data noted "000000"}
#' }
#' Each group of individuals from the same population are separated by a line
#' with only the character value "Pop" in firts column.
#' First lines include loci names and file name
#' @references \insertRef{khimoun2017landscape}{graph4lg}
#' @source \url{https://datadryad.org/resource/doi:10.5061/dryad.v8324}
#' @examples
#' data("data_pc_gpop")
#' str(data_pc_gpop)
"data_pc_gpop"

################################################################################

#' data_pc_loci genetic dataset
#'
#' Setophaga plumbea genetic dataset
#' Data collected in French Guadeloupe by A. Khimoun and S. Garnier
#' 432 individuals, 20 populations
#' 12 microsatellite loci (3 digits coding)
#' 169 different alleles
#'
#' @format An object of class 'loci' and 'data.frame' with the columns :
#' \describe{
#' \item{population}{Population name}
#' \item{Other columns}{12 loci's columns with microsatellites' data with
#' 3 digits coding, alleles separated by "/", and missing data noted "NA/NA"}
#' }
#' Row names correspond to individuals' ID
#' @references \insertRef{khimoun2017landscape}{graph4lg}
#' @source \url{https://datadryad.org/resource/doi:10.5061/dryad.v8324}
#' @examples
#' data("data_pc_loci")
#' length(unique(data_pc_loci$population))
"data_pc_loci"

################################################################################

#' data_pc_str genetic dataset
#'
#' Setophaga plumbea genetic dataset
#' Data collected in French Guadeloupe by A. Khimoun and S. Garnier
#' 432 individuals, 20 populations
#' 12 microsatellite loci (3 digits coding)
#' 169 different alleles
#'
#' @format An object of class 'data.frame' ready to be saved as a STRUCTRE
#' software file. The resulting text file would include the following columns :
#' \describe{
#' \item{1st column}{Individual names}
#' \item{2nd column}{Population names}
#' \item{Other columns}{12 loci's columns with microsatellites' data with
#' 3 digits coding, and missing data noted "-9"}
#' }
#' Each individual is displayed on two rows, with one allele on each row.
#' There are as many rows as twice the number of individuals.
#' @references \insertRef{khimoun2017landscape}{graph4lg}
#' @source \url{https://datadryad.org/resource/doi:10.5061/dryad.v8324}
#' @examples
#' data("data_pc_str")
"data_pc_str"


################################################################################

#' pts_pop_pc : details on Setophaga plumbea populations
#'
#' Setophaga plumbea dataset
#' Data from French Guadeloupe
#' 20 populations sampled by A. Khimoun and S. Garnier
#'
#' @format An object of class 'data.frame' with the following columns :
#' \describe{
#' \item{id}{Population ID of the 20 populations, as used during the sampling
#' stage}
#' \item{Longitude}{Sampling site's (population) longitude (WGS 84)}
#' \item{Latitude}{Sampling site's (population) latitude (WGS 84)}
#' \item{loc}{Character string indicating whether the population is located
#' on "Basse-Terre" (BT) or on "Grande-Terre" (GT) (2 main Guadeloupe's islands)}
#' \item{id_publi}{Population ID of the 20 populations, as used in the
#' paper of Khimoun et al. (2017).}
#' }
#' @references \insertRef{khimoun2017landscape}{graph4lg}
#' @source \url{https://datadryad.org/resource/doi:10.5061/dryad.v8324}
#' There are as many rows as there are sampled populations.
#' In the genetic dataset, some populations were removed, and other were
#' grouped, based on location and numbers of sampled individuals.
#' @examples
#' data("pts_pop_pc")
#' str(pts_pop_pc)
"pts_pop_pc"


################################################################################

##### SIMULATED DATA ###########################################################

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
#' on the cost distance between populations, calculated on a simulated resistance
#' surface (raster). Mutations were not possible. There were initially 600
#' alleles in total (many disappeared because of drift). Population stayed constant
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
#' \item{x}{Site's longitude (RGF93)}
#' \item{y}{Site's latitude (RGF93)}
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

