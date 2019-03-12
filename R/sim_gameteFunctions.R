#' Simulate crossover positions
#'
#' Simulate crossover positions along a chromatid bundle. \strong{This function will become an internal function}.
#'
#' Simulate the locations of chiasmata along a chromatid bundle according to the model proposed by Voorrips.  Voorrips proposes to use a gamma distribution with shape 2.63 and rate parameter 2*2.63 to model the distance between successive chiasmata.  To use this model, Voorrips notes that we must incorporate a burn-in process for the first chiasmata location since "in the case of chiasmata interference [...] the position of the next chiasmata location is affected by the position of the previous chiasmata."   The burn-in is accomplished by
#' \enumerate{
#' \item Setting the burn-in distance.  From Voorrips code, and verified by my own investigations a burn-in of 1000 cM seems sufficent.
#' \item Generate location of first chiasmata (from burn-start) from an exponential distribution with mean 50 cM. (Why exponential and not gamma? I think this is to do with stationary renewal processes)
#' \item Generate successive chiasmata from a gamma(shape = 2.63, rate = 2*2.63) distribution until a chiasmata location exceeds the chromosome start position, this is used as the position of the first chiasmata.
#'}
#' After we simulate the position of the first chiasmata, successive chiasmata are drawn from a  gamma(shape = 2.63, rate = 2*2.63) distribution until a proposed chiasmata exceeds the end of the chromosome.
#'
#'
#' NOTE: To simulate chiasmata formation \strong{without chiasmata interference} (i.e. Haldane's model) we assume that the distance between successive chiasmata is distributed as an exponential distribution with mean 50 cM.  This can be accomplished by setting \code{burn_in = 0} and \code{gamma_params = c(1, 2)}.
#'
#' @references Roeland E. Voorrips, Chris A Maliepaard. (2012), \emph{The simulation of meiosis in diploid and tetraploid organisms using various genetic models}. BMC Bioinformatics, 13:248.
#'
#' @inheritParams sim_RVstudy
#' @param chrom_map Data.frame with 1 row and 2 columns. The two columns represent the start and stop positions (in cM) over which to simulate recombination.
#'
#' @importFrom stats rexp
#' @importFrom stats rgamma
#'
#' @return A list of chiasmata postions.
#' @keywords internal
sim_chiasmataPositions <- function(chrom_map,
                                   burn_in = 1000,
                                   gamma_params = c(2.63, 2.63/0.5)){

  # generate first chiasmata position using burn in process
  # suggested in Voorrips 2012, see distToFirstChiasma in Multivalent.java
  # and ranExp, ranGamma, and ranDistInterference in Tools.java
  #10 as suggested by Voorrips, and then multiply by 100 to convert to cM
  burnDist <- burn_in*0.5
  try_pos <- rexp(1, rate = 2)*100

  while(length(try_pos[which(try_pos > burnDist)]) == 0){
    try_pos <- c(try_pos, try_pos[length(try_pos)] +
                   cumsum(rgamma(5 + (burnDist*gamma_params[2])/(gamma_params[1]*100),
                                 shape = gamma_params[1],
                                 rate = gamma_params[2])*100)
                 )
  }

  chiasmata_pos <- try_pos[min(which(try_pos > burnDist))] - burnDist + chrom_map[1,1]

  #simulate chiasmata along the length of the chromosome
  while (length(chiasmata_pos[which(chiasmata_pos >= chrom_map[1, 2] )]) == 0){
    chiasmata_pos <- c(chiasmata_pos,
                       chiasmata_pos[length(chiasmata_pos)] +
                         cumsum(rgamma(max(5, 5*round((diff(as.numeric(chrom_map)))/50)),
                                       shape = gamma_params[1],
                                       rate = gamma_params[2])*100)
                       )
  }

  chiasmata_pos <- chiasmata_pos[which(chiasmata_pos < chrom_map[1, 2] )]

  return(chiasmata_pos)

}



#' Simulate recombination among chromatids.
#'
#' Simulate recombination among a bundle of four chromatids.  \strong{This function will likely become an internal function}.
#'
#' Given the possible chiasmata positions returned from \code{sim_chiasmataPostions}, we randomly select two non-sister chromatids to participate in each recombination event.  We assume no chromatid interference so that the non-sister chromatids participating in a crossover event are independent of those chosen in previous crossover events.
#'
#' After simulating recombination among the bundle, we simulate meiosis I and II, by assigning a group identifier to each hapliod: refercence (Thompson 2000)
#' \itemize{
#' \item (Meiosis I: Single Cell to Two Cells) After recombination we assign bivalents to one of the two daughter cells with equal probability.  Remember recombination has already occurred, so we identify sister chromatids by the their centromeres (location specified by user).  This process is occurs independently for different chromosomes.
#' \item (Meiosis II: Each cell from meiosis I splits into two gametes) Each pair of homologous chromosomes are separated into two gametes with equal probability and independently from the assortment of non-homologous chromosome.
#' }
#'
#' @param num_chiasmata Numeric. The number of chiasmata simulated among the chromatid bundle.
#' @param allele_IDs List of length 2. The identification numbers for the respective paternal and maternal alleles of the individual for whom we wish to simulate recombination. (Can accomodate numeric or string entries)
#'
#' @return haploid_mat. A matrix with rows representing recombined haplotypes along with an identifier that defines which group each haploid will be associated with after meiosis II.
#'
#' @references Thompson, E. (2000). \emph{Statistical Inference from Genetic Data on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and Statistics, 6, I-169. Retrieved from http://www.jstor.org.proxy.lib.sfu.ca/stable/4153187
#'
#' @seealso \code{\link{sim_chiasmataPositions}}
#'
#' @keywords internal
sim_haploidFormation <- function(num_chiasmata,
                                 allele_IDs) {

  #NOTE: At this point we have not yet simulated recombination,
  #we only know the total number of chiasmata as well as how many
  #occur before the centromere.  Their positions are not actually
  #needed for this step.  In this step we are choosing which of the
  #sister chromatids will participate in each crossover.
  #When this function is used by sim_gameteFormation we will re-append
  #the information regarding the positions of the chiasmata.

  #each column in haploid_mat represents the alleles on either side of a chiasmata
  haploid_mat <- matrix(rep(allele_IDs, each = 2*(num_chiasmata + 1)),
                        nrow = 4, byrow = T)

  #choose non-sister chromatids to participate in each chiasmata event
  choose_gam = matrix(c(sample(c(1, 2), num_chiasmata, replace = T),
                        sample(c(3, 4), num_chiasmata, replace = T)), ncol = 2)

  #swap allele sequences at sucessive chiasmata locations from left to right
  if(num_chiasmata > 0){
    for (i in 1:num_chiasmata) {
      haploid_mat[choose_gam[i, ], c(1: i)] <- haploid_mat[rev(choose_gam[i, ]), c(1: i)]
    }
  }

  #store as dataframe
  haploid_mat <- as.data.frame(haploid_mat)

  #assign gamete group (this process is equivalent to meiosis I and II
  #after we select a single gamete for inheritance, see notes.)
  haploid_mat$gamete_grp <- sample(c("A", "B", "C", "D"), size = 4, replace = FALSE)

  return(haploid_mat)

}

#' Simulate formation of gametes.
#'
#' Simulate formation of gametes. \strong{This function will likely become an internal function}.
#'
#' @inheritParams sim_haploidFormation
#' @inheritParams sim_chiasmataPositions
#' @inheritParams sim_RVstudy
#'
#' @param chrom_map Data.frame.  A data.frame consisting of three columns: column 1 contains the chromosome numbers, column 2 start postion of chromosome (in cM), column 3 end position of chromosome (in cM).
#'
#' @return  A list containing the following:
#' @return \code{chrom_haps} A list of dataframes, each dataframe is a set of four recombined haplotypes for a single chromosome (in the order specified in \code{chrom_map}), each with a gamete group identifier column.
#' @return \code{gamete_group} A list of lists, each list contains the crossover positions for a single chromosome (in the order specified in \code{chrom_map}).
#' @keywords internal
sim_gameteFormation <- function(chrom_map, allele_IDs,
                                burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)) {

  #simulate chiasmata positions
  chrom_chiasmataPos <- lapply(c(1:nrow(chrom_map)),
                               function(x){
                                 sim_chiasmataPositions(chrom_map = chrom_map[x, -1],
                                                        burn_in, gamma_params)
                                 })

  #simulate haploid and gamete formation
  chrom_haps <- lapply(c(1:nrow(chrom_map)),
                       function(x){
                         sim_haploidFormation(num_chiasmata = length(chrom_chiasmataPos[[x]]),
                                              allele_IDs)
                         })

  fun_return <- list(chrom_haps, chrom_chiasmataPos)

  return(fun_return)

}

