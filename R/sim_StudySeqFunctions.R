#' Dummy sample function
#'
#' @param x A vector from which to sample
#' @param ... All other arguments in sample
#'
#' @return The randomly sampled items
#' @keywords internal
resample <- function(x, ...){
  x[sample.int(length(x), ...)]
}



#' Reduce haplos to contain non-cSNV data
#'
#' @inheritParams sim_FGenos
#'
#' @return The reduced haplotype matrix.
#'
#' @keywords internal
condition_haplos_no_cSNV <- function(haplos, RV_pool_loc){
  #determine which of the haplotypes carry a causal SNV
  RV_haps <- lapply(RV_pool_loc, function(x){
    which(haplos[, x] == 1)
  })

  #RV rows
  RV_rows <- Reduce(union, RV_haps)

  return(haplos[-RV_rows, ])
}

#' Find haplotypes that do not carry any cRVs
#'
#' @inheritParams sim_FGenos
#'
#' @return The reduced haplotype matrix.
#'
#' @keywords internal
find_no_cSNV_rows <- function(haplos, RV_pool_loc){
  #determine which of the haplotypes carry a causal SNV
  RV_haps <- lapply(RV_pool_loc, function(x){
    which(haplos[, x] == 1)
  })

  #RV rows
  any_RV_rows <- Reduce(union, RV_haps)

  #no RV rows
  no_RV_rows <- setdiff(1:nrow(haplos), any_RV_rows)

  return(no_RV_rows)
}

#' Draw Founder Genotypes from Haplotype Distribution Given Familial Risk Variant
#'
#' \strong{For internal use.}
#'
#' @param founder_ids Numeric list. The ID numbers of all founders.
#' @param RV_founder Numeric. The ID number(s) of the founder who introduced a cRV.
#' @param founder_pat_allele  Numeric list. The paternally inherited allele for each founder.
#' @param founder_mat_allele  Numeric list. The maternally inherited allele for each founder.
#' @param haplos sparseMatrix.  The sparseMatrix of genomes returned by \code{read_slimOut}.
#' @param RV_col_loc Numeric. The column location of the familial RV in haplos.
#' @param RV_pool_loc The column locations of each SNV (in haplos)  for SNVs in the pool of cRVs.
#'
#' @return list of familial founder genotypes
#'
#' @keywords internal
sim_FGenos <- function(founder_ids, RV_founder,
                       founder_pat_allele, founder_mat_allele,
                       haplos, RV_col_loc, RV_pool_loc) {

  #Determine which haplotypes carry the familial rare variant and which do not
  #Determine which haplotypes carry the familial cRV
  RV_hap_loc <- which(haplos[, RV_col_loc] == 1)

  #Determine which haplotypes do not carry ANY crv in the pool
  no_CRVrows <- find_no_cSNV_rows(haplos, RV_pool_loc)



  #here we handle the fully sporatic families
  #i.e. families that do not segregate any cSNVs
  #In this case, the haplotypes for ALL founders
  #is sampled from no_CRVhaps
  if(length(RV_founder) == 0){
    #sample all founder data from this pool
    founder_genos <- haplos[sample(x = no_CRVrows,
                                   size = 2*length(founder_ids),
                                   replace = TRUE), ]
  } else {

    #sample the paternally inherited founder haplotypes
    pat_inherited_haps <- sapply(founder_pat_allele, function(x){
      if(x == 0){
        resample(x = no_CRVrows, size = 1)
      } else {
        resample(x = RV_hap_loc, size = 1)
      }})

    #sample the maternally inherited founder haplotypes
    mat_inherited_haps <- sapply(founder_mat_allele, function(x){
      if(x == 0){
        resample(x = no_CRVrows, size = 1)
      } else {
        resample(x = RV_hap_loc, size = 1)
      }})

    #pull the sampled haplotypes from the haplos matrix
    founder_genos <- haplos[c(pat_inherited_haps, mat_inherited_haps), ]
  }

  #create IDs to associate founders to rows in founder_genos
  founder_genos_ID <- rep(founder_ids, 2)

  #re-order so that founder haplotypes appear in order
  founder_genos <- founder_genos[order(founder_genos_ID), ]
  founder_genos_ID <- founder_genos_ID[order(founder_genos_ID)]


  return(list(founder_genos, founder_genos_ID))
}

#' Remove unmutated markers from data
#'
#' Remove any markers for which all founders, in the study, are homozygous for the wild-type allele.  Since we do not model de novo mutations, it is not possible for non-founders to develop mutations at these loci.
#'
#' @param f_haps The founder haplotypes data. This is a list of family lists. By family, this contains the haplotypes for each founder (first item), and a list of ID numbers (second item) which is used to map the haplotype to the person to whom it belongs.
#' @param SNV_map data.frame. Catalogs the SNV data contained in the familial haplotypes.
#'
#' @return A list (by family) of haplotype matrices and ID vectors and the reduce marker data set.
#' @importFrom Matrix colSums
#' @keywords internal
remove_allWild <- function(f_haps, SNV_map){

  #determine which columns are all zero in founder haplotype data.
  #These are markers at which no one in ped_files will carry a SNV
  #(i.e. all family members will carry the wild type allele) since no
  #founders have mutated alleles to pass on.
  #We will reduce the size of the data by removing these superfluous

  #Determine all-wild columns by family
  wild_col <- lapply(f_haps, function(x){
    which(colSums(x[[1]]) == 0)
  })

  #determine all wild columns overall, i.e. for the study
  remove_cols <- Reduce(intersect, wild_col)

  #remove all wild columns from founder haplotypes
  red_haps <- lapply(f_haps, function(x){
    list(x[[1]][, -remove_cols],
         x[[2]])
  })

  #remove all wild columns from SNV_map
  #RECALL: rows of SNV_map are columns in haplotype data
  SNV_map <- SNV_map[-remove_cols, ]
  SNV_map$colID <- seq(1:nrow(SNV_map))
  row.names(SNV_map) = NULL

  return(list(red_haps, SNV_map))
}

#' Simulate sequence data for a sample of pedigrees
#'
#' Simulate single-nucleotide variant (SNV) data for a sample of pedigrees.
#'
#' The \code{sim_RVstudy} function is used to simulate single-nucleotide variant (SNV) data for a sample of pedigrees.  Please note: this function is NOT appropriate for users who wish to simulate genotype conditional on phenotype.  Instead, \code{sim_RVstudy} employs the following algorithm.
#'
#' \enumerate{
#' \item For each pedigree, we sample a single \strong{causal rare variant (cRV)} from a pool of SNVs specified by the user.
#' \item Upon identifying the familial cRV we sample founder haplotypes from haplotype data conditional on the founder's cRV status at the familial cRV locus.
#' \item Proceeding forward in time, from founders to more recent generations, for each parent/offspring pair we:
#' \enumerate{
#' \item simulate recombination and formation of gametes, according to the model proposed by Voorrips and Maliepaard (2012), and then
#' \item perform a conditional gene drop to model inheritance of the cRV.
#' }}
#'
#' It is important to note that due to the forwards-in-time algorithm used by \code{sim_RVstudy}, \strong{certain types of inbreeding and/or loops cannot be accommodated}. Please see examples.
#'
#' For a detailed description of the model employed by \code{sim_RVstudy}, please refer to section 6 of the vignette.
#'
#' The data frame of pedigrees, \code{ped_files}, supplied to \code{sim_RVstudy} must contain the variables:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{FamID} \tab numeric \tab family identification number\cr
#' \code{ID} \tab numeric \tab individual identification number\cr
#' \code{sex} \tab numeric \tab sex identification variable: \code{sex = 0} for males, and \code{sex = 1} females. \cr
#' \code{dadID} \tab numeric \tab identification number of father \cr
#' \code{momID} \tab numeric \tab identification number of mother \cr
#' \code{affected} \tab logical \tab disease status indicator: set \code{affected = TRUE} if individual has disease.\cr
#' \code{DA1} \tab numeric \tab paternally inherited allele at the cRV locus: \cr
#'  \tab \tab  \code{DA1 = 1} if the cRV is inherited, and \code{0} otherwise. \cr
#'  \code{DA2} \tab numeric \tab maternally inherited allele at the cRV locus: \cr
#'  \tab \tab \code{DA2 = 1} if the cRV is inherited, and \code{0} otherwise.\cr
#' }
#'
#' If \code{ped_files} does not contain the variables \code{DA1} and \code{DA2} the pedigrees are assumed to be fully sporadic.  Hence, the supplied pedigrees will not segregate any of the SNVs in the user-specified pool of cRVs.
#'
#' Pedigrees simulated by the \code{\link{sim_RVped}} and \code{\link{sim_ped}} functions of the \code{SimRVPedigree} package are properly formatted for the \code{sim_RVstudy} function.  That is, the pedigrees generated by these functions contain all of the variables required for \code{ped_files} (including \code{DA1} and \code{DA2}).
#'
#' The data frame \code{SNV_map} catalogs the SNVs in \code{haplos}. The variables in \code{SNV_map} must be formatted as follows:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{colID} \tab numeric \tab associates the rows in \code{SNV_map} to the columns of \code{haplos}\cr
#' \code{chrom} \tab numeric \tab the chromosome that the SNV resides on\cr
#' \code{position} \tab numeric \tab is the position of the SNV in base pairs when argument \cr
#' \tab \tab \code{pos_in_bp = TRUE} or centiMorgan when \code{pos_in_bp = FALSE}\cr
#' \code{marker} \tab character \tab (Optional) a unique character identifier for the SNV. \cr
#' \tab \tab If missing this variable will be created from \code{chrom} and \code{position}. \cr
#' \code{pathwaySNV} \tab logical \tab (Optional) identifies SNVs located within the pathway of interest as \code{TRUE} \cr
#' \code{is_CRV} \tab logical \tab  identifies causal rare variants (cRVs) as \code{TRUE}.\cr
#' }
#'
#' Please note that when the variable \code{is_CRV} is missing from \code{SNV_map}, we sample a single SNV to be the causal rare variant for all pedigrees in the study, which is identified in the returned \code{famStudy} object.
#'
#' @param ped_files Data frame. A data frame of pedigrees for which to simulate sequence data, see details.
#' @param SNV_data SNVdata. An object of class \code{SNVdata} created by \code{\link{SNVdata}}.
#' @param affected_only Logical. When \code{affected_only = TRUE}, we only simulate SNV data for the disease-affected individuals and the family members that connect them along a line of descent.  When \code{affected_only = FALSE}, SNV data is simulated for the entire study. By default, \code{affected_only = TRUE}.
#' @param pos_in_bp Logical. This argument indicates if the positions in \code{SNV_map} are listed in base pairs.  By default, \code{pos_in_bp = TRUE}. If the positions in \code{SNV_map} are listed in centiMorgan please set \code{pos_in_bp = FALSE} instead.
#' @param remove_wild Logical.  When \code{remove_wild = TRUE} the data is reduced by removing SNVs which are not observed in any of the study participants; otherwise if \code{remove_wild} \code{= FALSE} no data reduction occurs.  By default, \code{remove_wild = TRUE}.
#' @param gamma_params Numeric list of length 2. The respective shape and rate parameters of the gamma distribution used to simulate distance between chiasmata.  By default, \code{gamma} \code{_params} \code{= c(2.63, 2*2.63)}, as discussed in Voorrips and Maliepaard (2012).
#' @param burn_in Numeric. The "burn-in" distance in centiMorgan, as defined by Voorrips and Maliepaard (2012), which is required before simulating the location of the first chiasmata with interference. By default, \code{burn_in = 1000}.
#' The burn in distance in cM. By default, \code{burn_in = 1000}.
#' @param SNV_map This argument has been deprecated. Users now supply objects of class \code{SNVdata} to argument \code{SNV_data}.
#' @param haplos This argument has been deprecated. Users now supply objects of class \code{SNVdata} to argument \code{SNV_data}.
#'
#' @return  A object of class \code{famStudy}.  Objects of class \code{famStudy} are lists that include the following named items:
#' @return \item{\code{ped_files}}{A data frame containing the sample of pedigrees for which sequence data was simulated.}
#' @return \item{\code{ped_haplos}}{A sparse matrix that contains the simulated haplotypes for each pedigree member in \code{ped_files}.}
#' @return \item{\code{haplo_map}}{A data frame that maps the haplotypes (i.e. rows) in \code{ped_haplos} to the individuals in \code{ped_files}.}
#' @return \item{\code{SNV_map}}{A data frame cataloging the SNVs in \code{ped_haplos}.}
#' @return Objects of class \code{famStudy} are discussed in detail in section 5.2 of the vignette.
#'
#' @references Roeland E. Voorrips and Chris A Maliepaard. (2012). \emph{The simulation of meiosis in diploid and tetraploid organisms using various genetic models}. BMC Bioinformatics, 13:248.
#'
#' @references Christina Nieuwoudt, Angela Brooks-Wilson, and Jinko Graham. (2019). \emph{SimRVSequences: an R package to simulate genetic sequence data for pedigrees.} <doi:10.1101/534552>.
#'
#' @export
#'
#' @seealso \code{\link{sim_RVped}}, \code{\link{read_slim}}, \code{\link{summary.famStudy}}
#'
#' @examples
#' library(SimRVSequences)
#'
#' #load pedigree, haplotype, and mutation data
#' data(study_peds)
#' data(EXmuts)
#' data(EXhaps)
#'
#' # create variable 'is_CRV' in EXmuts.  This variable identifies the pool of
#' # causal rare variants  from which to sample familial cRVs.
#' EXmuts$is_CRV = FALSE
#' EXmuts$is_CRV[c(26, 139, 223, 228, 472)] = TRUE
#'
#' # create object of class SNVdata
#' my_SNVdata <- SNVdata(Haplotypes = EXhaps,
#'                       Mutations = EXmuts)
#'
#' #supply required inputs to the sim_RVstudy function
#' seqDat = sim_RVstudy(ped_files = study_peds,
#'                      SNV_data = my_SNVdata)
#'
#'
#' # Inbreeding examples
#' # Due to the forward-in-time model used by sim_RVstudy certain types of
#' # inbreeding and/or loops *may* cause fatal errors when using sim_RVstudy.
#' # The following examples demonstrate: (1) imbreeding that can be accommodated
#' # under this model, and (2) when this limitation is problematic.
#'
#' # Create inbreeding in family 1 of study_peds
#' imb_ped1 <- study_peds[study_peds$FamID == 3, ]
#' imb_ped1[imb_ped1$ID == 18, c("momID")] = 7
#' plot(imb_ped1)
#'
#' # Notice that this instance of inbreeding can be accommodated by our model.
#' seqDat = sim_RVstudy(ped_files = imb_ped1,
#'                      SNV_data = my_SNVdata)
#'
#' # Create different type of inbreeding in family 1 of study_peds
#' imb_ped2 <- study_peds[study_peds$FamID == 3, ]
#' imb_ped2[imb_ped1$ID == 8, c("momID")] = 18
#' plot(imb_ped2)
#'
#' # Notice that inbreeding in imb_ped2 will cause a fatal
#' # error when the sim_RVstudy function is executed
#' \dontrun{
#' seqDat = sim_RVstudy(ped_files = imb_ped2,
#'                      SNV_data = my_SNVdata)
#' }
#'
sim_RVstudy <- function(ped_files, SNV_data,
                        affected_only = TRUE,
                        remove_wild = TRUE,
                        pos_in_bp = TRUE,
                        gamma_params = c(2.63, 2.63/0.5),
                        burn_in = 1000,
                        SNV_map = NULL, haplos = NULL){

  if (!is.null(SNV_map) | !is.null(haplos)) {
    stop("Arguments 'SNV_map' and 'haplos' have been deprecated.\n Instead, please supply to argument 'SNV_data' an object of class SNVdata.  Execute help(SNVdata) for more information." )
  }

  if (!is.SNVdata(SNV_data)) {
    stop("Expecting SNV_data to be an object of class SNVdata")
  }

  #check to see if DA1 and DA2 are both missing, if so
  #assume fully sporadic and issue warning
  if (is.null(ped_files$DA1) & is.null(ped_files$DA2)) {
    ped_files$DA1 <- 0
    ped_files$DA2 <- 0
    warning("\n The variables DA1 and DA2 are missing from ped_files. \n Assuming fully sporadic ... \n...setting DA1 = DA2 = 0 for all pedigrees.")
  }


  #check ped_files for possible issues
  check_peds(ped_files)

  #assign generation number if not included in ped_file
  if(!"Gen" %in% colnames(ped_files)){
    ped_files$Gen <- unlist(lapply(unique(ped_files$FamID), function(x){assign_gen(ped_files[ped_files$FamID == x, ])}))
  }

  #quick and dirty assignment to see if new format creates problems down to line...
  #TODO: change all function definitions to accept new SNVdata object
  #NOTE: All of the same checks were performed on input data, the new SNVdata object
  SNV_map = SNV_data$Mutations
  haplos = SNV_data$Haplotypes

  #check to see that the sample contains affected relatives when the
  #affected_only setting is used
  if (affected_only & all(ped_files$affected  == FALSE)) {
    stop("\n There are no disease-affected relatives in this sample of pedigrees. \n To simulate data for pedigrees without disease-affected relatives use affected_only = FALSE.")
  }

  #collect list of FamIDs
  FamIDs <- unique(ped_files$FamID)

  #check for pedigree formatting issues
  for (i in FamIDs){
    check_ped(ped_files[ped_files$FamID == i, ])
  }

  #Reduce to affected-only pedigrees
  if (affected_only) {
    #reduce pedigrees to contain only disease-affected relative and
    #the individuals who connect them along a line of descent.
    Afams <- lapply(FamIDs, function(x){
      affected_onlyPed(ped_file = ped_files[which(ped_files$FamID == x),])
    })

    #combine the reduced pedigrees
    ped_files <- do.call("rbind", Afams)

    #check to see if any pedigrees were removed due to lack of
    #disease affected relatives and issue warning for removed pedigrees
    removed_peds <- setdiff(FamIDs, unique(ped_files$FamID))

    if (length(removed_peds) > 0){
      FamIDs <- unique(ped_files$FamID)
      warning("\n There are no disease-affected relatives in the pedigrees with FamID: ",
              paste0(removed_peds, collapse = ", "),
              "\n These pedigrees have been removed from ped_files.")
    }

  }

  #sampling from RV markers
  #to determine familial RV locus
  #NOTE: IF POSSIBLE SNV NOT SPECIFIED A SINGLE SNV IS CHOSEN AS
  #      CAUSAL FOR ALL FAMILIES IN STUDY.  DOES NOT CONSIDER ALLELE FREQUENCY.
  if (is.null(SNV_map$is_CRV)) {
    SNV_map$is_CRV = FALSE
    SNV_map$is_CRV[sample(1:nrow(SNV_map), size = 1)] = TRUE
    warning("The variable is_CRV is missing from SNV_map.",
            "\n ... randomly sampling one SNV to be the cRV for all pedigrees.")
  }

  #set the sampling probabilities for causal rare variants
  #When the derived allele frequencies are provided, we sample causal rare
  #variants according to their allele frequency.  When missing, we sample
  #cRVs with equal probability.
  if ("afreq" %in% colnames(SNV_map)) {
    sample_prob <- SNV_map$afreq[SNV_map$is_CRV]/sum(SNV_map$afreq[SNV_map$is_CRV])
  } else {
    sample_prob <- rep(1/sum(SNV_map$is_CRV), sum(SNV_map$is_CRV))
    warning("The variable afreq is missing from SNV_map. \n ...sampling cRVs with equal probability.")

  }

  #sample the familial cRV from the pool of potential cRVs with replacement.
  Fam_RVs <- sample(x = SNV_map$marker[SNV_map$is_CRV],
                    size = length(FamIDs),
                    prob = sample_prob,
                    replace = TRUE)

  #Given the location of familial risk variants, sample familial founder
  #haplotypes from conditional haplotype distribution
  f_genos <- lapply(c(1:length(FamIDs)), function(x){
    sim_FGenos(founder_ids = ped_files$ID[which(ped_files$FamID == FamIDs[x]
                                                & is.na(ped_files$dadID))],
               RV_founder = ped_files$ID[which(ped_files$FamID == FamIDs[x]
                                               & is.na(ped_files$dadID)
                                               & (ped_files$DA1 + ped_files$DA2) != 0)],
               founder_pat_allele = ped_files$DA1[which(ped_files$FamID == FamIDs[x]
                                                        & is.na(ped_files$dadID))],
               founder_mat_allele = ped_files$DA2[which(ped_files$FamID == FamIDs[x]
                                                        & is.na(ped_files$dadID))],
               haplos, RV_col_loc = which(SNV_map$marker == Fam_RVs[x]),
               RV_pool_loc = SNV_map$colID[SNV_map$is_CRV])
  })

  #If desired by user, we now reduce the size of the data by removing
  #markers not carried by any member of our study.
  if (remove_wild) {
    reduced_dat <- remove_allWild(f_haps = f_genos, SNV_map)
    f_genos <- reduced_dat[[1]]
    SNV_map <- reduced_dat[[2]]
  }

  #create chrom_map, this is used to determine the segments over
  #which we will simulate genetic recombination
  chrom_map <- create_chrom_map(SNV_map)

  #convert from base pairs to centiMorgan
  if (pos_in_bp) {
    options(digits = 9)
    chrom_map$start_pos <- convert_BP_to_cM(chrom_map$start_pos)
    chrom_map$end_pos <- convert_BP_to_cM(chrom_map$end_pos)

    SNV_map$position <- convert_BP_to_cM(SNV_map$position)
  }

  #simulate non-founder haploypes via conditional gene drop
  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_seq(ped_file = ped_files[ped_files$FamID == FamIDs[x], ],
            founder_genos = f_genos[[x]],
            SNV_map, chrom_map,
            RV_marker = Fam_RVs[x],
            burn_in, gamma_params)
  })

  ped_haplos <- do.call("rbind", lapply(ped_seqs, function(x){x$ped_genos}))
  haplo_map <- do.call("rbind", lapply(ped_seqs, function(x){x$geno_map}))

  #convert back to base pairs if we converted to CM
  if (pos_in_bp) {
    options(digits = 9)
    SNV_map$position <- convert_CM_to_BP(SNV_map$position)
  }

  return(famStudy(list(ped_files = ped_files, ped_haplos = ped_haplos, haplo_map = haplo_map, SNV_map = SNV_map)))
}
