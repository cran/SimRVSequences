#' Constructor function for an object of class famStudy
#'
#' @param study_data The list of items returned by the \code{sim_RVstudy} function.
#'
#' @return an object of class \code{famStudy}.
#' @keywords internal
famStudy <- function(study_data) {
  class(study_data) <- c("famStudy", class(study_data))
  return(study_data)
}

#' Check to see if object is of class famStudy
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{ped}.
#'
#' @keywords internal
is.famStudy <- function(x) {
  return(inherits(x, "famStudy"))
}

#' Summary function for objects of class famStudy
#'
#' Summary function for objects of class \code{famStudy}, i.e. objects returned by the \code{\link{sim_RVstudy}} function.
#'
#' The \code{summary.famStudy} function returns a list containing two items.  The first item, \code{fam_allele} \code{_count}, is a matrix that contains counts of the SNVs shared by the disease-affected relatives in each pedigree. This matrix will contain a row of counts for each pedigree in the supplied \code{famSutdy} object.  The first column in \code{fam_allele_count} is named \code{FamID} and identifies each pedigree by their family identification number.  The remaining columns in \code{fam_allele_count} are named according to the respective marker names of the shared SNVs.
#'
#' The second item returned by \code{summary.famStudy} is a data frame named \code{pathway_count}, which catalogs the SNVs shared among disease-affected study participants.  This data frame contains the following variables:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{chrom} \tab numeric \tab chromosome identification number \cr
#' \code{position} \tab numeric \tab the position of the SNV \cr
#' \code{marker} \tab character \tab a unique character identifier for the SNV \cr
#' \code{total} \tab numeric \tab the number of SNV copies observed in disease-affected study participants. \cr
#' \code{is_crv} \tab logical \tab  identifies causal rare variants (cRVs) as \code{TRUE} \cr
#' \code{pathwaySNV} \tab logical \tab identifies SNVs located within the pathway of interest as \code{TRUE}. \cr
#' }
#'
#' Please note, the variable \code{pathwaySNV} is omitted when missing from the \code{SNV_map} data frame in the \code{famStudy} object. See \code{\link{sim_RVstudy}} for more details.
#'
#'
#' @param object An object of class \code{famStudy}, returned by the \code{sim_RVstudy} function.
#' @param ... additional arguments passed to other methods.
#' @importFrom Matrix colSums
#'
#' @return \item{\code{fam_allele_count} }{A matrix that contains counts of the SNVs shared by the disease-affected relatives in each pedigree.}
#' @return \item{\code{pathway_count} }{A data frame that catalogs the SNVs shared among disease-affected study participants. See details.}
#'
#' @seealso \code{\link{sim_RVstudy}}
#' @export
#'
#' @examples
#' library(SimRVSequences)
#'
#' # load pedigree, haplotype, and mutation data
#' data(study_peds)
#' data(EXmuts)
#' data(EXhaps)
#'
#' # create variable is_CRV in EXmuts to identify the causal
#' # rare variants from which to sample familial cRVs.
#' EXmuts$is_CRV = FALSE
#' EXmuts$is_CRV[c(26, 139, 223, 228, 472)] = TRUE
#'
#' # supply required inputs to the sim_RVstudy function
#' seqDat = sim_RVstudy(ped_files = study_peds,
#'                      SNV_data = SNVdata(Haplotypes = EXhaps,
#'                                         Mutations = EXmuts))
#'
#' # to count the number of SNVs shared by the disease-affected
#' # relatives in each pedigree, supply the output returned by
#' # sim_RVstudy to the summary function
#' summary(seqDat)
#'
#'
summary.famStudy <- function(object, ...){
  if (!is.famStudy(object)) {
    stop("\n Expecting a object of class famStudy returned by the sim_RVstudy function.")
  }

  if (all(object$ped_files$affected == FALSE)) {
    warning("There are no disease-affected study participants for whom to obtain summary data.")
  }

  Fids <- sort(unique(object$ped_files$FamID))

  #count the number of SNVs shared by the disease affected relatives by family
  aff_allele_counts <- lapply(Fids, function(x){
    affected_allele_count(ped_haps = object$ped_haplos[object$haplo_map$FamID == x, ],
                          hap_map  = object$haplo_map[object$haplo_map$FamID == x, ],
                          ped_file = object$ped_files[object$ped_files$FamID == x, ])
  })

  fam_allele_count <- do.call(rbind, aff_allele_counts)
  fam_allele_count <- cbind(Fids, fam_allele_count)
  colnames(fam_allele_count) <- c("FamID", object$SNV_map$marker)

  #create a data frame that stores sharing among study members by SNV
  pathway_count <- data.frame(chrom = object$SNV_map$chrom,
                              position = object$SNV_map$position,
                              marker = object$SNV_map$marker,
                              total = as.numeric(colSums(fam_allele_count[, -1])),
                              is_CRV = object$SNV_map$is_CRV)

  if (!is.null(object$SNV_map$pathwaySNV)) {
    pathway_count$pathwaySNV = object$SNV_map$pathwaySNV
  }

  #reduce to the SNVs carried by at least one affected study participant
  pathway_count <- pathway_count[which(pathway_count$total != 0), ]

  #remove SNVs not carried by affecteds
  fam_allele_count <- fam_allele_count[, which(colSums(fam_allele_count) != 0)]
  return(list(fam_allele_count = fam_allele_count,
              pathway_count = pathway_count))
}



#' Determine total number of alleles shared by affecteds in a family
#'
#' @param ped_haps sparse matrix.  The familial haplotype data.
#' @param hap_map data frame. Mapping data: maps individuals in \code{ped_file} to the haplotypes in \code{ped_haps}.
#' @param ped_file data frame.  Pedigree file, a for single family; i.e. not an entire study of ped files.
#'
#' @importFrom Matrix colSums
#'
#' @return A list of sharing counts, in the same order as the SNVs in \code{ped_haps}.
#'
#' @keywords internal
affected_allele_count <- function(ped_haps, hap_map, ped_file){
  #determine the locations (rows) of the affecteds in ped_haps
  aff_IDs <- ped_file$ID[ped_file$affected]
  aff_rows <- which(hap_map$ID %in% aff_IDs)

  total_count <- colSums(ped_haps[aff_rows, ])
  return(total_count)
}


