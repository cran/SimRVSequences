#' Create recombination map
#'
#' Create a recombination map that can be used with SLiM (Haller and Messer 2017).
#'
#' The Eidos program SLiM (Haller and Messer 2017) is a versatile forwards-in-time evolutionary simulator.  SLiM simulates recombination hotspots by way of a user-specified recombination map.  This recombination map may be utilized to simulate mutations over unlinked regions (i.e. in different chromosomes) or in linked but non-contiguous regions (i.e in exon-only data).  The \code{create_slimMap} function may be used to generate the recombination map required by SLiM to simulate exon-only SNV data.
#'
#' We expect that \code{exon_df} does not contain any overlapping segments.  Prior to supplying the exon data to \code{create_slimMap} users must combine overlapping exons into a single observation.  The \code{\link{combine_exons}} function may be used to accomplish this task.
#'
#' The argument \code{exon_df} must contain the following variables:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{chrom} \tab numeric \tab chromosome identification number\cr
#' \code{exonStart} \tab numeric \tab the position of the first base pair in the exon\cr
#' \code{exonStop} \tab numeric \tab the position of the last base pair in the exon\cr
#' }
#'
#' The data frame returned by \code{create_slimMap} contains variables required by SLiM to simulate exon-only data.  Additionally, the returned data frame also includes variables that are required to re-map mutations to their correct positions when importing SLiM data to \code{R}.  The variables contained in the returned data frame are described as follows.
#' \describe{
#' \item{\code{chrom}}{The chromosome number.}
#' \item{\code{segLength}}{The length of the segment in base pairs.  We assume that segments contain the positions listed in \code{exonStart} and \code{exonEnd}.  Therefore, for a combined exon segment, \code{segLength} is calculated as \code{exonEnd - exonStart + 1}.}
#' \item{\code{recRate}}{The per-site per-generation recombination rate.  Following Harris and Nielson (2016), segments between exons on the same chromosome are simulated as a single base pair with \code{rec_rate} equal to recombination rate multiplied by the number of base pairs in the segment.  For each chromosome, a single site is created between the last exon on the previous chromosome and the first exon of the current chromosome.  This site will have recombination rate 0.5 to accommodate unlinked chromosomes.}
#' \item{\code{mutRate}}{The per-site per-generation mutation rate.  Since we are interested in exon-only data, the mutation rate outside exons is set to zero.}
#' \item{\code{exon}}{A logical variable that is \code{TRUE} if the segment is an exon and \code{FALSE} otherwise.}
#' \item{\code{simDist}}{The simulated exon length, in base pairs. When \code{exon = TRUE}, \code{simDist = segLength}; however, when \code{exon = FALSE}, \code{simDist = 1} since segments between exons on the same chromosome are simulated as a single base pair.}
#' \item{\code{endPos}}{The simulated end position, in base pairs, of the segment.}
#' }
#'
#' Only three of the variables returned by \code{create_slimMap} are required by SLiM to simulate exon-only data: \code{recRate}, \code{mutRate}, and \code{endPos}.  The other variables seen in the output above are used by the \code{\link{read_slim}} function to re-map mutations to their correct positions when importing SLiM data to \code{R}.
#'
#' Please note: SLiM is written in a scripting language called Eidos. Unlike an \code{R} array, the first position in an Eidos array is 0.  Therefore, users must shift the variable \code{endPos} forward 1 unit before supplying this variable to SLiM. See example.
#'
#' @param exon_df Data frame. A data frame that contains the positions of each exon to simulate.  This data frame must contain the variables \code{chrom}, \code{exonStart}, and \code{exonEnd}.  See details.
#' @param mutation_rate Numeric.  The per-site per-generation mutation rate, assumed to be constant across the genome. By default, \code{mutation_rate= 1E-8}, as in Harris and Nielson (2016).
#' @param recomb_rate Numeric.  The per-site per-generation mutation rate, assumed to be constant across the genome. By default, \code{mutation_rate= 1E-8}, as in Harris and Nielson (2016)
#'
#' @return A recombination map that may be used in conjunction with SLiM (Haller and Messer 2017).  See details and example.
#'
#' @references Benjamin Haller and Phillip W. Messer (2017). \emph{Slim 2: Flexible, interactive forward genetic simulations}. Molecular Biology and Evolution; 34(1), pp. 230-240.
#' @references Kelly Harris and Rasmus Nielsen (2016). \emph{The genetic cost of neanderthal introgression}. Genetics, 203(2): pp. 881-891.
#'
#' @export
#'
#' @seealso \code{\link{combine_exons}}
#'
#' @examples
#' #load hg_exons data
#' data(hg_exons)
#'
#' #since the exons in hg_exons have already been combined into
#' #overlapping exons, we supply hg_exons to create_slimMap
#' slimMap <- create_slimMap(hg_exons)
#' head(slimMap)
#'
#' # restrict output to the variables required by SLiM
#' slimMap <- slimMap[, c("recRate", "mutRate", "endPos")]
#'
#' # shift endPos up by one unit
#' slimMap$endPos <- slimMap$endPos - 1
#'
#' # print first four rows of slimMap
#' head(slimMap, n = 4)
#'
create_slimMap <- function(exon_df, mutation_rate = 1E-8, recomb_rate = 1E-8){
  #split into dataset for each chromosome
  bychr <- lapply(sort(unique(exon_df$chrom)), function(x){
    exon_df[exon_df$chrom == x, ]
  })

  #compile data on introns
  int_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               segLength = c((x$exonStart[1] - 1), (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               no = c(1:nrow(x)),
               recRate = recomb_rate*c(0.5/recomb_rate, (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               mutRate = rep(0, nrow(x)),
               type = rep("intron", nrow(x)),
               exon = rep(FALSE, nrow(x)),
               simDist = rep(1, nrow(x)),
               stringsAsFactors = TRUE)
  })
  #set the recombination rate to zero for the first chromosome
  #Note: the first recombination rate is 0.5 for all other
  #chromosomes so that successive chromosomes are unlinked
  int_dist[[1]]$recRate[1] = 0

  #compile data on exons
  ex_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               segLength = c(x$exonEnd - x$exonStart + 1),
               no = c(1:nrow(x)),
               recRate = rep(recomb_rate, nrow(x)),
               mutRate = rep(mutation_rate, nrow(x)),
               type = rep("exon", nrow(x)),
               exon = rep(TRUE, nrow(x)),
               simDist = c(x$exonEnd - x$exonStart + 1),
               stringsAsFactors = TRUE)
  })

  #combine intron and exon data
  ie_dist <- lapply(1:length(int_dist), function(x){
    rbind(int_dist[[x]], ex_dist[[x]])
  })

  #order by type and number so that they appear in proper order
  ie_dist <- lapply(ie_dist, function(x){
    x[order(x$no, x$type), ]
  })

  #recombine datasets and return
  rc_map <- do.call(rbind, ie_dist)
  rc_map$endPos <- cumsum(rc_map$simDist)
  row.names(rc_map) = NULL
  return(rc_map[, -c(3, 6)])
}

#' Re-map slim mutations
#'
#' Intended for internal use
#'
#' @param mutationDF The Mutation data frame returned by reMap_mutations
#' @param recomb_map The recombination map provided to slim
#'
#' @return A re-mapped mutation data frame
#' @keywords internal
reMap_mutations <- function(mutationDF, recomb_map){
  #split into data for different chromosomes, because it
  #makes myhead hurt to think about this as 1 chromosome
  bychr <- lapply(sort(unique(recomb_map$chrom)), function(x){
    recomb_map[recomb_map$chrom == x, ]
  })

  #subset by introns for each chromosome
  bychr_int <- lapply(bychr, function(x){
    x[!x$exon, ]
  })

  #get chrom starts and stops
  chr_start <- unlist(lapply(bychr, function(x){ x$endPos[1] }))
  chr_end <- unlist(lapply(bychr, function(x){ x$endPos[nrow(x)] }))

  #determine which chromosome each mutation falls on
  mutationDF$chrom <- cut(mutationDF$position,
                          breaks = c(1, chr_end),
                          labels = FALSE)

  #renumber for chromosomes included in recomb_map
  mutationDF$chrom <- sort(unique(recomb_map$chrom))[mutationDF$chrom]


  #create separate data frames for each chromosome
  mut_by_chrom <- lapply(sort(unique(mutationDF$chrom)), function(x){
    mutationDF[mutationDF$chrom == x, ]
  })


  for(i in 1:length(mut_by_chrom)){
    #shift mutations and map so that the first intron starts at
    #position 1 in each chromosome
    mut_by_chrom[[i]]$position <- mut_by_chrom[[i]]$position - (chr_start - 1)[i]
    bychr[[i]]$newEndPos <- bychr[[i]]$endPos - (chr_start - 1)[i]

    #determine which exon each mutation falls in
    mut_by_chrom[[i]]$ex_num <- (cut(mut_by_chrom[[i]]$position,
                                     breaks = bychr[[i]]$newEndPos,
                                     labels = FALSE) + 1)/2

    #get cumulative intron distance
    bychr_int[[i]]$cumDist <- cumsum(bychr_int[[i]]$segLength)

    mut_by_chrom[[i]]$position <- mut_by_chrom[[i]]$position +
      bychr_int[[i]]$cumDist[mut_by_chrom[[i]]$ex_num] -
      mut_by_chrom[[i]]$ex_num

  }

  mut_dat <- do.call(rbind, mut_by_chrom)
  return(mut_dat)
}

#' Import SLiM data to R
#'
#'To import SLiM data into \code{R}, we provide the \code{read_slim} function, which has been tested for SLiM versions 2.0-3.1. \strong{The \code{read_slim} function is only appropriate for single-nucleotide variant (SNV) data produced by SLiM's outputFull() method.}  We do not support output in MS or VCF data format, i.e. produced by outputVCFsample() or outputMSSample() in SLiM.
#'
#' In addition to reducing the size of the data, the argument \code{keep_maf} has practicable applicability.  In family-based studies, common SNVs are generally filtered out prior to analysis.  Users who intend to study common variants in addition to rare variants may need to run chromosome specific analyses to allow for allocation of large data sets in \code{R}.
#'
#' The argument \code{recomb_map} is used to remap mutations to their actual locations and chromosomes.  This is necessary when data has been simulated over non-contiguous regions such as exon-only data.  If \code{\link{create_slimMap}} was used to create the recombination map for SLiM, simply supply the output of \code{create_slimMap} to \code{recomb_map}.  If \code{recomb_map} is not provided we assume that the SNV data has been simulated over a contiguous segment starting with the first base pair on chromosome 1.
#'
#' The data frame \code{pathway_df} allows users to identify SNVs located within a pathway of interest.  When supplied, we expect that \code{pathwayDF} does not contain any overlapping segments.  \emph{All overlapping exons in \code{pathway_df} MUST be combined into a single observation.  Users may combine overlapping exons with the \code{\link{combine_exons}} function.}
#'
#' When \code{TRUE}, the logical argument \code{recode_recurrent} indicates that recurrent SNVs should be recorded as a single observation.  SLiM can model many types of mutations; e.g. neutral, beneficial, and deleterious mutations.  When different types of mutations occur at the same position carriers will experience different fitness effects depending on the carried mutation.  However, when mutations at the same location have the same fitness effects, they represent a recurrent mutation.  Even so, SLiM stores recurrent mutations separately and calculates their prevalence independently.  When the argument \code{recode_recurrent = TRUE} we store recurrent mutations as a single observation and calculate the derived allele frequency based on their combined prevalence.  This convention allows for both reduction in storage and correct estimation of the derived allele frequency of the mutation.  Users who prefer to store recurrent mutations from independent lineages as unique entries should set \code{recode_recurrent = FALSE}.
#'
#'An object of class \code{\link{SNVdata}}, which inherits from a \code{list} and contains:
#' The \code{read_slim} function returns an object of class \code{\link{SNVdata}}, which inherits from a \code{list} and contains the following two items:
#' \enumerate{
#' \item \code{Haplotypes} A sparse matrix of class dgCMatrix (see \code{\link{dgCMatrix-class}}). The columns in {Haplotypes} represent distinct SNVs, while the rows represent individual haplotypes. We note that this matrix contains two rows of data for each diploid individual in the population: one row for the maternally ihnherited haplotype and the other for the paternally inherited haplotype.
#' \item \code{Mutations} A data frame cataloging SNVs in \code{Haplotypes}. The variables in the \code{Mutations} data set are described as follows:
#' \describe{
#' \item{\code{colID}}{Associates the rows, i.e. SNVs, in \code{Mutations} to the columns of \code{Haplotypes}.}
#' \item{\code{chrom}}{The chromosome that the SNV resides on.}
#' \item{\code{position}}{The position of the SNV in base pairs.}
#' \item{\code{afreq}}{The derived allele frequency of the SNV.}
#' \item{\code{marker}}{A unique character identifier for the SNV.}
#' \item{\code{type}}{The mutation type, as specified in the user's slim simulation.}
#' \item{\code{pathwaySNV}}{Identifies SNVs located within the pathway of interest as \code{TRUE}.}
#' }}
#'
#' Please note: the variable \code{pathwaySNV} will be omitted when \code{pathway_df} is not supplied to \code{read_slim}.
#'
#' @param file_path character.  The file path or URL of the .txt output file created by the outputFull() method in SLiM.
#' @param keep_maf numeric. The largest allele frequency for retained SNVs, by default \code{keep_maf} \code{= 0.01}.  All variants with allele frequency greater than \code{keep_maf} will be removed. Please note, removing common variants is recommended for large data sets due to the limitations of data allocation in R. See details.
#' @param recomb_map data frame. (Optional) A recombination map of the same format as the data frame returned by \code{\link{create_slimMap}}. See details.
#' @param pathway_df data frame. (Optional) A data frame that contains the positions for each exon in a pathway of interest.  See details.
#' @param recode_recurrent logical. When \code{TRUE} recurrent SNVs are cataloged a single observation;  by default, \code{recode_recurrent = TRUE}. See details.
#'
#' @return  An object of class \code{\link{SNVdata}}, which inherits from a \code{list} and contains:
#' @return \item{\code{Haplotypes} }{A sparse matrix of haplotypes. See details.}
#' @return \item{\code{Mutations}}{A data frame cataloging SNVs in \code{Haplotypes}. See details.}
#' @importFrom Matrix sparseMatrix
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#'
#' @export
#'
#' @references Haller, B., Messer, P. W. (2017). \emph{Slim 2: Flexible, interactive forward genetic simulations}. Molecular Biology and Evolution; 34(1), pp. 230-240.
#' @references Douglas Bates and Martin Maechler (2018). \strong{Matrix: Sparse and Dense Matrix Classes and Methods}.
#' \emph{R package version 1.2-14}. https://CRAN.R-project.org/package=Matrix
#'
#' @seealso \code{\link{create_slimMap}}, \code{\link{combine_exons}}, \code{\link{dgCMatrix-class}}
#'
#' @examples
#' # Specify the URL of the example output data simulated by SLiM.
#' file_url <-
#' 'https://raw.githubusercontent.com/cnieuwoudt/Example--SLiMSim/master/example_SLIMout.txt'
#' s_out <- read_slim(file_url)
#'
#' class(s_out)
#' str(s_out)
#'
#'
#' # As seen above, read_slim returns an object of class SNVdata,
#' # which  contians two items.  The first is a sparse matrix
#' # named Haplotypes, which contains the haplotypes for each indiviual in the
#' # simulation.  The second item is a data set named Mutations, which catalogs
#' # the mutations in the Haplotypes matrix.
#'
#' # View the first 5 lines of the mutation data
#' head(s_out$Mutations, n = 5)
#'
#' # view the first 20 mutations on the first 10 haplotypes
#' s_out$Haplotypes[1:10, 1:20]
#'
#'
read_slim <- function(file_path,
                      keep_maf = 0.01,
                      recomb_map = NULL,
                      pathway_df = NULL,
                      recode_recurrent = TRUE){
  #NOTE: Time to read file ~19 secs
  message("Reading Slim File")
  exDat = readLines(file_path)

  #The default output for Slim  is a .txt file with the following headings:
  # - OUT: Contains number of generations, type (i.e. "A" for autosome), and
  #   the file name
  #
  # - Populations:
  #   next line is pop description, "p1", pop size, and
  #   type (i.e. "H" hermaphroditic)
  #
  # - Mutations:
  #   each mutation gets a separate line with (1) temp id, (2) permanent id,
  #   (3) mutation type,  (4) base position, (5) selection coefficient, 0 for
  #   neutral model, (6) dominance coefficient (0.5 for neutral), (7) id of
  #   sub-population in which the mutation arose, (8) the generation in which
  #   the mutation arose, and (9) the prevalence of the mutation.
  #
  # - Individuals:
  #   lists which genomes belong to which individual, redundant since this
  #   information is also listed at the beginning of each genome.
  #
  # - Genomes:
  #   each line lists one of the genomes for an individual, i.e. individual 1's
  #   genomes are contained in lines 1 & 2, individual 2's genomes are contained
  #   in lines 3 & 4, etc. Each line begins with the genome id, followed by the
  #   type (i.e. "A" for autosome), and a list of mutations.  The mutations are
  #   identified by the temporary ID (tempID) specified in the Mutations output.

  #find heading locations
  PopHead <- which(exDat == "Populations:")
  MutHead <- which(exDat == "Mutations:")
  IndHead <- which(exDat == "Individuals:")
  GenHead <- which(exDat == "Genomes:")

  popCount <- as.numeric(unlist(strsplit(exDat[PopHead + 1], split = " "))[2])

  #-----------#
  # Mutations #
  #-----------#
  #NOTE: time to create Mutations data set ~4 secs
  message("Creating Mutations Dataset")

  #extract mutation data from slim's Mutation output
  #only retaining the tempID, position, and prevalence of each mutation
  MutOut <- do.call(rbind, strsplit(exDat[(MutHead + 1):(IndHead - 1)], split = " ", fixed = TRUE))
  MutData <- data.frame(tempID = as.numeric(MutOut[, 1]),
                        type = MutOut[, 3],
                        position = as.numeric(MutOut[, 4]),
                        #selCoef = as.numeric(MutOut[, 5]),
                        #domCoef = as.numeric(MutOut[, 6]),
                        #pop = MutOut[, 7],
                        #genNo = as.numeric(MutOut[, 8]),
                        prevalence = as.numeric(MutOut[, 9]),
                        stringsAsFactors = TRUE)

  #add 1 to temp ID so that we can easily associate mutations
  #to columns by default slim's first tempID is 0, not 1.
  MutData$tempID <- MutData$tempID + 1
  #First position in slim is 0, not 1
  MutData$position <- MutData$position + 1

  #calculate the derived allele frequency
  #If recode option is TRUE we calculate the derived
  #allele frequecy by first summing the prevelance
  #for the identical mutations and then dividing
  #by the population size
  if (recode_recurrent) {
    calc_afreq <- MutData %>%
      group_by(.data$type, .data$position) %>%
      mutate(afreq = sum(.data$prevalence)/(2*popCount))
    MutData$afreq <- calc_afreq$afreq
  } else {
    MutData$afreq <- MutData$prevalence/(2*popCount)
  }

  #order Mutation dataset by tempID, so that (later) we can order
  #the mutations in each haplotypes by increasing genomic position
  MutData <- MutData[order(MutData$tempID), ]

  #Identify the future sparseMatrix column ID of variants with <= keep_maf.
  #Variants with large maf are assigned column ID 0, i.e. thrown away.
  #Variants with sufficiently small maf are assigned increasing colIDs.
  #These are like new tempIDs, they do not reflect genomic position.
  MutData$colID <- cumsum(MutData$afreq <= keep_maf)*(MutData$afreq <= keep_maf)

  #Using the identified colID, create dataframe of rare mutations only
  RareMutData <- MutData[MutData$colID > 0, ]

  #-----------#
  # Genotypes #
  #-----------#
  #NOTE: jpos is the most computationally expensive task (uses strsplit)
  #this chuck takes ~2.2 minutes to run (for genome-wide, exon-only data)
  #This is an improvement from the old time: ~ 6 mins
  message("Creating Sparse Haplotypes Matrix")

  #determine future row and column position of each mutation listed in genomes
  #row will correspond to person, column will correspond to the tempID of the
  #mutation
  jpos <- lapply(1:(2*popCount), function(x){
    extract_tempIDs(mutString = exDat[GenHead + x],
                    rarePos = MutData$colID)
  })

  ipos <- lapply(1:length(jpos), function(x){
    rep(x, length(jpos[[x]]))
  })

  #create sparse matrix containing mutations(columns) for each individual(row)
  GenoData <- sparseMatrix(i = unlist(ipos),
                           j = unlist(jpos),
                           x = rep(1, length(unlist(jpos))),
                           dims = c(2*popCount, nrow(RareMutData)))

  #order by genomic postion of rare mutation
  GenoData <- GenoData[, order(RareMutData$position)]
  RareMutData <- RareMutData[order(RareMutData$position),]

  #Re-format tempID so that it corresponds to the column
  # (in GenoData) that the mutation is stored in
  RareMutData$colID <- 1:nrow(RareMutData)
  RareMutData <- RareMutData[, -1] #remove the old tempID

  #-----------------------------#
  # Re-code Identical Mutations #
  #-----------------------------#
  # Assuming that all mutations are of the same type, different mutations at
  # the same site are actually identical mutations from different lineages.
  # For simplicity, we recode these mutations so that they are only cataloged
  # once.
  if (recode_recurrent) {
    message("Recoding Identical Mutations")

    #For each mutation type, check to see if there
    #are any identical mutations to re-code
    if (any(sapply(unique(RareMutData$type),
                   function(x){any(duplicated(RareMutData$position[RareMutData$type == x]))}))) {

      #determine which mutation types have identical mutations that need
      #to be recoded
      rc_types <- unique(RareMutData$type)[sapply(unique(RareMutData$type),
                             function(x){any(duplicated(RareMutData$position[RareMutData$type == x]))})]

      #we have to use a loop to recode identical mutations by type
      #because recoding involves combining and removing columns from the
      #Haplotype matrix and then recoding the colIDs in the Mutations data frame.
      for (i in 1:length(rc_types)) {
        com_id_muts <- combine_identicalmutations(mutmap = RareMutData,
                                                  hapmat = GenoData,
                                                  mut_type = rc_types[i])
        RareMutData <- com_id_muts[[1]]
        GenoData <- com_id_muts[[2]]
      }
    }
  }




  #------------------#
  # Re-map Mutations #
  #------------------#
  if (!is.null(recomb_map)) {
    message("Remapping Mutations")
    RareMutData <- reMap_mutations(mutationDF = RareMutData,
                                   recomb_map)
  } else {
    RareMutData$chrom <- 1
  }

  row.names(RareMutData) = NULL

  #Create unique marker names
  RareMutData$marker <- make.unique(paste0(RareMutData$chrom, sep = "_", RareMutData$position))

  #reduce RareMutData, to the columns we actually need
  #really should clean this up soon
  RareMutData <- RareMutData[, c("colID", "chrom", "position",
                                 "afreq", "marker", "type")]

  # RareMutData <- RareMutData[, c("colID", "chrom", "position",
  #                                "afreq", "marker", "type",
  #                                "selCoef", "domCoef", "pop")]

  #----------------------#
  # Identify Pathway RVs #
  #----------------------#
  #if pathway data has been supplied, identify pathway SNVs
  if (!is.null(pathway_df)) {
    message("Identifying Pathway SNVs")
    RareMutData <- identify_pathwaySNVs(markerDF = RareMutData,
                                        pathwayDF = pathway_df)
  }


  return(SNVdata(Haplotypes = GenoData, Mutations = RareMutData))
}



#' Determine i and j positions of mutations for sparse matrix
#'
#' @param mutString character. String containing mutations
#' @param rarePos numeric vector. position of variations
#'
#' @return data.frame with x and y positions of mutations
#' @keywords internal
extract_tempIDs <- function(mutString, rarePos){
  tids <- as.numeric(strsplit(mutString, split = " ", fixed = TRUE)[[1]][-c(1:2)]) + 1
  #Subset colID by tids (tempID) and retain the colIDs that are non-zero,
  #i.e. the rare varaints
  rarePos[tids][rarePos[tids] > 0]
}

#' Combine identical mutations
#'
#' Assuming that all mutations are of the same type (in SLIM simulation), different mutations at the same site are actually identical mutations from different lineages. This function re-codes these mutations so that they are only cataloged once.
#'
#' @inheritParams read_slim
#' @param mutmap data.frame. The SNV_map with identical mutations
#' @param hapmat sparseMatrix. The sparseMatrix of haplotypes
#' @param mut_type character. the name of the mutation type for which to re-code identical mutations
#' @importFrom Matrix rowSums
#'
#' @return  A list containing:
#' @return \item{\code{hapmat}}{A sparse matrix of haplotypes. See details.}
#' @return \item{\code{mutmap}}{A data frame cataloging SNVs in \code{hapmap}.}
#'
#' @keywords internal
combine_identicalmutations <- function(mutmap, hapmat, mut_type){

  #store the information for this mutation type
  mut_dat <- mutmap[mutmap$type == mut_type, ]

  # identify the positions at which identical mutations
  # from different lineages exist for this mutation type.
  im_pos = unique(mut_dat$position[duplicated(mut_dat$position)])

  # find the column locations for the identical mutations
  # for this mutation type
  col_loc <- lapply(im_pos, function(x){
    mutmap$colID[which(mutmap$position == x & mutmap$type == mut_type)]
  })

  #find the combined SNV data
  comb_mut <- lapply(col_loc, function(x){
    rowSums(hapmat[, x])
  })

  keep_SNVcol <- sapply(col_loc, function(x){x[1]})

  #replace with combined haplotype data
  hapmat[, keep_SNVcol] <- do.call(cbind, comb_mut)

  #determine the superfluous columns and remove them, since we have
  #already accounted for them by combining the columns for this SNV
  #NOTE: Do NOT use sapply, will not wirk if a SNV has more than 1 replicate
  remove_cols <- unlist(lapply(col_loc, function(x){
    x[-1]
  }))

  #remove the columns of the identical SNVs
  hapmat <- hapmat[, -remove_cols]

  #remove the corresponding rows from mutmap and re-lable the colID
  mutmap <- mutmap[-remove_cols, ]
  mutmap$colID <- 1:nrow(mutmap)

  return(list(mutmap, hapmat))
}
