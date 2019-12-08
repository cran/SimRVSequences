#' Get parent and offspring information from a pedigree
#'
#' \strong{For internal use.}
#'
#' @param ped_file Data.frame. The pedigree file, must have same format as pedigree simulated with \code{sim_RVped}
#'
#' @return A list containing the parent's paternal and maternal alleles at the disease locus, and the RV status of the offspring
#' @keywords internal
#' @importFrom reshape2 melt
get_parOffInfo <- function(ped_file){

  mdata <- melt(ped_file[which(!is.na(ped_file$dadID)),
                         which(colnames(ped_file) %in%
                                 c("ID", "dadID", "momID", "Gen"))],
                id = c("ID", "Gen"))
  colnames(mdata) = c("offspring_ID", "Gen", "parent", "parent_ID")


  # mdata$Off_RVstatus <- apply(as.data.frame(mdata[, 1]), 1, function(x){
  #   sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) %in% c("DA1", "DA2"))])
  # })

  mdata$Off_RVstatus <- sapply(1:nrow(mdata), function(x){
    ifelse(mdata$parent[x] == "dadID",
           ped_file$DA1[ped_file$ID == mdata$offspring_ID[x]],
           ped_file$DA2[ped_file$ID == mdata$offspring_ID[x]])})


  mdata$Par_DA1 <- sapply(mdata[, 4], function(x){
    sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) == "DA1")])
  })

  # mdata$Par_DA2 <- apply(as.data.frame(mdata[, 4]), 1, function(x){
  #   sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) == "DA2")])
  # })

  mdata$Par_DA2 <- sapply(mdata[, 4], function(x){
    sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) == "DA2")])
  })

  return(mdata)
}


#' Reduce chiasmata vector to crossovers that this gamete participated in based on the allele vector.
#'
#' For internal use.
#'
#' For example, if coded allele vector was c(1, 1, 2), then this gamete did not participate in the first crossover,
#' hence the list of chiasmata event locations, say c(40, 80), would be reduced to c(80), and the coded allele vector
#' would be reduced to c(1, 2). That is, there is only one crossover from haplotype 1 to 2 at position 80.
#'
#' We don't actually return the coded vector at this point.  We only need the code for the first haplotype, after applying this function every crossover will be meaningful. That is the coded vector will never contain repeats after this function runs.
#'
#' @param gamete_haplo Numeric vector. The inherited haplotype.
#' @param chias_locations  Numeric vector.  Chiasmata locations.
#'
#' @return The locations of crossovers
#' @keywords internal
reduce_to_events <- function(gamete_haplo, chias_locations){
  if(sum(gamete_haplo == gamete_haplo[1]) == length(gamete_haplo)){
    cross_loc = numeric(0)
  } else {
    keep_ind <- c()
    i = 1

    #find position of last match to first element of gamete_haplo
    keep_ind[i] <- Position(function(x){x != gamete_haplo[1]}, gamete_haplo) - 1
    d = keep_ind[i]

    while (d < length(gamete_haplo)) {
      rhap <- gamete_haplo[c((keep_ind[i] + 1) : length(gamete_haplo))]
      if ( sum(rhap == rhap[1]) == length(rhap) | length(rhap) == 1 ) {
        #if the remaining items to check have length 1 or are
        #all the same element we are done; i.e. break while condition
        d <- length(gamete_haplo) + 1
      } else {
        keep_ind[i + 1] <- keep_ind[i] +
          (Position(function(x){x != rhap[1]}, rhap) - 1)

        d <- keep_ind[i + 1]
        i = i + 1
      }
    }
    cross_loc <- chias_locations[keep_ind]
  }
  return(cross_loc)
}

#' Determine if input is an odd number
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @keywords internal
is_odd <- function(x) {x %% 2 != 0}


#' Determine if input is an integer
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @keywords internal
is_int <- function(x) {x %% 1 == 0}

#' Remove unaffected relatives
#'
#' Remove unaffected relatives
#'
#' @param ped_file data.frame. A pedigree.
#'
#' @return \code{retA_ped} A pedigree containing only affected members, obligate carriers, and founders.
#' @keywords internal
affected_onlyPed = function(ped_file){

  #assign individuals with unknown affection status to FALSE,
  #since we will not simulate sequence data for these
  #individuals unless necessary

  if (any(is.na(ped_file$affected))) {
    ped_file$affected[which(is.na(ped_file$affected))] = FALSE
  }

  #create new ped file with affecteds only
  retA_ped <- ped_file[ped_file$affected, ]

  if (nrow(retA_ped) == 0) {
    #warning(paste0("Removing pedigree with with FamID ", sep = "", ped_file$FamID[1],  ": No disease-affected relatives. \n"))
    return(retA_ped)
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(retA_ped$dadID,
                               retA_ped$ID[which(retA_ped$sex == 0)])
      readd_dad <- retA_ped$dadID[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(retA_ped$momID,
                               retA_ped$ID[which(retA_ped$sex == 1)])
      readd_mom <- retA_ped$momID[miss_mom]
      readd_mom <- unique(readd_mom[!is.na(readd_mom)])

      #check to see if we need to readd anyone
      if (length(c(readd_dad, readd_mom)) == 0) {
        d <- 1
      } else {
        #Now pull the rows containing the required parents
        # from the original ped_file
        readd <- ped_file[which(ped_file$ID %in% c(readd_dad, readd_mom)), ]

        #combine with affected ped file
        retA_ped <- rbind(retA_ped, readd)
      }
    }
  }

  return(retA_ped)
}

#' Convert from basepairs to centimorgan
#'
#' Convert from basepairs to centimorgan
#'
#' @param pos_BP Numeric.  The position in basepairs.
#'
#' @return pos_CM The postion in centiMorgans
#' @keywords internal
convert_BP_to_cM <- function(pos_BP){ pos_BP/1000000 }

#' Convert from centiMorgan to basepairs
#'
#' Convert from centiMorgan to basepairs
#'
#' @param pos_CM Numeric.  The position in centiMorgan.
#'
#' @return pos_BP The postion in basepairs
#' @keywords internal
convert_CM_to_BP <- function(pos_CM){ pos_CM*1000000 }


#' Assign generation number based on oldest founder
#'
#' @param x an object of class ped
#'
#' @return a list of generation numbers for pedigree members, in the order listed in \code{x}.
#' @importFrom kinship2 kindepth
#' @importFrom SimRVPedigree ped2pedigree
#' @importFrom SimRVPedigree new.ped
#' @importFrom methods is
#' @keywords internal
assign_gen <- function(x){
  # create a ped object
  # this will also check to see that
  # all required fields are present
  # i.e. FamID, ID, dadID, momID, sex, and affected
  if (!is(x, "ped")) x <- new.ped(x)

  Gen <- NA
  mates <- cbind(x$dadID, x$momID)
  #remove all rows with only zeros, these are founders
  mates <- unique(mates)
  mates <- mates[!is.na(mates[, 1]), ]

  #get kindepth and set Gen to kindepth when kindepth is non-zero
  kd <- kindepth(ped2pedigree(x))
  Gen[kd != 0] <- kd[kd != 0]

  if(is(mates, "matrix")){
    for(i in 1:nrow(mates)){
      mate_gens <-  kd[x$ID %in% mates[i, ]]
      Gen[x$ID %in% mates[i, ]] <- max(mate_gens)
    }
  } else {
    Gen[x$ID %in% mates] <- 0
  }

  return(Gen + 1)
}
