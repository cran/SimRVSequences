#' Construct offspring sequence from parental allele vector
#'
#' For internal use.
#'
#' Prior to running this function, for each chromosome, we have a coded sequence (i.e. c(1, 2, 1)) and a list of crossover points (i.e. c(45, 100)). This data would indicate that the offspring inherts all genetic data from the start postion to postion 45 of the parent's 1st haplotype, then all genetic data after position 45 up to postion 100 from the 2nd haplotype, and all genetic data after postion 100 from the 1st haplotype.  This function completes that task.
#'
#' @param parental_genotypes The parental genotype sequence information.
#' @param CSNV_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosome that the marker resides on, column 3: the position of the marker in cM.
#' @param inherited_haplotype The coded haplotype sequence, which indiciates the haplotype that the data will come from.
#' @param chiasmata_locations A list of crossover locations.
#' @param REDchrom_map Data.frame.  The chromosome map, reduced to the chromosome in question.
#'
#' @importFrom Matrix colSums
#'
#' @return offspring_seq. The genetic data that the offspring inherits from this parent.  This will be a recombined sequence.
#' @keywords internal
reconstruct_fromHaplotype <- function(parental_genotypes,
                                      CSNV_map,
                                      inherited_haplotype,
                                      chiasmata_locations,
                                      REDchrom_map){

  # reduce chiasmata_locations to the chiasmata that
  # the haplotypes participated in
  cross_loc <- reduce_to_events(as.numeric(inherited_haplotype), chiasmata_locations)

  if (class(parental_genotypes) == "numeric"){
    #If parental haplotype is of length 1, we return the inherited haplotype
    offspring_seq <- parental_genotypes[as.numeric(inherited_haplotype[1])]

  } else if (all(colSums(parental_genotypes) == 0) | length(cross_loc) == 0) {
    #if parental haplotypes types do not contain any SNVs (for this chrom)
    #we can return either of the two original parental haplotypes
    #as they will remain unchanged after recombination
    #
    #If parental haplotypes do not participate in recombination events
    #we return the inherited haplotype
    offspring_seq <- parental_genotypes[as.numeric(inherited_haplotype[1]), ]

  } else if (length(cross_loc) > 0){
    #This is the non-trival case. i.e. more than 1 SNV on the chrom
    #with crossovers

    #determine which inherited haplotype participated in chiasmata
    switch_alleles_loc <- c(REDchrom_map$start_pos,
                            cross_loc,
                            REDchrom_map$end_pos + 1) #1 added here in case of marker at end of chromosome

    #determine first allele ID in inherited_haplotype
    #set the offspring's haplotype sequence to the other (original) parental haplotype
    offspring_seq <- parental_genotypes[ifelse(inherited_haplotype[1, 1] == 1, 2, 1), ]

    #store the first allele ID in the haplotype sequence
    switch_alle <- inherited_haplotype[1, 1]

    #cycle through all swaps
    for(i in 1:(length(switch_alleles_loc) %/% 2)){
      start_switch <- switch_alleles_loc[2*i - 1]
      end_switch <- switch_alleles_loc[2*i]

      #determine the columns that need to be swapped
      swap_cols <- which(CSNV_map$position >= start_switch & CSNV_map$position < end_switch)

      #Occasionally, there is no marker data to swap since markers may be far apart
      #when this occurs length(swap_cols) = 0. When this is the case we do nothing.
      if (length(swap_cols) > 0) {
        #switch alleles between crossovers
        offspring_seq[swap_cols] <- parental_genotypes[switch_alle, swap_cols]
      }
    }

  }

  return(offspring_seq)
}

#' Simulate sequence data for a pedigree
#'
#' @inheritParams sim_gameteInheritance
#' @inheritParams sim_RVstudy
#' @param ped_file Data frame. A single pedigree. Must match format of pedigree simulated by sim_RVped
#' @param RV_marker character. The marker name of the RV locus.
#' @param SNV_map Data frame. A data frame that catalogs the SNVs in \code{haplos}.  If the \code{\link{read_slim}} function was used to import SLiM data to \code{R}, the data frame \code{Mutations} is of the proper format for \code{SNV_map}.  However, users must add the variable \code{is_CRV} to this data frame, see details.
#' @param founder_genos Dataframe.  A dataframe with rows corresponding to founders, and columns corresponding to markers.  Markers must be listed in same order as \code{SNV_map}.
#'
#' @return offspring_sequences
#' @keywords internal
sim_seq <- function(ped_file, founder_genos,
                    SNV_map, chrom_map, RV_marker,
                    burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)){

  #Get parent/offspring information
  #i.e. for each offspring find RV_status,
  #parent IDs, and parent alleles at RV locus
  PO_info <- get_parOffInfo(ped_file)
  PO_info <- PO_info[order(PO_info$Gen, PO_info$offspring_ID),]

  ped_genos <- founder_genos[[1]]
  ped_geno_IDs <- founder_genos[[2]]

  #determine the chromosome number and location of the familial RV locus
  #then store as a dataframe with chrom in the first column
  RVL <- SNV_map[which(SNV_map$marker == RV_marker),
                 which(colnames(SNV_map) %in% c("chrom", "position"))]

  if(colnames(RVL[1]) != "chrom"){
    RVL <- RVL[, c(2, 1)]
  }

  #for each offspring simulate transmission of parental data
  for (i in 1:nrow(PO_info)) {
    #simulate recombination events for this parent offspring pair
    loop_gams <- sim_gameteInheritance(RV_locus = RVL,
                                       parent_RValleles = PO_info[i, c(6, 7)],
                                       offspring_RVstatus = PO_info[i, 5],
                                       chrom_map,
                                       allele_IDs = c(1, 2),
                                       burn_in, gamma_params)

    #construct offspring's inherited material from this parent
    loop_seq <- lapply(c(1:nrow(chrom_map)),
                       function(x){
                         reconstruct_fromHaplotype(parental_genotypes =
                                                     ped_genos[which(ped_geno_IDs == PO_info[i, 4]),
                                                                        which(SNV_map$chrom == chrom_map$chrom[x])],
                                                   CSNV_map = SNV_map[which(SNV_map$chrom == chrom_map$chrom[x]),],
                                                   inherited_haplotype = loop_gams$haplotypes[[x]],
                                                   chiasmata_locations = loop_gams$cross_locations[[x]],
                                                   REDchrom_map = chrom_map[x, ])
                       })

    #append ID for this haplotype to the list of IDs
    ped_geno_IDs <- c(ped_geno_IDs, PO_info[i, 1])

    #TODO: preallocate rows in ped_genos based on number of non-founders
    #append this haplotype to the other familial haplotypes
    ped_genos <- rbind(ped_genos, unlist(loop_seq))
  }

  #Determine if this is a sporadic pedigree
  printed_FamRV <- ifelse(all(ped_file[, c("DA1", "DA2")] == 0), "no_CRV", RV_marker)

  #create a data.frame to store identifying info
  geno_map <- data.frame(FamID = rep(ped_file$FamID[1], length(ped_geno_IDs)),
                         ID = ped_geno_IDs,
                         affected =  rep(FALSE, length(ped_geno_IDs)),
                         FamCRV = rep(printed_FamRV, length(ped_geno_IDs)))

  #identify affected individuals
  geno_map$affected[geno_map$ID %in% ped_file$ID[ped_file$affected]] <- TRUE

  #Return the genomes matrix and a data.frame continating identifying
  #information for the of IDs to identify the
  #family member to whom
  return(list(ped_genos = ped_genos, geno_map = geno_map))
}
