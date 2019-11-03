#' Create chromosome map from marker map
#'
#' INTENDED FOR INTERNAL USE ONLY
#'
#' @param SNV_map Data frame. A data frame that catalogs the SNVs in \code{haplos}.  If the \code{\link{read_slim}} function was used to import SLiM data to \code{R}, the data frame \code{Mutations} is of the proper format for \code{SNV_map}.  However, users must add the variable \code{is_CRV} to this data frame, see details.

#'
#' @return a dataframe catalouging the start and stop positions, in base pairs, for each chromosome.  We use this information to determine what regions to simulate recombination over.
#' @keywords internal
create_chrom_map <- function(SNV_map){
  chrom_map <- do.call(rbind, lapply(sort(unique(SNV_map$chrom)), function(x){
    c(x, range(SNV_map$position[SNV_map$chrom == x]))
  }))

  chrom_map <- as.data.frame(chrom_map)
  colnames(chrom_map) = c("chrom", "start_pos", "end_pos")
  return(chrom_map)
}
