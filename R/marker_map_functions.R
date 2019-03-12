#' Create chromosome map from marker map
#'
#' INTENDED FOR INTERNAL USE ONLY
#'
#' @inheritParams sim_RVstudy
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
