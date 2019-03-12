#' Identify SNVs located in a specified pathway
#'
#' @param markerDF A data frame containing SNV data.   See details.
#' @param pathwayDF A data frame containing pathway data.  See details.
#'
#' @return A dataframe with the same format as \code{markerDF}, but with a variable \code{pathwaySNV} marked FALSE for all variants located outside exons in \code{pathwayDF}, or with derived allele frequency greater than \code{carrier_prob}.
#' @keywords internal
identify_pathwaySNVs <- function(markerDF, pathwayDF){

  if (is.null(markerDF$pathwaySNV)) {
    markerDF$pathwaySNV <- TRUE
  }

  #Mark pathwaySNV FALSE for any variants that do not fall
  #within the exons catalouged in pathwayDF
  markerDF <- do.call(rbind, lapply(unique(markerDF$chrom), function(x){
    identify_pathwayRVs_byChrom(pathwayDF[pathwayDF$chrom == x, ],
                        markerDF[markerDF$chrom == x, ])
  }))

  return(markerDF)

}


#' Identify variants located in a defined pathway
#'
#' @param path_by_chrom The pathway data for the chromosome under consideration
#' @param marker_map_by_chrom The marker_map for the chromosome under consideration
#'
#' @return marker_map_by_chrom with pathwaySNV identified
#' @keywords internal
identify_pathwayRVs_byChrom <- function(path_by_chrom, marker_map_by_chrom){
  if(nrow(path_by_chrom) == 0){
    marker_map_by_chrom$pathwaySNV <- FALSE
  } else {
    #since we will want to include variants that occur at the first base pair
    #location, subtracting 1 from start positions
    #Similarly, since we want the first mutation in marker_map_by_chrom
    #to have a bin, subtracting 1 from this position
    if(min(marker_map_by_chrom$position) != min(path_by_chrom$exonStart)){
      cbreaks <- sort(c((min(marker_map_by_chrom$position) - 1),
                        (path_by_chrom$exonStart - 1),
                        unique(c(path_by_chrom$exonEnd,
                                 max(marker_map_by_chrom$position)))))

      keep_bins <- seq(2, length(cbreaks) - 1, by = 2)
    } else {
      cbreaks = sort(c((path_by_chrom$exonStart - 1),
                       unique(c(path_by_chrom$exonEnd,
                                max(marker_map_by_chrom$position)))))
      keep_bins <- seq(1, length(cbreaks) - 1, by = 2)
    }

    if(any(duplicated(cbreaks))){
      stop("Expecting non-overlapping exonStart and exonEnd positions. \n Please combine overlapping segments into a single entry.")}

    marker_map_by_chrom$pathwaySNV <- cut(marker_map_by_chrom$position,
                                         breaks = cbreaks,
                                         labels = FALSE)
    marker_map_by_chrom$pathwaySNV <- marker_map_by_chrom$pathwaySNV %in% keep_bins
  }

  return(marker_map_by_chrom)
}
