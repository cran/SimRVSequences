#' Combine overlapping exons
#'
#' Combine overlapping exons into a single observation
#'
#' @param exon_data data.frame. This data frame must include named variables: \code{chrom}, a chromosome identifer; \code{exonStart}, the first position of the exon in base pairs; and \code{exonEnd}, the last position of the exon in base pairs.
#'
#' @return A data frame of combined exon segments.  This data frame includes the variables: \code{chrom}, a chromosome identifier; \code{exonStart}, the first position of the combined exon segment in base pairs; and \code{exonEnd}, the last position of the combined exon segment in base pairs.
#' @export
#'
#' @examples
#' # create a data frame that contains the
#' # the variables: chrom, exonStart, and exonEnd
#' exDat <- data.frame(chrom     = c(1, 1, 1, 2, 2, 2),
#'                     exonStart = c(1, 2, 5, 1, 3, 3),
#'                     exonEnd   = c(3, 4, 7, 4, 5, 6))
#'
#' exDat
#'
#' # supply exDat to combine_exons to combine
#' # overlapping exon segments
#' combine_exons(exDat)
#'
combine_exons <- function(exon_data){

  #check to see if exon_data contains the three named
  #columns we require to combine exons
  if(any(!(c("chrom", "exonStart", "exonEnd") %in% colnames(exon_data)))){
    stop("exon_data does not include named columns: 'chrom', 'exonStart', and 'exonEnd'.")
  }

  #just in case, order
  exon_data <- exon_data[order(exon_data$chrom, exon_data$exonStart, exon_data$exonEnd), ]

  cexons <- do.call(rbind, lapply(unique(exon_data$chrom), function(x){
    combine_exons_by_chrom(chrom = x,
                           start_stop_dat = exon_data[exon_data$chrom == x, c("exonStart", "exonEnd")])
  }))

 return(as.data.frame(cexons))

}


#' Combine exons within a chromosome
#'
#' @param chrom the chromosome number
#' @param start_stop_dat the exon start and stop data, i.e. two columns of a dataframe or matrix.  Start positons should be contained in column 1 and stop positions in column 2.
#'
#' @importFrom intervals interval_union
#' @importFrom intervals Intervals
#'
#' @return a matrix with combined exons for a single chromosome
#' @keywords internal
combine_exons_by_chrom <- function(chrom, start_stop_dat){

  #combine exons in this chromosome using the interval_union
  #and Intervals functions provided by the intervals package
  ex_mat <- interval_union(Intervals(start_stop_dat))@.Data


  #check to see if any exons are directly next to one another,
  #if so shift end point so that they get combined into a single segment
  adjacent_exons <- which((ex_mat[-1, 1] - ex_mat[-nrow(ex_mat), 2]) == 1)
  test_var <- length(adjacent_exons) > 0
  while (test_var) {
    ex_mat[adjacent_exons, 2] <- ex_mat[(adjacent_exons + rep(1, length(adjacent_exons))), 2]
    ex_mat <- interval_union(Intervals(ex_mat))@.Data
    adjacent_exons <- which((ex_mat[-1, 1] - ex_mat[-nrow(ex_mat), 2]) == 1)
    test_var <- length(adjacent_exons) > 0
  }

  #add chromosome variable
  ex_mat <- cbind(rep(chrom, nrow(ex_mat)), ex_mat)

  #add column names
  colnames(ex_mat) <- c('chrom', 'exonStart', 'exonEnd')

  return(ex_mat)
}
