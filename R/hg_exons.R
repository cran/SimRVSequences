#' Human exon data
#'
#' This data set catalogs combined exon segments from the 22 human autosomes.
#'
#' The \code{hg_exons} data set catalogs the positions of exons residing in the 22 human autosomes. The data contained in \code{hg_exons} was collected from the hg 38 reference genome with the UCSC Genome Browser's Table Brower Tool.  In \code{hg_exons} overlapping exons have been combined into a single observation. When exons from genes with different NCBI accession numbers have been combined the variable \code{NCBIref} will contain multiple accession numbers separated by commas.  We note that different accession numbers may exist for transcript variants of the same gene.
#'
#'
#' @docType data
#'
#' @references Karolchik, D., Hinrichs, A. S., Furey, T. S., Roskin, K. M., Sugnet, C. W., Haussler, D., and Ken, W. J. (2004). The UCSC Table Browser data retrieval tool. \emph{Nucleic Acids Res}. Accessed on 6 February 2018.
#' @references Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, T. H., Zahler, A. M., and Haussler, D. (2002). The human genome browser at UCSC. \emph{Genome Res}, 12(6):996-1006.
#'
#' @format A data set with 223565 rows and 4 variables:
#' \describe{
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{exonStart}{Numeric. The position of the first base pair in the combined exon segment.}
#'   \item{exonStop}{Numeric. The position of the last base pair in the combined exon segment.}
#'   \item{NCBIref}{Character. The NCBI reference sequence accession number of the gene(s) in which the exon(s) reside.}
#' }
"hg_exons"
